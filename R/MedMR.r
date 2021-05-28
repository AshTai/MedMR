#' Multiply robust estimation of natural indirect effects with multiple ordered mediators.
#'
#' MedMR is used to assess the natural indirect effects with multiple mediators by multiply robust estimation, G-computation, and inverse probability weighting.
#'
#' @name MedMR
#' @author An-Shun Tai \email{anshuntai@nctu.edu.tw} Chen-Chien Lee\email{xu3t6fu0@gmail.com}
#' @param data a data frame, where the column is the variable and row is the sample.
#' @param exposure variable name of the exposure.
#' @param mediators variable names of mediators.
#' @param outcome variable name of the outcome.
#' @param confounder variable names of confounders.
#' @param m_type a description of the mediator model. This should be one of "d" (logistic model), "normal" (Gaussian model), "exp" (exponential model), and "poisson" (poisson model).
#' @param y_type a description of the outcome model. "discrete" uses logistic model, and "continuous" uses Gaussian model.
#' @param method a description of the estimation method. "IPW" stands for the inverse probability weighting; "G" stands for G-computation; "Robust" stands for multiply robust estimation. Default is "Robust".
#' @param boot a logical value. if 'TRUE' nonparametric bootstrap will be used for SD estimation. Default is 'FALSE'.
#' @param boot_ratio a ratio of the random samples of individuals from the original data.
#' @param boot_rep a number of replicates
#' @param seed_num seed. Default is 123.
#' @param double_mont a number of Monte Carlo simulations in multiply robust estimation.
#' @param single_mont a number of Monte Carlo simulations in G-computation.
#' @param b_rep if true resampling with replacement will be used; if false resampling without replacement will be used.
#' @return A list of estimated direct, indirect, and total effects.
#' TE is the total effect. Seq0 is the natural direct effect.
#' Seq1,...SeqK are the natural indirect effects of mediators.
#' @export
#' @examples
#' data("HCCexample")
#' ff <- MedMR(data=HCCexample,exposure="a",mediators=c("m1","m2"),outcome="y",
#' confounder=c("gender","age","smoke","alcohol"),
#' m_type=c("d","d"),y_type = "discrete",
#' method = "Robust",
#' scale = "difference",
#' boot = T,
#' boot_ratio = 0.8,
#' boot_rep = 20,
#' seed_num = 123,
#' double_mont = 1e5,
#' single_mont = 2e5,
#' b_rep = TRUE)


MedMR = function (data,exposure,mediators,outcome,confounder,
                             m_type,y_type = "discrete",
                             method = "Robust",
                             scale = "difference",
                             boot = FALSE,
                             boot_ratio = 0.8,
                             boot_rep = 10,
                             seed_num = 123,
                             double_mont = 1e5,
                             single_mont = 2e5,
                             b_rep = TRUE)
{
  require(progress)
  set.seed(seed_num)
  n_sample = dim(data)[1]
  n_effect = length(mediators) +1
  mrPF = function (data_in)
  {
    k = length(mediators)
    all_status = list()
    all_status[[1]] = rep(0,(k+1))

    for(i in c(2:(k+2)))
    {
      all_status[[i]] = c(rep(1,i-1),rep(0,k+2-i))
    }

    n_sample = dim(data_in)[1]

    # Sampling function
    model_sample = function(size,parameter,type = "normal")
    {
      out = switch (type,
                    normal = rnorm(size,parameter[,1],parameter[,2]),
                    poisson = rpois(size,parameter[,1]),
                    exp = rexp(size,(1/parameter[,1])))
      return(out)
    }

    # Prob function
    model_prob = function(samp,parameter,type = "normal")
    {
      out = switch (type,
                    normal = dnorm(samp,parameter[,1],parameter[,2]),
                    poisson = dpois(samp,parameter[,1]),
                    exp = dexp(samp,(1/parameter[,1])))
      return(out)
    }

    # MODEL M

    M_model = list()
    M_sd = numeric(k)
    fmla = paste(c(mediators[1],paste(c(exposure,confounder),collapse="+")),collapse ="~")
    if (m_type[1] != "d")
    {
      model1 = lm(fmla,data = data_in)
      M_model[[1]] = model1
      M_sd[1] = summary(model1)$sigma
    }else
    {
      model1 = glm(fmla,data = data_in,family = binomial(link = "logit"))
      M_model[[1]] = model1
      M_sd[1] = NA
    }

    for (i in c(2:k))
    {
      fmla = paste(c(mediators[i],paste(c(exposure,mediators[1:(i-1)],confounder),collapse="+")),collapse ="~")
      if (m_type[i] != "d")
      {
        model1 = lm(fmla,data = data_in)
        M_model[[i]]=model1
        M_sd[i] = summary(model1)$sigma
      }else
      {
        model1 = glm(fmla,data = data_in,family = binomial(link = "logit"))
        M_model[[i]] = model1
        M_sd[i] = NA
      }
    }


    # MODEL Y

    Y_model = list()
    fmla = paste(c(outcome,paste(c(exposure,mediators,confounder),collapse="+")),collapse ="~")

    if (y_type == "discrete")
    {
      Y_model = glm(fmla,data = data_in,family = binomial(link = "logit"))
      Y_sd = summary(Y_model)$sigma
    }else
    {
      Y_model = lm(fmla,data = data_in)
      Y_sd = summary(Y_model)$sigma
    }

    # MODEL A

    A_model = list()

    for (i in c(1:(k-1)))
    {
      fmla = paste(c(exposure,paste(c(mediators[1:i],confounder),collapse="+")),collapse ="~")
      model1 = glm(fmla,data = data_in,family = binomial(link = "logit"))
      A_model[[i]]=model1
    }

    fmla = paste(c(exposure,paste(confounder,collapse="+")),collapse ="~")
    A_model_pure = glm(fmla,data = data_in,family = binomial(link = "logit"))

    # IPW part
    IPW = function (status)
    {
      a = status[1]
      e = status[-1]
      up = 1
      for (i in c(1:k))
      {
        if(m_type[i] != "d")
        {
          data_pre = data_in
          data_pre[exposure] = e[i]
          mean1 = predict(M_model[[i]],data_pre)
          parameter1 = cbind(mean1,M_sd[i])
          mi = as.matrix(data_in[mediators[i]])
          up = model_prob(mi,parameter1,type = m_type[i]) * up
        }else
        {
          data_pre = data_in
          data_pre[exposure] = e[i]
          mean1 = predict.glm(M_model[[i]],data_pre,type = "response")
          mean1_c = 1 - mean1
          mi = ifelse(data_in[mediators[i]]==1,mean1,mean1_c)
          up = up * mi
        }

      }

      down = 1

      for (i in c(1:k))
      {
        if(m_type[i] != "d")
        {
          data_pre = data_in
          data_pre[exposure] = a
          mean1 = predict(M_model[[i]],data_pre)
          parameter1 = cbind(mean1,M_sd[i])
          mi = as.matrix(data_in[mediators[i]])
          down = model_prob(mi,parameter1,type = m_type[i]) * down
        }else
        {
          data_pre = data_in
          data_pre[exposure] = a
          mean1 = predict.glm(M_model[[i]],data_pre,type = "response")
          mean1_c = 1 - mean1
          mi = ifelse(data_in[mediators[i]]==1,mean1,mean1_c)
          down = down * mi
        }
      }

      ac_prob = predict.glm(A_model_pure,data_in,type = "response")
      if (a ==0)
      {
        ac_prob = 1 - ac_prob
      }
      down = down * ac_prob

      Ind = ifelse(data_in[exposure]==a,1,0)

      out = (Ind * up) / down

      return (out)
    }

    # Robust + Reg part
    Reg_Robust = function (status)
    {
      set.seed(seed_num)
      c_sample_pos = sample(c(1:n_sample),double_mont,replace = TRUE)
      c_sample = data_in[confounder][c_sample_pos,]
      a_sample = data_in[exposure][c_sample_pos,]
      # status_in = status
      inner_sample = function (c_samp,status_in)
      {
        a = status_in[1]
        e = status_in[-1]
        v = matrix(rep(0,k),nrow = double_mont,ncol = k)
        colnames(v) = mediators
        v = as.data.frame(v)

        for (i in c(1:k))
        {
          if (m_type[i] != "d")
          {
            data_pre = cbind(e[i],v,c_samp,row.names = NULL)
            names(data_pre) = paste(c(exposure,mediators,confounder))
            mean1 = predict(M_model[[i]],data_pre)
            sd1 = M_sd[i]
            parameter1 = cbind(mean1,sd1)
            v1 = model_sample(double_mont,parameter1,type = m_type[i])
            v[mediators[i]] = v1
          }else
          {
            data_pre = cbind(e[i],v,c_samp,row.names = NULL)
            names(data_pre) = paste(c(exposure,mediators,confounder))
            mean1 = predict.glm(M_model[[i]],data_pre,type = "response")
            v1 = runif(double_mont,0,1)
            v1 = ifelse(v1 < mean1,1,0)
            v[mediators[i]] = v1
          }

        }

        data_sample = cbind(a,v,c_samp,row.names = NULL)
        names(data_sample) = paste(c(exposure,mediators,confounder))

        if (y_type == "discrete")
        {
          EY = predict.glm(Y_model,data_sample,type = "response")
        }else
        {
          EY = predict(Y_model,data_sample)
        }


        return (EY)
      }
      outer_sample = function(a_samp,c_samp,status_in)
      {
        a = status_in[1]
        e = status_in[-1]
        v = matrix(rep(0,k),nrow = double_mont,ncol = k)
        colnames(v) = mediators
        v = as.data.frame(v)

        for (i in c(1:k))
        {
          if (m_type[i] != "d")
          {
            data_pre = cbind(a_samp,v,c_samp,row.names = NULL)
            names(data_pre) = paste(c(exposure,mediators,confounder))
            mean1 = predict(M_model[[i]],data_pre)
            sd1 = M_sd[i]
            parameter1 = cbind(mean1,sd1)
            v1 = model_sample(double_mont,parameter1,type = m_type[i])
            v[mediators[i]] = v1
          }else
          {
            data_pre = cbind(a_samp,v,c_samp,row.names = NULL)
            names(data_pre) = paste(c(exposure,mediators,confounder))
            mean1 = predict.glm(M_model[[i]],data_pre,type = "response")
            v1 = runif(double_mont,0,1)
            v1 = ifelse(v1 < mean1,1,0)
            v[mediators[i]] = v1
          }

        }

        data_sample = cbind(a_samp,v,c_samp,row.names = NULL)
        names(data_sample) = paste(c(exposure,mediators,confounder))

        up = 1
        for (i in c(1:(k-1)))
        {
          aprob = predict.glm(A_model[[i]],data_sample,type = "response")
          if (e[i] != 1)
          {
            aprob = 1-aprob
          }
          up = up * aprob
        }

        down = predict.glm(A_model_pure,data_sample,type = "response")
        if (e[1] != 1)
        {
          down = 1 - down
        }

        for (i in c(1:(k-1)))
        {
          aprob = predict.glm(A_model[[i]],data_sample,type = "response")

          if (e[(i+1)] != 1)
          {
            aprob = 1-aprob
          }

          down = down * aprob

        }
        Ind = ifelse(data_sample[exposure] == e[k],1,0)
        return( (Ind * up) / down )

      }
      in_res = inner_sample(c_sample,status)
      out_res = outer_sample(a_sample,c_sample,status)
      comb_res = in_res * out_res
      return(mean(comb_res))
    }

    # Pure Reg

    P_Reg = function(status)
    {
      set.seed(seed_num)
      a = status[1]
      e = status[-1]
      pos_samp = sample(c(1:n_sample),single_mont,replace = TRUE)
      v = matrix(rep(0,k),nrow = single_mont,ncol = k)
      colnames(v) = mediators
      v = as.data.frame(v)
      c_row = data_in[confounder][pos_samp,]

      for (i in c(1:k))
      {
        if (m_type[i] != "d")
        {
          data_pre = cbind(e[i],v,c_row,row.names = NULL)
          names(data_pre) = paste(c(exposure,mediators,confounder))
          mean1 = predict(M_model[[i]],data_pre)
          sd1 = M_sd[i]
          parameter1 = cbind(mean1,sd1)
          v1 = model_sample(single_mont,parameter1,type = m_type[i])
          v[mediators[i]] = v1
        }else
        {
          data_pre = cbind(e[i],v,c_row,row.names = NULL)
          names(data_pre) = paste(c(exposure,mediators,confounder))
          mean1 = predict.glm(M_model[[i]],data_pre,type = "response")
          v1 = runif(single_mont,0,1)
          v1 = ifelse(v1 < mean1,1,0)
          v[mediators[i]] = v1
        }
      }

      data_sample = cbind(a,v,c_row,row.names = NULL)
      names(data_sample) = paste(c(exposure,mediators,confounder))

      if (y_type == "discrete")
      {
        EY = predict.glm(Y_model,data_sample,type = "response")
      }else
      {
        EY = predict(Y_model,data_sample)
      }
      return (mean(EY))
    }

    RBW = function(status)
    {
      a = status[1]
      e = status[-1]

      up = 1
      for (i in c(1:(k-1)))
      {
        aprob = predict.glm(A_model[[i]],data_in,type = "response")
        if (e[i] != 1)
        {
          aprob = 1-aprob
        }
        up = up * aprob
      }

      down = predict.glm(A_model_pure,data_in,type = "response")
      if (e[1] != 1)
      {
        down = 1 - down
      }

      for (i in c(1:(k-1)))
      {
        aprob = predict.glm(A_model[[i]],data_in,type = "response")

        if (e[(i+1)] != 1)
        {
          aprob = 1-aprob
        }

        down = down * aprob

      }

      Ind = ifelse(data_in[exposure] == e[k],1,0)
      return( (Ind * up) / down )
    }

    EY = function(status)
    {
      a = status[1]
      data_pre = data_in
      data_pre[exposure] = a
      if (y_type == "discrete")
      {
        EY_pre = predict.glm(Y_model,data_pre,type = "response")
      }else
      {
        EY_pre = predict(Y_model,data_pre)
      }
      return(EY_pre)
    }

    if (method =="Robust")
    {
      para_est = function(status)
      {
        IPW_W = IPW(status)
        RBW_W = RBW(status)
        Reg_Robust_V = Reg_Robust(status)
        REG_V = P_Reg(status)
        EY_V = EY(status)
        Y_V = data_in[outcome]
        RB_est = mean(as.matrix(IPW_W * (Y_V - EY_V) + (RBW_W * EY_V) )) + REG_V - Reg_Robust_V
        REG_est = REG_V
        IPW_est = mean(as.matrix(IPW_W * (Y_V)))
        est_result = list("Robust" = RB_est,"G_computation" = REG_est,"IPW" = IPW_est)
        return(est_result)
      }

      effect_name = paste0("Seq",sep = as.character(c(0:k)))

      IPW_effect = numeric(k+1)
      Robust_effect = numeric(k+1)
      G_computation_effect = numeric(k+1)
      all_parameter = list()

      for (i in c(1:(k+2)))
      {
        all_parameter[[i]] = para_est(all_status[[i]])
      }

      for(i in c(1:(k+1)))
      {
        if (scale == "ratio")
        {
          IPW_effect[i] = all_parameter[[i+1]]$IPW / all_parameter[[i]]$IPW
          G_computation_effect[i] = all_parameter[[i+1]]$G_computation / all_parameter[[i]]$G_computation
          Robust_effect[i] = all_parameter[[i+1]]$Robust / all_parameter[[i]]$Robust
        }else
        {
          IPW_effect[i] = all_parameter[[i+1]]$IPW - all_parameter[[i]]$IPW
          G_computation_effect[i] = all_parameter[[i+1]]$G_computation - all_parameter[[i]]$G_computation
          Robust_effect[i] = all_parameter[[i+1]]$Robust - all_parameter[[i]]$Robust
        }

      }

      if (scale == "ratio")
      {
        IPW_TE =  all_parameter[[k+2]]$IPW/all_parameter[[1]]$IPW
        IPW_effect = c(IPW_effect,IPW_TE)
        names(IPW_effect) = c(effect_name,"TE")
        G_computation_TE =  all_parameter[[k+2]]$G_computation/all_parameter[[1]]$G_computation
        G_computation_effect = c(G_computation_effect,G_computation_TE)
        names(G_computation_effect) = c(effect_name,"TE")
        Robust_TE =  all_parameter[[k+2]]$Robust/all_parameter[[1]]$Robust
        Robust_effect = c(Robust_effect,Robust_TE)
        names(Robust_effect) = c(effect_name,"TE")
      }else
      {
        IPW_TE =  all_parameter[[k+2]]$IPW-all_parameter[[1]]$IPW
        IPW_effect = c(IPW_effect,IPW_TE)
        names(IPW_effect) = c(effect_name,"TE")
        G_computation_TE =  all_parameter[[k+2]]$G_computation-all_parameter[[1]]$G_computation
        G_computation_effect = c(G_computation_effect,G_computation_TE)
        names(G_computation_effect) = c(effect_name,"TE")
        Robust_TE =  all_parameter[[k+2]]$Robust-all_parameter[[1]]$Robust
        Robust_effect = c(Robust_effect,Robust_TE)
        names(Robust_effect) = c(effect_name,"TE")
      }

      output = list("IPW_effect" = IPW_effect,"G_computation_effect" = G_computation_effect,"Robust_effect" = Robust_effect)
      return(output)
    }
    if (method =="IPW")
    {
      para_est = function(status)
      {
        IPW_W = IPW(status)
        Y_V = data_in[outcome]
        IPW_est = mean(as.matrix(IPW_W * (Y_V)))
        return(IPW_est)
      }

      effect_name = paste0("Seq",sep = as.character(c(0:k)))

      IPW_effect = numeric(k+1)
      all_parameter = list()

      for (i in c(1:(k+2)))
      {
        all_parameter[[i]] = para_est(all_status[[i]])
      }

      for(i in c(1:(k+1)))
      {
        if (scale == "ratio")
        {
          IPW_effect[i] = all_parameter[[i+1]] / all_parameter[[i]]
        }else
        {
          IPW_effect[i] = all_parameter[[i+1]] - all_parameter[[i]]
        }
      }
      if (scale == "ratio")
      {
        IPW_TE =  all_parameter[[k+2]]/all_parameter[[1]]
        IPW_effect = c(IPW_effect,IPW_TE)
        names(IPW_effect) = c(effect_name,"TE")
      }else
      {
        IPW_TE =  all_parameter[[k+2]]-all_parameter[[1]]
        IPW_effect = c(IPW_effect,IPW_TE)
        names(IPW_effect) = c(effect_name,"TE")
      }
      out = list()
      out$IPW_effect = IPW_effect
      return(out)
    }
    if (method =="G")
    {
      para_est = function(status)
      {
        REG_est = P_Reg(status)
        return(REG_est)
      }

      effect_name = paste0("Seq",sep = as.character(c(0:k)))

      G_computation_effect = numeric(k+1)
      all_parameter = list()

      for (i in c(1:(k+2)))
      {
        all_parameter[[i]] = para_est(all_status[[i]])
      }

      for(i in c(1:(k+1)))
      {
        if (scale == "ratio")
        {
          G_computation_effect[i] = all_parameter[[i+1]] / all_parameter[[i]]
        }else
        {
          G_computation_effect[i] = all_parameter[[i+1]] - all_parameter[[i]]
        }
      }
      if (scale == "ratio")
      {
        G_computation_TE =  all_parameter[[k+2]]/all_parameter[[1]]
        G_computation_effect = c(G_computation_effect,G_computation_TE)
        names(G_computation_effect) = c(effect_name,"TE")
      }else
      {
        G_computation_TE =  all_parameter[[k+2]]-all_parameter[[1]]
        G_computation_effect = c(G_computation_effect,G_computation_TE)
        names(G_computation_effect) = c(effect_name,"TE")
      }
      out = list()
      out$G_computation_effect = G_computation_effect
      return(out)
    }
  }

  if (boot == TRUE)
  {
    pb = progress_bar$new(format = " [:bar] Bootsratp_replicate: :current / :total :percent Estimted_remaining_time: :eta",
                          clear = FALSE,
                          width = 100,
                          total = boot_rep)
    set.seed(seed_num)
    if (method == "Robust")
    {
      all = mrPF(data)
      IPW_sd = matrix(0,ncol = boot_rep,nrow = n_effect+1)
      Reg_sd = matrix(0,ncol = boot_rep,nrow = n_effect+1)
      RB_sd = matrix(0,ncol = boot_rep,nrow = n_effect+1)
      all_pos = list()
      for (i in c(1:boot_rep))
      {
        pos = sample(c(1:n_sample),size = floor(boot_ratio*n_sample),replace = b_rep)
        all_pos[[i]] = pos
      }
      for(i in c(1:boot_rep))
      {
        res = mrPF(data[all_pos[[i]],])
        IPW_sd[,i] = res$IPW_effect
        Reg_sd[,i] = res$G_computation_effect
        RB_sd[,i] = res$Robust_effect
        pb$tick()
      }



      effect_name = c(paste0("Seq",sep = as.character(c(0:length(mediators)))),"TE")
      IPW_sd_t1 = apply(IPW_sd,1,sd)
      names(IPW_sd_t1) = effect_name
      Reg_sd_t1 = apply(Reg_sd,1,sd)
      names(Reg_sd_t1) = effect_name
      RB_sd_t1 = apply(RB_sd,1,sd)
      names(RB_sd_t1) = effect_name
      IPW_pvalue = 2*pnorm(-abs(all$IPW_effect),0,IPW_sd_t1)
      Reg_pvalue = 2*pnorm(-abs(all$G_computation_effect),0,Reg_sd_t1)
      RB_pvalue = 2*pnorm(-abs(all$Robust_effect),0,RB_sd_t1)

      names(IPW_pvalue) = effect_name
      names(Reg_pvalue) = effect_name
      names(RB_pvalue) = effect_name

      output = list()
      output$IPW_effect = all$IPW_effect
      output$G_computation_effect = all$G_computation_effect
      output$Robust_effect = all$Robust_effect
      output$IPW_sd = IPW_sd_t1
      output$G_computation_sd = Reg_sd_t1
      output$Robust_sd = RB_sd_t1
      output$IPW_pvalue = IPW_pvalue
      output$G_computation_pvalue = Reg_pvalue
      output$Robust_pvalue = RB_pvalue
      output$IPW_Parametric_CI = matrix(cbind(all$IPW_effect-1.96*IPW_sd_t1,all$IPW_effect+1.96*IPW_sd_t1),
                                        ncol = 2,
                                        dimnames = list(effect_name,c("2.5%","97.5%")))
      output$G_computation_Parametric_CI = matrix(cbind(all$G_computation_effect-1.96*Reg_sd_t1,all$G_computation_effect+1.96*Reg_sd_t1),
                                                  ncol = 2,
                                                  dimnames = list(effect_name,c("2.5%","97.5%")))
      output$Robust_Parametric_CI = matrix(cbind(all$Robust_effect-1.96*RB_sd_t1,all$Robust_effect+1.96*RB_sd_t1),
                                           ncol = 2,
                                           dimnames = list(effect_name,c("2.5%","97.5%")))

      output$IPW_nonParametric_CI = matrix(t(apply(IPW_sd,1,function(x) quantile(x,prob = c(0.025,0.975)))),
                                           ncol = 2,
                                           dimnames = list(effect_name,c("2.5%","97.5%")))
      output$G_computation_nonParametric_CI = matrix(t(apply(Reg_sd,1,function(x) quantile(x,prob = c(0.025,0.975)))),
                                                     ncol = 2,
                                                     dimnames = list(effect_name,c("2.5%","97.5%")))
      output$Robust_nonParametric_CI = matrix(t(apply(RB_sd,1,function(x) quantile(x,prob = c(0.025,0.975)))),
                                              ncol = 2,
                                              dimnames = list(effect_name,c("2.5%","97.5%")))

      return(output)
    }
    if (method == "IPW")
    {
      IPW_sd = matrix(0,ncol = boot_rep,nrow = n_effect+1)
      all = mrPF(data)
      all_pos = list()
      for (i in c(1:boot_rep))
      {
        pos = sample(c(1:n_sample),size = floor(boot_ratio*n_sample),replace = b_rep)
        all_pos[[i]] = pos
      }
      for (i in c(1:boot_rep))
      {
        res = mrPF(data[all_pos[[i]],])
        IPW_sd[,i] = res[[1]]
        # cat('Bootstrap ',100*(i/boot_rep),'%\n')
        pb$tick()

      }

      effect_name = c(paste0("Seq",sep = as.character(c(0:length(mediators)))),"TE")
      IPW_sd_t1 = apply(IPW_sd,1,sd)
      names(IPW_sd_t1) = effect_name

      IPW_pvalue = 2*pnorm(-abs(all[[1]]),0,IPW_sd_t1)
      names(IPW_pvalue) = effect_name

      output = list()
      output$IPW_effect = all[[1]]
      output$IPW_sd = IPW_sd_t1
      output$IPW_pvalue = IPW_pvalue
      output$IPW_Parametric_CI = matrix(cbind(all[[1]]-1.96*IPW_sd_t1,all[[1]]+1.96*IPW_sd_t1),
                                        ncol = 2,
                                        dimnames = list(effect_name,c("2.5%","97.5%")))
      output$IPW_nonParametric_CI = matrix(t(apply(IPW_sd,1,function(x) quantile(x,prob = c(0.025,0.975)))),
                                           ncol = 2,
                                           dimnames = list(effect_name,c("2.5%","97.5%")))

      return(output)
    }
    if (method == "G")
    {
      all = mrPF(data)
      Reg_sd = matrix(0,ncol = boot_rep,nrow = n_effect+1)
      all_pos = list()
      for (i in c(1:boot_rep))
      {
        pos = sample(c(1:n_sample),size = floor(boot_ratio*n_sample),replace = b_rep)
        all_pos[[i]] = pos
      }
      for (i in c(1:boot_rep))
      {
        res = mrPF(data[all_pos[[i]],])[[1]]
        Reg_sd[,i] = res
        pb$tick()

      }

      effect_name = c(paste0("Seq",sep = as.character(c(0:length(mediators)))),"TE")

      Reg_sd_t1 = apply(Reg_sd,1,sd)
      names(Reg_sd_t1) = effect_name
      Reg_pvalue = 2*pnorm(-abs(all[[1]]),0,Reg_sd_t1)
      names(Reg_pvalue) = effect_name

      output = list()
      output$G_computation_effect = all[[1]]
      output$G_computation_sd = Reg_sd_t1
      output$G_computation_pvalue = Reg_pvalue
      output$G_computation_Parametric_CI = matrix(cbind(all[[1]]-1.96*Reg_sd_t1,all[[1]]+1.96*Reg_sd_t1),
                                                  ncol = 2,
                                                  dimnames = list(effect_name,c("2.5%","97.5%")))
      output$G_computation_nonParametric_CI = matrix(t(apply(Reg_sd,1,function(x) quantile(x,prob = c(0.025,0.975)))),
                                                     ncol = 2,
                                                     dimnames = list(effect_name,c("2.5%","97.5%")))

      return(output)
    }
  }else
  {
    return(mrPF(data))
  }
}
