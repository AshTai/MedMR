# MedMR
 Multiply robust estimation of natural indirect effects with multiple ordered mediators.
 
To install this package, please copy and paste following codes into your R session:

install.packages("devtools")
library(devtools)
install_github("AshTai/MedMR")

Example

library(MedMR)
data("HCCexample")
ff <- MedMR(data=HCCexample,exposure="a",mediators=c("m1","m2"),outcome="y",
confounder=c("gender","age","smoke","alcohol"),
m_type=c("d","d"),y_type = "discrete",
method = "Robust",
scale = "difference",
boot = T,
boot_ratio = 0.8,
boot_rep = 20,
seed_num = 123,
double_mont = 1e5,
single_mont = 2e5,
b_rep = TRUE)

An-Shun Tai (anshuntai@nctu.edu.tw) https://anshuntai.weebly.com
