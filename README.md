# MedMR
 Multiply robust estimation of natural indirect effects with multiple ordered mediators.
 To install this package, please copy and paste following codes into your R session:

1. install.packages("devtools")
2. library(devtools)
3. install_github("AshTai/MedMR")


## Example
```R
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
```


## Contact information
An-Shun Tai ([anshuntai@nctu.edu.tw](mailto:anshuntai@nctu.edu.tw))
https://anshuntai.weebly.com
