
# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(DESeq2)
library(survival)

# src ---------------------------------------------------------------------

source(file='src/doparallel.R', local = TRUE)

# Load data ---------------------------------------------------------------

total416.os.se <- readr::read_rds(file = 'data/rda/total416.os.se.rds.gz')
total434.pfs.se <- readr::read_rds(file = 'data/rda/total434.pfs.se.rds.gz')
Surv(time=total416.os.se$duration, event=total416.os.se$event)
f1 <- survfit(formula = Surv(time = duration, event = event) ~ 1, data = total416.os.se@colData)
names(f1)

plot(f1, xlab='Months', ylab = 'Overall survival probability')
sd <- survdiff(formula = Surv(time = duration, event = event) ~ stage, data = total416.os.se@colData)
1-pchisq(sd$chisq, length(sd$n) - 1)
