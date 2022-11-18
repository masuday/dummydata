library(car)      # for Anova
library(emmeans)  # for lsmeans
library(sasLM)    # alternatives

data <- read.table("lsm.dat")
colnames(data) <- c("A","cov","dg","B")
data$A <- as.factor(data$A)
data$B <- as.factor(data$B)

# standard method
options(contrasts=c("contr.sum","contr.poly"))
result <- lm(dg ~ A + B + A*B + cov, data=data)
summary(result)
car::Anova(result, type="III")
emmeans::lsmeans(result, ~ A)
emmeans::lsmeans(result, ~ B)
emmeans::lsmeans(result, ~ A*B)
emmeans::lsmeans(result, ~ cov)

# alternative method
sasLM::GLM(dg ~ A + B + A*B + cov, data)
sasLM::aov3(dg ~ A + B + A*B + cov, data)
sasLM::LSM(dg ~ A + B + A*B + cov, data)
