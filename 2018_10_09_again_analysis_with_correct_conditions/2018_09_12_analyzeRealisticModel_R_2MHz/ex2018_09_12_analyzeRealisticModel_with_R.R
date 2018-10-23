rate_EMCLs<-read.csv(file("result/2018_08_12_analyzeRealisticModel/rate_EMCLs.csv"),header = FALSE)
ave_SOS<-read.csv(file("result/2018_08_12_analyzeRealisticModel/aveSOS.csv"),header = FALSE)
rate_EMCLs<-t(rate_EMCLs)
ave_SOS<-t(ave_SOS)
m <- lm(ave_SOS ~ rate_EMCLs)
windows()
par(ps=18,mex=1.2)
plot(ave_SOS ~ rate_EMCLs,
     xlab = 'EMCL[%]',
     ylab = '•½‹Ï‰¹‘¬’l[m/s]')
abline(m)
savePlot(filename = "result/2018_09_12_analyzeRealisticModel_with_R/linear_regression.png",type ="png" )
dev.off()
par(ps=18,mex=1.2)
plot(ave_SOS ~ rate_EMCLs,
     xlab = 'EMCL[%]',
     ylab = '•½‹Ï‰¹‘¬’l[m/s]')
abline(m)
m
summary(m)
coef(m)
confint(m)
resid(m)
deviance(m)
anova(m)

