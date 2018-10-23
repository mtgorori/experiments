rate_IMCLs_reg <- as.vector(as.matrix(rate_IMCLs)) 
rate_EMCLs_reg <- as.vector(as.matrix(rate_EMCLs))
meanSOS_reg <- as.vector(as.matrix(meanSOS))
maxSOS_reg <- as.vector(as.matrix(maxSOS))
minSOS_reg <- as.vector(as.matrix(minSOS))
stdSOS_reg <- as.vector(as.matrix(stdSOS))
skewnessSOS_reg <- as.vector(as.matrix(skewnessSOS))
kurtosisSOS_reg <- as.vector(as.matrix(kurtosisSOS))

multi_reg <- lm(rate_IMCLs_reg ~ rate_EMCLs + meanSOS_reg + maxSOS_reg + minSOS_reg + stdSOS_reg + skewnessSOS_reg + kurtosisSOS_reg)
summary(multi_reg)

