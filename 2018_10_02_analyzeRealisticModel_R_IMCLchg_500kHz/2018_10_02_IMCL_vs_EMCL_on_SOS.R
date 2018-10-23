# library(ggplot2)
# library(ggsci)
# g <- ggplot(t(rate_IMCLs[,1]),t(aveSOS2[,1]),mapping =aes(x = IMCL占有率,y = 平均音速値))
# g <- g + geom_line()
# plot(g)

#IMCL11種用意．平均音速をプロット．個別に
for (i in 1:25) {
  windows()
  par(ps=18,mex=1.2,pty='s')
  plot(t(rate_IMCLs[,i]),t(aveSOS2[,i]),
       xlab = 'IMCL[%]',
       ylab = '平均音速値[m/s]',
       main = sprintf('EMCL = %0.2f %%',rate_EMCLs[1,i]))
  savefilename <- sprintf('//Azlab-fs01/東研究室/個人work/竹内(ひ)/result/2018_10_02_analyzeRealisticModel_R_IMCLchg_500kHz/IMCL_vs_aveSOS_%d.png',i)
  savePlot(filename = savefilename,type="png")
  dev.off()  
}

#IMCL11種用意．平均音速をプロット．一括で．
windows()
par(mfrow=c(5,5),pty='s')
for (i in 1:25) {
  plot(t(rate_IMCLs[,i]),t(aveSOS2[,i]),
       xlab = 'IMCL[%]',
       ylab = '平均音速値[m/s]',
       main = sprintf('EMCL = %0.2f %%',rate_EMCLs[1,i]))
}
savefilename <- sprintf('//Azlab-fs01/東研究室/個人work/竹内(ひ)/result/2018_10_02_analyzeRealisticModel_R_IMCLchg_500kHz/IMCL_vs_aveSOS_all.png')
savePlot(filename = savefilename,type="png")
dev.off()

#上のやり方では細かすぎるので，一つのグラフに重ねてプロットする．
windows()
par(ps=18,mex=1.2,pty='s')
matplot(rate_IMCLs[,1:25]*0.1,aveSOS2[,1:25],
        col = 1,
        type = 'o',
        pch = 1:25,
        xlab = 'IMCL[%]',
        ylab = '平均音速値[m/s]')
savefilename <- sprintf('//Azlab-fs01/東研究室/個人work/竹内(ひ)/result/2018_10_02_analyzeRealisticModel_R_IMCLchg_500kHz/IMCL_vs_aveSOS_sum.png')
savePlot(filename = savefilename,type="png")
dev.off()

#横軸：EMCL，縦軸：平均音速値，一つのグラフに重ねてプロットする．
windows()
par(ps=18,mex=1.2,pty='s')
matplot(t(rate_EMCLs),t(aveSOS2),
        col = 1,
        type = 'o',
        pch = 1:11,
        xlab = 'EMCL[%]',
        ylab = '平均音速値[m/s]')
savefilename <- sprintf('//Azlab-fs01/東研究室/個人work/竹内(ひ)/result/2018_10_02_analyzeRealisticModel_R_IMCLchg_500kHz/EMCL_vs_aveSOS_sum.png')
savePlot(filename = savefilename,type="png")
dev.off()

#横軸：IMCL，縦軸：EMCL，一つのグラフに重ねてプロットする．
windows()
par(ps=18,mex=1.2,pty='s')
matplot(t(rate_IMCLs),t(rate_EMCLs),
        col = 1,
        type = 'o',
        pch = 1:11,
        xlab = 'IMCL[%]',
        ylab = 'EMCL[%]')
savefilename <- sprintf('//Azlab-fs01/東研究室/個人work/竹内(ひ)/result/2018_10_02_analyzeRealisticModel_R_IMCLchg_500kHz/IMCL_vs_EMCL_sum.png')
savePlot(filename = savefilename,type="png")
dev.off()

#water-fallプロットもしくはメッシュプロットで，IMCLとEMCLが平均音速値に与える影響を見てみる．
#使い方に時間を掛けすぎている．matlabで行う．
library(scatterplot3d)
window()
par(ps=18,mex=1.2,pty='s')
scatter
scatterplot3d(t(rate_EMCLs),t(rate_IMCLs)*0.1,t(aveSOS2))

# 横軸にIMCL率，縦軸に平均音速値をとり，全サンプル分プロットする．
# これらプロットデータに対して回帰直線を求めて回帰係数を求める．
delta_SOS_by_IMCL1 <- matrix(nrow = 25,ncol = 1)
for (ii in 1:25) {
x_tmp = t(rate_IMCLs[,ii])*0.1
y_tmp = t(aveSOS2[,ii])
result <- lm(t(y_tmp) ~ t(x_tmp))
delta_SOS_by_IMCL1[ii,1] <- as.vector(result$coefficients[2])
}
ave_delta_SOS_by_IMCL1 <- mean(delta_SOS_by_IMCL1)

# 横軸にIMCL率，縦軸に平均音速値をとり，全サンプル分プロットする．
# これらプロットデータに対して回帰直線を求めて回帰係数を求める．
delta_SOS_by_EMCL1 <- matrix(nrow = 11,ncol = 1)
for (ii in 1:11) {
  x_tmp = rate_EMCLs[ii,]
  y_tmp = aveSOS2[ii,]
  result <- lm(t(y_tmp) ~ t(x_tmp))
  delta_SOS_by_EMCL1[ii,1] <- as.vector(result$coefficients[2])
}
ave_delta_SOS_by_EMCL1 <- mean(delta_SOS_by_EMCL1)