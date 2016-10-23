setwd("/Users/Scott/Google Drive/Hurwitz Lab/cuffnorm-out")

tlr4 <- read.delim("tlr4-pcr-data.tab")

order<-c('s+h-','s-h-','s+h+','s-h+')

tlr4$Series <- factor(tlr4$Series,order)

boxplot(tlr4$Values ~ tlr4$Series,col='white')

#S-h- vs. S+h-
t.test(tlr4[tlr4$Series=='s+h-',]$Values,tlr4[tlr4$Series=='s-h-',]$Values)

#S-h- vs. S+h+
t.test(tlr4[tlr4$Series=='s-h-',]$Values,tlr4[tlr4$Series=='s+h+',]$Values)

#S+h+ vs. S-h+
t.test(tlr4[tlr4$Series=='s+h+',]$Values,tlr4[tlr4$Series=='s-h+',]$Values)

#S+h- vs. S-h+
t.test(tlr4[tlr4$Series=='s+h-',]$Values,tlr4[tlr4$Series=='s-h+',]$Values)

bfrag <- data.frame(series=c('s+h-','s-h-','s+h+','s-h+'),values=c(2997.993999,4566.114029,5776.721804,6506))
bfrag$series <- factor(bfrag$series,order)

#scaling to min of minimum tlr4 expression
min(tlr4$Values) #0.311
min(bfrag$values) #2997.994
scale_adjust = min(tlr4$Values)-min(bfrag$values)

par(mar = c(5,5,2,5))
boxplot(tlr4$Values ~ tlr4$Series,col='white',ylab='TLR4 relative expression',xlab='Genotype (SMAD3 + or -) and Phenotype (H.hep + or -)')
par(new=T)
plot(x=c(0.5,1.5,2.5,3.5),bfrag$values+scale_adjust, type='l', col='red', axes=F, xlab=NA, ylab=NA,xlim=c(0,4))
axis(side = 4)
mtext(side = 4, line = 3, 'B.fragilis scaled counts')
legend("topleft", legend=c('B.fragilis counts','TLR4 expression'), lty=c(1,0), pch=c(NA, 22), col=c("red", "black"))

#Correlation####
by_tlr4<-by(tlr4$Values,tlr4[,"Series"],median)
#by_tlr4<-by(tlr4$Values,tlr4[,"Series"],mean)
medians_tlr4<-data.frame(series=c('s+h-','s-h-','s+h+','s-h+'),values=c(by_tlr4[[1]],by_tlr4[[2]],by_tlr4[[3]],by_tlr4[[4]]))
cor(medians_tlr4$values,bfrag$values,method = "p")
cor.test(medians_tlr4$values,bfrag$values,method = "p")

b_model <- lm(bfrag$values ~ medians_tlr4$values)
plot(medians_tlr4$values,bfrag$values,ylab='B.fragilis normalized counts',xlab='TLR4 relative expression')
abline(3467,238,col='blue')
legend("bottomright", legend=c('Line of regression','Correlated values', 'pvalue = 0.09862'), lty=c(1,0,0), pch=c(NA, 1,-1), col=c("blue", "black", "black"))

#now for nfkb ####
nfkb <- read.delim("nfkb-pcr-data.tab")

order<-c('s+h-','s-h-','s+h+','s-h+')

nfkb$Series <- factor(nfkb$Series,order)

boxplot(nfkb$Values ~ nfkb$Series,col='white')

#S-h- vs. S+h-
t.test(nfkb[nfkb$Series=='s+h-',]$Values,nfkb[nfkb$Series=='s-h-',]$Values)

#S-h- vs. S+h+
t.test(nfkb[nfkb$Series=='s-h-',]$Values,nfkb[nfkb$Series=='s+h+',]$Values)

#S+h+ vs. S-h+
t.test(nfkb[nfkb$Series=='s+h+',]$Values,nfkb[nfkb$Series=='s-h+',]$Values)

#S+h- vs. S-h+
t.test(nfkb[nfkb$Series=='s+h-',]$Values,nfkb[nfkb$Series=='s-h+',]$Values)

#S+h- vs. S+h+
t.test(nfkb[nfkb$Series=='s+h-',]$Values,nfkb[nfkb$Series=='s+h+',]$Values)

#scaling to min of minimum nfkb expression
min(nfkb$Values) #0.4322686
min(bfrag$values) #2997.994
scale_adjust = min(nfkb$Values)-min(bfrag$values)
par(mar = c(5,5,2,5))
boxplot(nfkb$Values ~ nfkb$Series,col='white',ylab='NFkB relative expression',xlab='Genotype (SMAD3 + or -) and Phenotype (H.hep + or -)')
par(new=T)
plot(x=c(0.5,1.5,2.5,3.5),bfrag$values+scale_adjust, type='l', col='red', axes=F, xlab=NA, ylab=NA,xlim=c(0,4))
axis(side = 4)
mtext(side = 4, line = 3, 'B.fragilis scaled counts')
legend("topleft", legend=c('B.fragilis counts','NFkB expression'), lty=c(1,0), pch=c(NA, 22), col=c("red", "black"))


clusterof2 <- rbind(nfkb,tlr4)
clusterof2$Series <- factor(clusterof2$Series,order)
boxplot(Values ~ Series, data = clusterof2)


#Correlation2####
by_clusterof2<-by(clusterof2$Values,clusterof2[,"Series"],median)
#by_clusterof2<-by(clusterof2$Values,clusterof2[,"Series"],mean)
medians_clusterof2<-data.frame(series=c('s+h-','s-h-','s+h+','s-h+'),values=c(by_clusterof2[[1]],by_clusterof2[[2]],by_clusterof2[[3]],by_clusterof2[[4]]))
cor(medians_clusterof2$values,bfrag$values,method = "p")
cor.test(medians_clusterof2$values,bfrag$values,method = "p")

b_model <- lm(bfrag$values ~ medians_clusterof2$values)
plot(medians_clusterof2$values,bfrag$values,ylab='B.fragilis normalized counts',xlab='TLR4+nfkb median expression')
abline(1620,1323,col='blue')
legend("bottomright", legend=c('Line of regression','Correlated values', 'pvalue = 0.01459'), lty=c(1,0,0), pch=c(NA, 1,-1), col=c("blue", "black", "black"))

#now for tlr2 ####
tlr2 <- read.delim("tlr2-pcr-data.tab")

order<-c('s+h-','s-h-','s+h+','s-h+')

tlr2$Series <- factor(tlr2$Series,order)

boxplot(tlr2$Values ~ tlr2$Series,col='white')

clusterof3 <- rbind(nfkb,tlr4,tlr2)
clusterof3$Series <- factor(clusterof3$Series,order)
boxplot(Values ~ Series, data = clusterof3)

#Correlation3####
by_clusterof3<-by(clusterof3$Values,clusterof3[,"Series"],median)
#by_clusterof3<-by(clusterof3$Values,clusterof3[,"Series"],mean)
medians_clusterof3<-data.frame(series=c('s+h-','s-h-','s+h+','s-h+'),values=c(by_clusterof3[[1]],by_clusterof3[[2]],by_clusterof3[[3]],by_clusterof3[[4]]))
cor(medians_clusterof3$values,bfrag$values,method = "p")
cor.test(medians_clusterof3$values,bfrag$values,method = "p")

b_model <- lm(bfrag$values ~ medians_clusterof3$values)
plot(medians_clusterof3$values,bfrag$values,ylab='B.fragilis normalized counts',xlab='TLR4/TLR2/nfkb relative expression')
abline(2205,1128,col='blue')
legend("bottomright", legend=c('Line of regression','Correlated values', 'pvalue = 0.1247'), lty=c(1,0,0), pch=c(NA, 1,-1), col=c("blue", "black", "black"))
