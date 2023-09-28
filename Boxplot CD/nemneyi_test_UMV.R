#install.packages("PMCMR")

library(tsutils)
library(PMCMR)
#interval plot at confidence level 0.95 mean+-CI
label<-c("CLAMI", "CLAMI+",	"ACL",	"TCLP","SQRT","CBRT","UMV",	"KMS",	"NGC",	"MBC",	"HC","NB",	"SVM",	"KNN",	"RF",	"C50")

cd_fm<-read.csv("bw_fm_cls.csv",check.names=FALSE)
cd_fm<-as.matrix(cd_fm)
cd_acc<-read.csv("bw_acc_cls.csv",check.names=FALSE)
cd_acc<-as.matrix(cd_acc)
cd_mcc<-read.csv("bw_mcc_cls.csv",check.names=FALSE)
cd_mcc<-as.matrix(cd_mcc)

#fmeasure
# png("cd_fm.png", width = 5, height = 4, units = 'in', res = 600)
# nr_fm<-nemenyi(cd_fm,conf.level=0.95,plottype = c("mcb"),ylab="Mean Scores (F-measure)",names=label)  #after friedman test
# dev.off()
pdf('cd_fm.pdf')
par(mar = c(4,4,4,0.3))    #c(bottom,left,top,right)).
nr_fm<-nemenyi(cd_fm,conf.level=0.95,plottype = c("mcb"),ylab="Mean Scores (F-measure)",names=label)  #after friedman test
dev.off()

#Accuracy
# png("cd_acc.png", width = 5, height = 4, units = 'in', res = 600)
# nr_acc<-nemenyi(cd_acc,conf.level=0.95,plottype = c("mcb"),ylab="Mean Scores (Accuracy)",names=label)
# dev.off()
pdf('cd_acc.pdf')
par(mar = c(4,4,4,0.3))    #c(bottom,left,top,right)).
nr_acc<-nemenyi(cd_acc,conf.level=0.95,plottype = c("mcb"),ylab="Mean Scores (Accuracy)",names=label)
dev.off()

#MCC
# png("cd_mcc.png", width = 5, height = 4, units = 'in', res = 600)
# nr_mcc<-nemenyi(cd_mcc,conf.level=0.95,plottype = c("mcb"),ylab="Mean Scores (MCC)",names=label)
# dev.off()
pdf('cd_mcc.pdf')
par(mar = c(4,4,4,0.3))    #c(bottom,left,top,right)).
nr_mcc<-nemenyi(cd_mcc,conf.level=0.95,plottype = c("mcb"),ylab="Mean Scores (MCC)",names=label)
dev.off()

pv_acc<-posthoc.friedman.nemenyi.test(cd_acc)
pv_fm<-posthoc.friedman.nemenyi.test(cd_fm)
pv_mcc<-posthoc.friedman.nemenyi.test(cd_mcc)


View(pv_acc)
pv_tclp_phfnt<-matrix(c(pv_acc[["p.value"]],pv_fm[["p.value"]], pv_mcc[["p.value"]]), ncol = 5, byrow = TRUE)
colnames(pv_tclp_phfnt)<-c('SL','KMS','CLA','CLAMI','TCL','TCLP')
rownames(pv_tclp_phfnt)<-c('KMS','CLA','CLAMI','TCL','TCLP')
pv_tclp_phfnt<-as.table(pv_tclp_phfnt)
friedman.test(cd_acc)


