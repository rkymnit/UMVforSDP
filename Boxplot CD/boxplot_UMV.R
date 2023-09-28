#FP boxplot after labeling
label<-c("CLAMI", "CLAMI+",	"ACL",	"TCLP","SQRT","CBRT","UMV",	"KMS",	"NGC",	"MBC",	"HC","NB",	"SVM",	"KNN",	"RF",	"C50")
#label<-c("CLAMI", "CLAMI+",	"ACL",	"TCLP","SQRT","CBRT","UMV",	"KMS",	"NGC",	"MBC",	"HC","NB",	"SVM",	"KNN",	"RF",	"C50")


fm_cls<-read.csv("bw_fm_cls.csv",check.names=FALSE)
#cd_fm<-as.matrix(cd_fm)
acc_cls<-read.csv("bw_acc_cls.csv",check.names=FALSE)
#cd_acc<-as.matrix(cd_acc)
mcc_cls<-read.csv("bw_mcc_cls.csv",check.names=FALSE)
#cd_mcc<-as.matrix(cd_mcc)
#1.Accuracy 
# acc_cls<-read.csv("bw_acc_cls.csv",check.names=FALSE)  #use this for reading symbol of header
# png("bplot_acc_cls.png", width = 15, height = 4, units = 'in', res = 600)
# boxplot(acc_cls,ylab="Accuracy (%)",names=label) # Make plot
# dev.off()


# #1.fm
# fm_cls<-read.csv("bw_fm_cls.csv")
# png("bplot_fm_cls.png", width = 16, height = 4, units = 'in', res = 600)
# boxplot(fm_cls,ylab="F-measure",names=label) # Make plot
# dev.off()
# 
# 
# #1.mcc
# mcc_cls<-read.csv("bw_mcc_cls.csv")
# png("bplot_mcc_cls.png", width = 16, height = 4, units = 'in', res = 600)
# boxplot(mcc_cls,ylab="MCC",names=label)
# dev.off()

#################################################

pdf('bplot_acc_cls.pdf',width = 16)#, height = 3)
par(mar = c(4.2,4.1,0.3,0.3))    #c(bottom,left,top,right)).
boxplot(acc_cls,ylab="Accuracy (%)",names=label,las=2, col=c("gray"),border="brown")
dev.off()

pdf('bplot_fm_cls.pdf',width = 16)
par(mar = c(4.2,4.1,0.3,0.3))    #c(bottom,left,top,right)).
boxplot(fm_cls,ylab="F-measure",names=label,las=2, col=c("gray"),border="brown")
dev.off()


pdf('bw_mcc_cls.pdf',width = 16)
par(mar = c(4.2,4.1,0.3,0.3))    #c(bottom,left,top,right)).
boxplot(mcc_cls,ylab="MCC",names=label,las=2, col=c("gray"),border="brown")
dev.off()

