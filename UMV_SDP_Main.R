
set.seed(786)
packageslist <- list("tidyverse", "cluster", "factoextra", "e1071", "caret", "gridExtra", "fBasics", "dplyr", "plyr", "ggplot2", "reshape2", "VGAM", "corrplot","mlbench", "MASS", "RColorBrewer", "RANN", "ROCR", "pROC", "MLmetrics", "imputeTS", "xgboost","mccr"
                     ,"kohonen","vegan","mclust","cclust","FCPS","DMwR","ROSE")
load_packages <- lapply(packageslist, require, character.only = T)


getmode<-function(x){
  uniq<-unique(x)
  uniq[which.max(tabulate(match(x,uniq)))]
}

ONS=0.99    
hp1_const_sr_tclp<-ONS
hp1_const_sr_ltp<-ONS
hp1_const_sr_srtp<-ONS
hp1_const_sr_crtp<-ONS    #Means actual clami+

hp1_org<-read.csv("etherium_defect.csv")
str(hp1_org)
#hp1_org<-hp1_org[,1:6]
#hp1_org$bug<-as.factor(hp1_org$bug>0)
#levels(hp1_org$bug) <- c("NO", "YES")

hp1_org_ncol<-ncol(hp1_org)
hp1_org_ncol
hp1_org_nrow<-nrow(hp1_org)
hp1_org_nrow
hp1_org_count<-count(hp1_org$bug)
hp1_org_count
hp1_org_bugpr<-(hp1_org_count[2,2]/hp1_org_nrow)*100   #bug % module in dataset
hp1_org_bugpr

hp1_01b<-hp1_org
levels(hp1_01b$bug)<-0:1  #level bug as 0 and 1 for No and YES
str(hp1_01b)
# remove the label bug from hp1_org, it become hp1_wb
hp1_wb<-hp1_org[,-ncol(hp1_org)]

# hist(hp1_wb$CC)   #Right skewed
# hp1_skew_org<-skewness(hp1_wb)
# hp1_skew_org




hp1_wb_lt<-hp1_wb    #threshold using log10
hp1_wb_ltp<-hp1_wb
hp1_wb_tcl<-hp1_wb
hp1_wb_tclp<-hp1_wb
hp1_wb_srt<-hp1_wb
hp1_wb_srtp<-hp1_wb
hp1_wb_crt<-hp1_wb
hp1_wb_crtp<-hp1_wb




hp1_sl_cv<-hp1_org
# 
#cor1<-cor(hp1_wb_km)
####################################
## Threshold calculation start and implement tcl
#Check with calculate log without adding 1
hp1_wb_tcl1<-hp1_wb_tcl+1  #Add 1 to the metrics so that easily calculated the log of the metrics. because log 0 is not defined
hp1_wb_tcl1[,c(1:ncol(hp1_wb_tcl1))] <- log(hp1_wb_tcl1[,c(1:ncol(hp1_wb_tcl1))])   #take the log of all metrics

hp1_wb_tcl1_mean<- colMeans(hp1_wb_tcl1)
hp1_wb_tcl1_median<- apply(hp1_wb_tcl1,2,median)  #2 is margin means col, sd is function
hp1_wb_tcl1_mode<- apply(hp1_wb_tcl1,2,getmode)  #2 is margin means col, sd is function
hp1_wb_tcl1_sd<- apply(hp1_wb_tcl1,2,sd)  #2 is margin means col, sd is function

#===============================================================================================================
# now calculate threshold of each  metric for the logrithmic value of metrics given metric values thr=exp(T'), T'=mean+sd
hp1_wb_tcl1_th_m<-c()
for(i in 1:ncol(hp1_wb_tcl1))
{
  hp1_wb_tcl1_th_m[i]<-hp1_wb_tcl1_mean[i]+hp1_wb_tcl1_sd[i]
}
hp1_wb_tcl1_thr_exp<-exp(hp1_wb_tcl1_th_m) #Actual threshold used for tcl methods
hp1_wb_tcl1_threshold<-hp1_wb_tcl1_thr_exp  #Actual threshold
#=================Derived threshold completed====================================================
#Now for better ment mean median mode should be aomost equal to the threshold
hp1_wb_tcl_basic_stats<-basicStats(hp1_wb_tcl)
hp1_wb_tcl1_basic_stats<-basicStats(hp1_wb_tcl1)
hp1_wb_tcl_mean<- colMeans(hp1_wb_tcl)
hp1_wb_tcl_median<- apply(hp1_wb_tcl,2,median)  #2 is margin means col, sd is function
hp1_wb_tcl_mode<- apply(hp1_wb_tcl,2,getmode)
hp1_wb_tcl_skewness<-matrix(c(hp1_wb_tcl_basic_stats[15,],hp1_wb_tcl1_basic_stats[15,],hp1_wb_tcl_mean,hp1_wb_tcl_median,hp1_wb_tcl_mode,hp1_wb_tcl1_threshold),ncol = ncol(hp1_wb_tcl1), byrow = TRUE)
rownames(hp1_wb_tcl_skewness)<- c('Skewness_before_log','Skewness_after_log','mean','median','mode','threshold')
#write.csv(hp1_wb_tcl_skewness,"hp1_wb_tcl_skewness_bfrLog_aftrLog_mean_med_mode_thrs.csv")
#=====================================================================================
#Find the value of k for each instances and and label as YES or NO
hp1_wb_tcl_nRowsDf<-nrow(hp1_wb_tcl)
hp1_wb_tcl_nColsDf<-ncol(hp1_wb_tcl)
hp1_wb_tcl_thr<- hp1_wb_tcl1_threshold  #use threshold calculated using Log transform methods

#Now calculate the values of k for each instances
for(i in 1:hp1_wb_tcl_nRowsDf)
{
  loop=0
  for(j in 1:hp1_wb_tcl_nColsDf)
  {
    if(hp1_wb_tcl[i,j]>=hp1_wb_tcl_thr[j])
    {
      loop=loop+1
    }
    hp1_wb_tcl$k[i]<-loop
  }
}

#Lable top half as YES and bottom half as NO, find max of k values then divide by 2 and get floor value.
for(i in 1:hp1_wb_tcl_nRowsDf)
{
  if(hp1_wb_tcl$k[i]>floor(max(hp1_wb_tcl$k)/2))
  {
    hp1_wb_tcl$label_bug[i]='YES'
  }
  else {
    hp1_wb_tcl$label_bug[i]='NO'
  }
}
#View(hp1_wb_tcl)
#K value calculation and Labeling completed Now compute performance of TCL
hp1_wb_tcl_cm<-confusionMatrix(factor(hp1_org$bug), factor(hp1_wb_tcl$label_bug), mode="prec_recall")
#View(hp1_wb_tcl$label_bug)
hp1_wb_tcl_cm

hp1_wb_tcl_01<-hp1_wb_tcl       #store the actial data in temp data with labelin with 0 1 as bug
hp1_wb_tcl_01$label_bug<-ifelse(hp1_wb_tcl_01$label_bug == "YES", 1,0)
hp1_wb_tcl_mcc<-mccr(hp1_01b$bug,hp1_wb_tcl_01$label_bug)
hp1_wb_tcl_mcc



#Silhouette Plot and AvgSilWidth
hp1_wb_tcl$bnum<-ifelse(hp1_wb_tcl$label_bug=='YES',2,1)
hp1_wb_tcl_dist<-dist(hp1_wb_tcl[,1:6],method = "euclidean")^2
hp1_wb_tcl_sil=silhouette(hp1_wb_tcl$bnum,hp1_wb_tcl_dist)
set.seed(21)
pdf("hp1_wb_tcl_sil_plot.pdf")
plot(hp1_wb_tcl_sil)  #0.78 for eth)
dev.off()
hp1_wb_tcl_sil_summary<-summary(hp1_wb_tcl_sil, FUN=mean)


#Melting data according to label_bug
hp1_wb_tcl_melt<- melt(hp1_wb_tcl[,c(1,2,3,4,5,6,8)], id.var = "label_bug")
hp1_wb_tcl_melt$label_bug<-as.factor(hp1_wb_tcl_melt$label_bug)
pdf("hp1_wb_tcl_bp_label_bug_wise.pdf")
p <- ggplot(data = hp1_wb_tcl_melt, aes(x=variable, y=value)) + geom_boxplot(aes(fill=label_bug))
p + facet_wrap( ~ variable, scales="free")   #Boxplot as per requirment need to scale y axis
dev.off()

#============== TCL completed   ========================================


#Start TCLP
##========   Metric selection  started  ==============================================================
#Now perform Metric selection based on metric violation score MVS (metric violation score), vt (violation table)
#MVS of i metric=Ci/Fi,  Ci is no. of violation in ith metiric and Fi is the no. of metric values in the i-th metric(no. of instances)
#After labeling the instances now perform metric selection
hp1_wb_tcl_nRowsDf<-nrow(hp1_wb_tclp)
hp1_wb_tcl_nColsDf<-ncol(hp1_wb_tclp)
hp1_wb_tcl_vt<-matrix(nrow=hp1_wb_tcl_nRowsDf,ncol=hp1_wb_tcl_nColsDf)
for(i in 1:hp1_wb_tcl_nColsDf)
{
  hp1_wb_tcl_cy<-0  #no. of yes in coloumn cy
  hp1_wb_tcl_cn<-0
  for(j in 1:hp1_wb_tcl_nRowsDf)
  {
    if( hp1_wb_tcl$label_bug[j]=='YES')
    {
      r <- if(hp1_wb_tcl[j,i] <= hp1_wb_tcl_thr[i]) 1 else 0
      hp1_wb_tcl_cy<-hp1_wb_tcl_cy+r
      hp1_wb_tcl_vt[j,i]<-hp1_wb_tcl[j,i]-hp1_wb_tcl_thr[i]     #its for clami+ to write difference
    }
    else if( hp1_wb_tcl$label_bug[j]=='NO')
    {
      s<- if(hp1_wb_tcl[j,i] > hp1_wb_tcl_thr[i]) 1 else 0

      hp1_wb_tcl_cn<-hp1_wb_tcl_cn+s
      hp1_wb_tcl_vt[j,i]<-hp1_wb_tcl[j,i]-hp1_wb_tcl_thr[i]  #its difference table used in clami+

    }
  }
  hp1_wb_tcl[hp1_wb_tcl_nRowsDf+1,i]<-(hp1_wb_tcl_cy+hp1_wb_tcl_cn)  #the last row is added as the no. of violation row for each metrics. correct
}

#Now arrange the metrics in ascending order based on the value of last row and choose first ceilling(log2(n)) metrics for next phase.
#Firstly remove the last coloumn as label_bug for ordering metrics because its factor coloumn
hp1_wb_tcl2<-hp1_wb_tcl[,1:ncol(hp1_wb_tcl)-1]
hp1_wb_tcl2=hp1_wb_tcl2[,order(hp1_wb_tcl2[nrow(hp1_wb_tcl2),])]

# select the metrics with minimum mvs score
hp1_wb_tcl_mvs_min<-min(hp1_wb_tcl[hp1_wb_tcl_nRowsDf+1,1:hp1_wb_tcl_nColsDf])  #min of 8th row
#Select the metrics which has 8th row contain the hp1_wb_tcl_mvs_min value and store it in hp1_wb_tcl_mvs means dataset with select metrics with min mvs
hp1_wb_tcl_nmvs<- sum(hp1_wb_tcl_mvs_min == hp1_wb_tcl[hp1_wb_tcl_nRowsDf+1,1:hp1_wb_tcl_nColsDf])  #number of minimum vs
cat("number of selected metrics=",hp1_wb_tcl_nmvs)
hp1_wb_tcl_mvs<-matrix()
hp1_wb_tcl_mvs_row<-hp1_wb_tcl[hp1_wb_tcl_nRowsDf+1,]  #last row as MVS values of each metrics
hp1_wb_tcl_mvs_row_order<-hp1_wb_tcl_mvs_row[,order(hp1_wb_tcl_mvs_row)]


hp1_wb_tcl_wmvs<-hp1_wb_tcl[1:hp1_wb_tcl_nRowsDf,]   #Removing the last row means mvs row from hp1_wb_tcl, without mvs row
hp1_wb_tcl_thr_st<-c()  #st means selected threshold of selected metrics

#Select n% of metrics using hp1_wb_tcl_select_nm


hp1_wb_tcl_select_nm<-ceiling(log2(hp1_wb_tcl_nColsDf))  #Select first ceiling(log2(n)) metrics based on mvs
hp1_wb_tcl_ams<-matrix()   #hp1_wb_tcl after metric selection
for(i in 1:hp1_wb_tcl_select_nm)
{
  hp1_wb_tcl_ams<-cbind(hp1_wb_tcl_ams,hp1_wb_tcl2[i])   #hp1_wb_tcl2 used
}

hp1_wb_tcl_ams<-cbind(hp1_wb_tcl_ams[2:ncol(hp1_wb_tcl_ams)],label_bug=hp1_wb_tcl$label_bug)

str(hp1_wb_tcl_ams)
hp1_wb_tcl_ams_colnames<-colnames(hp1_wb_tcl_ams[,1:hp1_wb_tcl_select_nm])
hp1_wb_tcl_mvs_row_order_first3<-hp1_wb_tcl_mvs_row_order[,1:hp1_wb_tcl_select_nm]
#write.csv(as.data.frame(hp1_wb_tcl_mvs_row_order_first3),"hp1_wb_tcl_ams_first3.csv")
#write.csv(hp1_wb_tcl_mvs_row_order,'hp1_wb_tcl_mvs_row_order.csv')
##=======================================================================================
##metric selection completed now work on instance selection on 18.00, 03/01/2021
hp1_wb_tcl_ams_nColsDf<-ncol(hp1_wb_tcl_ams)-1
hp1_wb_tcl_ams_nRowsDf<-nrow(hp1_wb_tcl_ams)-1

for(i in 1:hp1_wb_tcl_ams_nRowsDf)
{
  hp1_DS_iy<-0  #no. of yes in coloumn cy
  hp1_DS_in<-0
  for(j in 1:hp1_wb_tcl_ams_nColsDf)
  {
    if( hp1_wb_tcl_ams$label_bug[i]=='YES')
    {
      r <- if(hp1_wb_tcl_ams[i,j] <= hp1_wb_tcl_ams[hp1_wb_tcl_ams_nRowsDf+1,j]) 1 else 0
      hp1_DS_iy<-hp1_DS_iy+r
      #hp1_DS_vt[j,i]<-hp1_wb_tcl_ams[j,i]-hp1_DS_thr_st[i]     #its for clami+ to write difference
    }
    else if( hp1_wb_tcl_ams$label_bug[i]=='NO')
    {
      s<- if(hp1_wb_tcl_ams[i,j] > hp1_wb_tcl_ams[hp1_wb_tcl_ams_nRowsDf+1,j]) 1 else 0

      hp1_DS_in<-hp1_DS_in+s
      #hp1_DS_ivt[j,i]<-hp1_wb_tcl_ams[j,i]-hp1_DS_thr[i]
    }
  }
  hp1_wb_tcl_ams[i,hp1_wb_tcl_ams_nColsDf+2]<-(hp1_DS_iy+hp1_DS_in)  #the 4th coloumn is the no. of violation row for each instances correct done 13.1 6/1/2021
}
#View(hp1_wb_tcl_ams)
count(hp1_wb_tcl_ams[,ncol(hp1_wb_tcl_ams)])  #instance violation score completed   data[,ncol(data)]
## UP to here it executed completely  its final ds after metric selection
################################################################################
#Now arrange the instances in ascending order based on the value of last coloumn and choose first n-n/12 metrics for next phase.
#Firstly remove the last coloumn as label_bug for ordering rows because its factor coloumn
hp1_wb_tcl_ams_row_ord<-hp1_wb_tcl_ams[order(hp1_wb_tcl_ams[,ncol(hp1_wb_tcl_ams)]),]
hp1_wb_tcl_ams_row_ord_nrow<-nrow(hp1_wb_tcl_ams)
hp1_wb_tcl_ams_row_ord<-hp1_wb_tcl_ams_row_ord[-hp1_wb_tcl_ams_row_ord_nrow,]
str(hp1_wb_tcl_ams_row_ord)

hp1_wb_tcl_ams_cnv<-count(hp1_wb_tcl_ams_row_ord$V5)   #count the number of violations existed in the last coloumn V5
hp1_wb_tcl_ams_cnv


#Now select top n% of instances using hp1_const_sr_tclp


hp1_wb_tcl_ams_row_ord<-head(hp1_wb_tcl_ams_row_ord,hp1_const_sr_tclp*hp1_wb_tcl_ams_row_ord_nrow)   #select no. of rows
hp1_wb_tcl_ams_is_nc<-count(hp1_wb_tcl_ams_row_ord$label_bug)  #count the number of yes and no
hp1_wb_tcl_ams_is_nc
hp1_wb_tcl_ams_is_nc_r<-hp1_wb_tcl_ams_is_nc[1,2]/hp1_wb_tcl_ams_is_nc[2,2]  #ratio
hp1_wb_tcl_ams_is_nc_r

hp1_wb_tcl_ams_row_ord <- hp1_wb_tcl_ams_row_ord[sample(nrow(hp1_wb_tcl_ams_row_ord)),]   #randomize instances
hp1_wb_tcl_ams_ais<-hp1_wb_tcl_ams_row_ord[,-ncol(hp1_wb_tcl_ams_row_ord)]                  #Remove last coloumn as ivs

hp1_wb_tcl_ams_ais_nrow_ratio<-cbind(nrow(hp1_wb_tcl_ams_ais),hp1_wb_tcl_ams_is_nc_r)
#write.csv(as.data.frame(hp1_wb_tcl_ams_ais_nrow_ratio),"hp1_wb_tcl_ams_ais_nrow_ratio.csv")
#View(hp1_wb_tcl_ams_ais)  #This is the final input data
######################################################################

#Now apply the SMOTE: Synthetic Minority Oversampling Technique To Handle Class Imbalancy In Binary Classification
count(hp1_wb_tcl_ams_ais$label_bug)
hp1_wb_tcl_ams_ais_rose<-ROSE(label_bug~.,data=hp1_wb_tcl_ams_ais,seed = 786)$data   #it is better than smote
count(hp1_wb_tcl_ams_ais_rose$label_bug)
set.seed(200)
hp1_wb_tcl_ams_ais_rose<-hp1_wb_tcl_ams_ais_rose[sample(nrow(hp1_wb_tcl_ams_ais_rose)),]

#=======================================================================================================================
# #data<-read.csv(file.choose())  *********change 1:#m,#m+2
hp1_wb_tcl_data<-hp1_wb_tcl_ams_ais_rose     #label name   "label_bug"
hp1_tclp_org<-hp1_org                      #label name "bug"
hp1_wb_tcl_predictors<-colnames(hp1_wb_tcl_data[,1:ncol(hp1_wb_tcl_data)-1])
#write.csv(hp1_wb_tcl_predictors,"hp1_wb_tcl_predictors.csv")
hp1_tclp_org_predictors<-colnames(hp1_tclp_org[,1:ncol(hp1_tclp_org)-1])
hp1_wb_tcl_outcomeName<-'label_bug'
hp1_tclp_org_outcomeName<-'bug'

fitControl <- trainControl(      #train control for TCLP
  method = "cv",
  number = 10,
  savePredictions = 'final',  summaryFunction = multiClassSummary,
  classProbs = T )

set.seed(3)
hp1_wb_tcl_model_rf<-train(hp1_wb_tcl_data[,hp1_wb_tcl_predictors],hp1_wb_tcl_data[,hp1_wb_tcl_outcomeName],method='rf',trControl=fitControl,tuneLength=5)
hp1_wb_tcl_model_rf[["resample"]]<-na.mean(hp1_wb_tcl_model_rf[["resample"]])
hp1_wb_tcl_model_rf
hp1_tclp_org$predrf_bug<-predict(object = hp1_wb_tcl_model_rf,hp1_tclp_org[,hp1_tclp_org_predictors],type='raw')   #predict on original data using pre trained model
hp1_wb_tcl_cm_rf<-caret::confusionMatrix(factor(hp1_tclp_org$bug),factor(hp1_tclp_org$predrf_bug), mode= "prec_recall")
hp1_wb_tcl_cm_rf
#View(hp1_tclp_org$predrf_bug)

hp1_wb_tclp_01<-hp1_tclp_org       #store the actual data in temp data with labelin with 0 1 as bug
hp1_wb_tclp_01$predrf_bug<-ifelse(hp1_wb_tclp_01$predrf_bug == "YES", 1,0)
hp1_wb_tclp_mcc<-mccr(hp1_01b$bug,hp1_wb_tclp_01$predrf_bug)



#Silhouette Plot and AvgSilWidth
hp1_tclp_org$bnum<-ifelse(hp1_tclp_org$predrf_bug=='YES',2,1)
hp1_tclp_org_dist<-dist(hp1_tclp_org[,1:6],method = "euclidean")^2
hp1_tclp_org_sil=silhouette(hp1_tclp_org$bnum,hp1_tclp_org_dist)
set.seed(21)
pdf("hp1_tclp_org_sil_plot.pdf")
plot(hp1_tclp_org_sil)  #0.78 for eth)
dev.off()
hp1_tclp_org_sil_summary<-summary(hp1_tclp_org_sil, FUN=mean)


#Melting data according to predrf_bug
hp1_tclp_org_melt<- melt(hp1_tclp_org[,c(1,2,3,4,5,6,8)], id.var = "predrf_bug")
hp1_tclp_org_melt$predrf_bug<-as.factor(hp1_tclp_org_melt$predrf_bug)
pdf("hp1_tclp_org_bp_predrf_bug_wise.pdf")
p <- ggplot(data = hp1_tclp_org_melt, aes(x=variable, y=value)) + geom_boxplot(aes(fill=predrf_bug))
p + facet_wrap( ~ variable, scales="free")   #Boxplot as per requirment need to scale y axis
dev.off()
#######################################################################################################


###############################################################################

#============== TCLP completed   ========================================

#Start TCL with square root


###############################################################################
#A function used in the hglm package for the inverse square root family. usage
#inverse.sqrt()  for inverse use ^2 or multiply itself sqrt producing 0.94 acc on eth

## Threshold calculation start and implement srt
#Check with calculate log without adding 1
hp1_wb_srt1<-hp1_wb_srt+1  #Add 1 to the metrics so that easily calculated the log of the metrics. because log 0 is not defined
hp1_wb_srt1[,c(1:ncol(hp1_wb_srt1))] <- sqrt(hp1_wb_srt1[,c(1:ncol(hp1_wb_srt1))])   #take the log of all metrics

hp1_wb_srt1_mean<- colMeans(hp1_wb_srt1)
hp1_wb_srt1_median<- apply(hp1_wb_srt1,2,median)  #2 is margin means col, sd is function
hp1_wb_srt1_mode<- apply(hp1_wb_srt1,2,getmode)  #2 is margin means col, sd is function
hp1_wb_srt1_sd<- apply(hp1_wb_srt1,2,sd)  #2 is margin means col, sd is function

#===============================================================================================================
# now calculate threshold of each  metric for the logrithmic value of metrics given metric values thr=exp(T'), T'=mean+sd
hp1_wb_srt1_th_m<-c()
for(i in 1:ncol(hp1_wb_srt1))
{
  hp1_wb_srt1_th_m[i]<-hp1_wb_srt1_mean[i]+hp1_wb_srt1_sd[i]
}
hp1_wb_srt1_thr_exp<-(hp1_wb_srt1_th_m)^2 #Actual threshold used for srt methods
hp1_wb_srt1_threshold<-hp1_wb_srt1_thr_exp  #Actual threshold
#=================Derived threshold completed====================================================
#Now for better ment mean median mode should be aomost equal to the threshold
hp1_wb_srt_basic_stats<-basicStats(hp1_wb_srt)
hp1_wb_srt1_basic_stats<-basicStats(hp1_wb_srt1)
hp1_wb_srt_mean<- colMeans(hp1_wb_srt)
hp1_wb_srt_median<- apply(hp1_wb_srt,2,median)  #2 is margin means col, sd is function
hp1_wb_srt_mode<- apply(hp1_wb_srt,2,getmode)
hp1_wb_srt_skewness<-matrix(c(hp1_wb_srt_basic_stats[15,],hp1_wb_srt1_basic_stats[15,],hp1_wb_srt_mean,hp1_wb_srt_median,hp1_wb_srt_mode,hp1_wb_srt1_threshold),ncol = ncol(hp1_wb_srt1), byrow = TRUE)
rownames(hp1_wb_srt_skewness)<- c('Skewness_before_log','Skewness_after_log','mean','median','mode','threshold')
write.csv(hp1_wb_srt_skewness,"hp1_wb_srt_skewness_bfrLog_aftrLog_mean_med_mode_thrs.csv")
#=====================================================================================
#Find the value of k for each instances and and label as YES or NO
hp1_wb_srt_nRowsDf<-nrow(hp1_wb_srt)
hp1_wb_srt_nColsDf<-ncol(hp1_wb_srt)
hp1_wb_srt_thr<- hp1_wb_srt1_threshold  #use threshold calculated using Log transform methods

#Now calculate the values of k for each instances
for(i in 1:hp1_wb_srt_nRowsDf)
{
  loop=0
  for(j in 1:hp1_wb_srt_nColsDf)
  {
    if(hp1_wb_srt[i,j]>=hp1_wb_srt_thr[j])
    {
      loop=loop+1
    }
    hp1_wb_srt$k[i]<-loop
  }
}

#Lable top half as YES and bottom half as NO, find max of k values then divide by 2 and get floor value.
for(i in 1:hp1_wb_srt_nRowsDf)
{
  if(hp1_wb_srt$k[i]>floor(max(hp1_wb_srt$k)/2))
  {
    hp1_wb_srt$label_bug[i]='YES'
  }
  else {
    hp1_wb_srt$label_bug[i]='NO'
  }
}
#K value calculation and Labeling completed Now compute performance of TCL
hp1_wb_srt_cm<-confusionMatrix(factor(hp1_org$bug), factor(hp1_wb_srt$label_bug), mode="prec_recall")
#View(hp1_wb_srt$label_bug)
hp1_wb_srt_cm

hp1_wb_srt_01<-hp1_wb_srt       #store the actial data in temp data with labelin with 0 1 as bug
hp1_wb_srt_01$label_bug<-ifelse(hp1_wb_srt_01$label_bug == "YES", 1,0)
hp1_wb_srt_mcc<-mccr(hp1_01b$bug,hp1_wb_srt_01$label_bug)
hp1_wb_srt_mcc



#Silhouette Plot and AvgSilWidth
hp1_wb_srt$bnum<-ifelse(hp1_wb_srt$label_bug=='YES',2,1)
hp1_wb_srt_dist<-dist(hp1_wb_srt[,1:6],method = "euclidean")^2
hp1_wb_srt_sil=silhouette(hp1_wb_srt$bnum,hp1_wb_srt_dist)
set.seed(21)
pdf("hp1_wb_srt_sil_plot.pdf")
plot(hp1_wb_srt_sil)  
dev.off()
hp1_wb_srt_sil_summary<-summary(hp1_wb_srt_sil, FUN=mean)


#Melting data according to label_bug
hp1_wb_srt_melt<- melt(hp1_wb_srt[,c(1,2,3,4,5,6,8)], id.var = "label_bug")
hp1_wb_srt_melt$label_bug<-as.factor(hp1_wb_srt_melt$label_bug)
pdf("hp1_wb_srt_bp_label_bug_wise.pdf")
p <- ggplot(data = hp1_wb_srt_melt, aes(x=variable, y=value)) + geom_boxplot(aes(fill=label_bug))
p + facet_wrap( ~ variable, scales="free")   #Boxplot as per requirment need to scale y axis
dev.off()

#============== TCL completed   ========================================


#Start TCLP
##========   Metric selection  started  ==============================================================
#Now perform Metric selection based on metric violation score MVS (metric violation score), vt (violation table)
#MVS of i metric=Ci/Fi,  Ci is no. of violation in ith metiric and Fi is the no. of metric values in the i-th metric(no. of instances)
#After labeling the instances now perform metric selection
hp1_wb_srt_nRowsDf<-nrow(hp1_wb_srtp)
hp1_wb_srt_nColsDf<-ncol(hp1_wb_srtp)
hp1_wb_srt_vt<-matrix(nrow=hp1_wb_srt_nRowsDf,ncol=hp1_wb_srt_nColsDf)
for(i in 1:hp1_wb_srt_nColsDf)
{
  hp1_wb_srt_cy<-0  #no. of yes in coloumn cy
  hp1_wb_srt_cn<-0
  for(j in 1:hp1_wb_srt_nRowsDf)
  {
    if( hp1_wb_srt$label_bug[j]=='YES')
    {
      r <- if(hp1_wb_srt[j,i] <= hp1_wb_srt_thr[i]) 1 else 0
      hp1_wb_srt_cy<-hp1_wb_srt_cy+r
      hp1_wb_srt_vt[j,i]<-hp1_wb_srt[j,i]-hp1_wb_srt_thr[i]     #its for clami+ to write difference
    }
    else if( hp1_wb_srt$label_bug[j]=='NO')
    {
      s<- if(hp1_wb_srt[j,i] > hp1_wb_srt_thr[i]) 1 else 0
      
      hp1_wb_srt_cn<-hp1_wb_srt_cn+s
      hp1_wb_srt_vt[j,i]<-hp1_wb_srt[j,i]-hp1_wb_srt_thr[i]  #its difference table used in clami+
      
    }
  }
  hp1_wb_srt[hp1_wb_srt_nRowsDf+1,i]<-(hp1_wb_srt_cy+hp1_wb_srt_cn)  #the last row is added as the no. of violation row for each metrics. correct
}

#Now arrange the metrics in ascending order based on the value of last row and choose first ceilling(log2(n)) metrics for next phase.
#Firstly remove the last coloumn as label_bug for ordering metrics because its factor coloumn
hp1_wb_srt2<-hp1_wb_srt[,1:ncol(hp1_wb_srt)-1]
hp1_wb_srt2=hp1_wb_srt2[,order(hp1_wb_srt2[nrow(hp1_wb_srt2),])]

# select the metrics with minimum mvs score
hp1_wb_srt_mvs_min<-min(hp1_wb_srt[hp1_wb_srt_nRowsDf+1,1:hp1_wb_srt_nColsDf])  #min of 8th row
#Select the metrics which has 8th row contain the hp1_wb_srt_mvs_min value and store it in hp1_wb_srt_mvs means dataset with select metrics with min mvs
hp1_wb_srt_nmvs<- sum(hp1_wb_srt_mvs_min == hp1_wb_srt[hp1_wb_srt_nRowsDf+1,1:hp1_wb_srt_nColsDf])  #number of minimum vs
cat("number of selected metrics=",hp1_wb_srt_nmvs)
hp1_wb_srt_mvs<-matrix()
hp1_wb_srt_mvs_row<-hp1_wb_srt[hp1_wb_srt_nRowsDf+1,]  #last row as MVS values of each metrics
hp1_wb_srt_mvs_row_order<-hp1_wb_srt_mvs_row[,order(hp1_wb_srt_mvs_row)]


hp1_wb_srt_wmvs<-hp1_wb_srt[1:hp1_wb_srt_nRowsDf,]   #Removing the last row means mvs row from hp1_wb_srt, without mvs row
hp1_wb_srt_thr_st<-c()  #st means selected threshold of selected metrics

#Select n% of metrics using hp1_wb_srt_select_nm


hp1_wb_srt_select_nm<-ceiling(log2(hp1_wb_srt_nColsDf))  #Select first ceiling(log2(n)) metrics based on mvs
hp1_wb_srt_ams<-matrix()   #hp1_wb_srt after metric selection
for(i in 1:hp1_wb_srt_select_nm)
{
  hp1_wb_srt_ams<-cbind(hp1_wb_srt_ams,hp1_wb_srt2[i])   #hp1_wb_srt2 used
}

hp1_wb_srt_ams<-cbind(hp1_wb_srt_ams[2:ncol(hp1_wb_srt_ams)],label_bug=hp1_wb_srt$label_bug)

str(hp1_wb_srt_ams)
hp1_wb_srt_ams_colnames<-colnames(hp1_wb_srt_ams[,1:hp1_wb_srt_select_nm])
hp1_wb_srt_mvs_row_order_first3<-hp1_wb_srt_mvs_row_order[,1:hp1_wb_srt_select_nm]
#write.csv(as.data.frame(hp1_wb_srt_mvs_row_order_first3),"hp1_wb_srt_ams_first3.csv")
#write.csv(hp1_wb_srt_mvs_row_order,'hp1_wb_srt_mvs_row_order.csv')
##=======================================================================================
##metric selection completed now work on instance selection on 18.00, 03/01/2021
hp1_wb_srt_ams_nColsDf<-ncol(hp1_wb_srt_ams)-1
hp1_wb_srt_ams_nRowsDf<-nrow(hp1_wb_srt_ams)-1

for(i in 1:hp1_wb_srt_ams_nRowsDf)
{
  hp1_DS_iy<-0  #no. of yes in coloumn cy
  hp1_DS_in<-0
  for(j in 1:hp1_wb_srt_ams_nColsDf)
  {
    if( hp1_wb_srt_ams$label_bug[i]=='YES')
    {
      r <- if(hp1_wb_srt_ams[i,j] <= hp1_wb_srt_ams[hp1_wb_srt_ams_nRowsDf+1,j]) 1 else 0
      hp1_DS_iy<-hp1_DS_iy+r
      #hp1_DS_vt[j,i]<-hp1_wb_srt_ams[j,i]-hp1_DS_thr_st[i]     #its for clami+ to write difference
    }
    else if( hp1_wb_srt_ams$label_bug[i]=='NO')
    {
      s<- if(hp1_wb_srt_ams[i,j] > hp1_wb_srt_ams[hp1_wb_srt_ams_nRowsDf+1,j]) 1 else 0
      
      hp1_DS_in<-hp1_DS_in+s
      #hp1_DS_ivt[j,i]<-hp1_wb_srt_ams[j,i]-hp1_DS_thr[i]
    }
  }
  hp1_wb_srt_ams[i,hp1_wb_srt_ams_nColsDf+2]<-(hp1_DS_iy+hp1_DS_in)  #the 4th coloumn is the no. of violation row for each instances correct done 13.1 6/1/2021
}
#View(hp1_wb_srt_ams)
count(hp1_wb_srt_ams[,ncol(hp1_wb_srt_ams)])  #instance violation score completed   data[,ncol(data)]
## UP to here it executed completely  its final ds after metric selection
################################################################################
#Now arrange the instances in ascending order based on the value of last coloumn and choose first n-n/12 metrics for next phase.
#Firstly remove the last coloumn as label_bug for ordering rows because its factor coloumn
hp1_wb_srt_ams_row_ord<-hp1_wb_srt_ams[order(hp1_wb_srt_ams[,ncol(hp1_wb_srt_ams)]),]
hp1_wb_srt_ams_row_ord_nrow<-nrow(hp1_wb_srt_ams)
hp1_wb_srt_ams_row_ord<-hp1_wb_srt_ams_row_ord[-hp1_wb_srt_ams_row_ord_nrow,]
str(hp1_wb_srt_ams_row_ord)

hp1_wb_srt_ams_cnv<-count(hp1_wb_srt_ams_row_ord$V5)   #count the number of violations existed in the last coloumn V5
hp1_wb_srt_ams_cnv


#Now select top n% of instances using hp1_const_sr_srtp


hp1_wb_srt_ams_row_ord<-head(hp1_wb_srt_ams_row_ord,hp1_const_sr_srtp*hp1_wb_srt_ams_row_ord_nrow)   #select no. of rows
hp1_wb_srt_ams_is_nc<-count(hp1_wb_srt_ams_row_ord$label_bug)  #count the number of yes and no
hp1_wb_srt_ams_is_nc
hp1_wb_srt_ams_is_nc_r<-hp1_wb_srt_ams_is_nc[1,2]/hp1_wb_srt_ams_is_nc[2,2]  #ratio
hp1_wb_srt_ams_is_nc_r

hp1_wb_srt_ams_row_ord <- hp1_wb_srt_ams_row_ord[sample(nrow(hp1_wb_srt_ams_row_ord)),]   #randomize instances
hp1_wb_srt_ams_ais<-hp1_wb_srt_ams_row_ord[,-ncol(hp1_wb_srt_ams_row_ord)]                  #Remove last coloumn as ivs

hp1_wb_srt_ams_ais_nrow_ratio<-cbind(nrow(hp1_wb_srt_ams_ais),hp1_wb_srt_ams_is_nc_r)
#write.csv(as.data.frame(hp1_wb_srt_ams_ais_nrow_ratio),"hp1_wb_srt_ams_ais_nrow_ratio.csv")
#View(hp1_wb_srt_ams_ais)  #This is the final input data
######################################################################

#Now apply the SMOTE: Synthetic Minority Oversampling Technique To Handle Class Imbalancy In Binary Classification
count(hp1_wb_srt_ams_ais$label_bug)
hp1_wb_srt_ams_ais_rose<-ROSE(label_bug~.,data=hp1_wb_srt_ams_ais,seed = 786)$data   #it is better than smote
count(hp1_wb_srt_ams_ais_rose$label_bug)
set.seed(200)
hp1_wb_srt_ams_ais_rose<-hp1_wb_srt_ams_ais_rose[sample(nrow(hp1_wb_srt_ams_ais_rose)),]

#=======================================================================================================================
# #data<-read.csv(file.choose())  *********change 1:#m,#m+2
hp1_wb_srt_data<-hp1_wb_srt_ams_ais_rose     #label name   "label_bug"
hp1_srtp_org<-hp1_org                      #label name "bug"
hp1_wb_srt_predictors<-colnames(hp1_wb_srt_data[,1:ncol(hp1_wb_srt_data)-1])
#write.csv(hp1_wb_srt_predictors,"hp1_wb_srt_predictors.csv")
hp1_srtp_org_predictors<-colnames(hp1_srtp_org[,1:ncol(hp1_srtp_org)-1])
hp1_wb_srt_outcomeName<-'label_bug'
hp1_srtp_org_outcomeName<-'bug'

fitControl <- trainControl(      #train control for TCLP
  method = "cv",
  number = 10,
  savePredictions = 'final',  summaryFunction = multiClassSummary,
  classProbs = T )

set.seed(3)
hp1_wb_srt_model_rf<-train(hp1_wb_srt_data[,hp1_wb_srt_predictors],hp1_wb_srt_data[,hp1_wb_srt_outcomeName],method='rf',trControl=fitControl,tuneLength=5)
hp1_wb_srt_model_rf[["resample"]]<-na.mean(hp1_wb_srt_model_rf[["resample"]])
hp1_wb_srt_model_rf
hp1_srtp_org$predrf_bug<-predict(object = hp1_wb_srt_model_rf,hp1_srtp_org[,hp1_srtp_org_predictors],type='raw')   #predict on original data using pre trained model
hp1_wb_srt_cm_rf<-caret::confusionMatrix(factor(hp1_srtp_org$bug),factor(hp1_srtp_org$predrf_bug), mode= "prec_recall")
hp1_wb_srt_cm_rf
#View(hp1_srtp_org$predrf_bug)

hp1_wb_srtp_01<-hp1_srtp_org       #store the actual data in temp data with labelin with 0 1 as bug
hp1_wb_srtp_01$predrf_bug<-ifelse(hp1_wb_srtp_01$predrf_bug == "YES", 1,0)
hp1_wb_srtp_mcc<-mccr(hp1_01b$bug,hp1_wb_srtp_01$predrf_bug)



#Silhouette Plot and AvgSilWidth
hp1_srtp_org$bnum<-ifelse(hp1_srtp_org$predrf_bug=='YES',2,1)
hp1_srtp_org_dist<-dist(hp1_srtp_org[,1:6],method = "euclidean")^2
hp1_srtp_org_sil=silhouette(hp1_srtp_org$bnum,hp1_srtp_org_dist)
set.seed(21)
pdf("hp1_srtp_org_sil_plot.pdf")
plot(hp1_srtp_org_sil)  # for eth)
dev.off()
hp1_srtp_org_sil_summary<-summary(hp1_srtp_org_sil, FUN=mean)


#Melting data according to predrf_bug
hp1_srtp_org_melt<- melt(hp1_srtp_org[,c(1,2,3,4,5,6,8)], id.var = "predrf_bug")
hp1_srtp_org_melt$predrf_bug<-as.factor(hp1_srtp_org_melt$predrf_bug)
pdf("hp1_srtp_org_bp_predrf_bug_wise.pdf")
p <- ggplot(data = hp1_srtp_org_melt, aes(x=variable, y=value)) + geom_boxplot(aes(fill=predrf_bug))
p + facet_wrap( ~ variable, scales="free")   #Boxplot as per requirment need to scale y axis
dev.off()
#######################################################################################################

#################################################################################

#start cube root and inverse


#===============================================================================



## Threshold calculation start and implement crt
#Check with calculate log without adding 1
hp1_wb_crt1<-hp1_wb_crt+1  #Add 1 to the metrics so that easily calculated the log of the metrics. because log 0 is not defined
hp1_wb_crt1[,c(1:ncol(hp1_wb_crt1))] <- (hp1_wb_crt1[,c(1:ncol(hp1_wb_crt1))])^(1/3)   #take the log of all metrics

hp1_wb_crt1_mean<- colMeans(hp1_wb_crt1)
hp1_wb_crt1_median<- apply(hp1_wb_crt1,2,median)  #2 is margin means col, sd is function
hp1_wb_crt1_mode<- apply(hp1_wb_crt1,2,getmode)  #2 is margin means col, sd is function
hp1_wb_crt1_sd<- apply(hp1_wb_crt1,2,sd)  #2 is margin means col, sd is function

#===============================================================================================================
# now calculate threshold of each  metric for the logrithmic value of metrics given metric values thr=exp(T'), T'=mean+sd
hp1_wb_crt1_th_m<-c()
for(i in 1:ncol(hp1_wb_crt1))
{
  hp1_wb_crt1_th_m[i]<-hp1_wb_crt1_mean[i]+hp1_wb_crt1_sd[i]
}
hp1_wb_crt1_thr_exp<-(hp1_wb_crt1_th_m)^3 #Actual threshold used for crt methods
hp1_wb_crt1_threshold<-hp1_wb_crt1_thr_exp  #Actual threshold
#=================Derived threshold completed====================================================
#Now for better ment mean median mode should be aomost equal to the threshold
hp1_wb_crt_basic_stats<-basicStats(hp1_wb_crt)
hp1_wb_crt1_basic_stats<-basicStats(hp1_wb_crt1)
hp1_wb_crt_mean<- colMeans(hp1_wb_crt)
hp1_wb_crt_median<- apply(hp1_wb_crt,2,median)  #2 is margin means col, sd is function
hp1_wb_crt_mode<- apply(hp1_wb_crt,2,getmode)
hp1_wb_crt_skewness<-matrix(c(hp1_wb_crt_basic_stats[15,],hp1_wb_crt1_basic_stats[15,],hp1_wb_crt_mean,hp1_wb_crt_median,hp1_wb_crt_mode,hp1_wb_crt1_threshold),ncol = ncol(hp1_wb_crt1), byrow = TRUE)
rownames(hp1_wb_crt_skewness)<- c('Skewness_before_log','Skewness_after_log','mean','median','mode','threshold')
#write.csv(hp1_wb_crt_skewness,"hp1_wb_crt_skewness_bfrLog_aftrLog_mean_med_mode_thrs.csv")
#=====================================================================================
#Find the value of k for each instances and and label as YES or NO
hp1_wb_crt_nRowsDf<-nrow(hp1_wb_crt)
hp1_wb_crt_nColsDf<-ncol(hp1_wb_crt)
hp1_wb_crt_thr<- hp1_wb_crt1_threshold  #use threshold calculated using Log transform methods

#Now calculate the values of k for each instances
for(i in 1:hp1_wb_crt_nRowsDf)
{
  loop=0
  for(j in 1:hp1_wb_crt_nColsDf)
  {
    if(hp1_wb_crt[i,j]>=hp1_wb_crt_thr[j])
    {
      loop=loop+1
    }
    hp1_wb_crt$k[i]<-loop
  }
}

#Lable top half as YES and bottom half as NO, find max of k values then divide by 2 and get floor value.
for(i in 1:hp1_wb_crt_nRowsDf)
{
  if(hp1_wb_crt$k[i]>floor(max(hp1_wb_crt$k)/2))
  {
    hp1_wb_crt$label_bug[i]='YES'
  }
  else {
    hp1_wb_crt$label_bug[i]='NO'
  }
}
#K value calculation and Labeling completed Now compute performance of TCL
hp1_wb_crt_cm<-confusionMatrix(factor(hp1_org$bug), factor(hp1_wb_crt$label_bug), mode="prec_recall")
#View(hp1_wb_crt$label_bug)
hp1_wb_crt_cm

hp1_wb_crt_01<-hp1_wb_crt       #store the actial data in temp data with labelin with 0 1 as bug
hp1_wb_crt_01$label_bug<-ifelse(hp1_wb_crt_01$label_bug == "YES", 1,0)
hp1_wb_crt_mcc<-mccr(hp1_01b$bug,hp1_wb_crt_01$label_bug)
hp1_wb_crt_mcc



#Silhouette Plot and AvgSilWidth
hp1_wb_crt$bnum<-ifelse(hp1_wb_crt$label_bug=='YES',2,1)
hp1_wb_crt_dist<-dist(hp1_wb_crt[,1:6],method = "euclidean")^2
hp1_wb_crt_sil=silhouette(hp1_wb_crt$bnum,hp1_wb_crt_dist)
set.seed(21)
pdf("hp1_wb_crt_sil_plot.pdf")
plot(hp1_wb_crt_sil)  #0.78 for eth)
dev.off()
hp1_wb_crt_sil_summary<-summary(hp1_wb_crt_sil, FUN=mean)


#Melting data according to label_bug
hp1_wb_crt_melt<- melt(hp1_wb_crt[,c(1,2,3,4,5,6,8)], id.var = "label_bug")
hp1_wb_crt_melt$label_bug<-as.factor(hp1_wb_crt_melt$label_bug)
pdf("hp1_wb_crt_bp_label_bug_wise.pdf")
p <- ggplot(data = hp1_wb_crt_melt, aes(x=variable, y=value)) + geom_boxplot(aes(fill=label_bug))
p + facet_wrap( ~ variable, scales="free")   #Boxplot as per requirment need to scale y axis
dev.off()

#============== TCL completed   ========================================


#Start TCLP
##========   Metric selection  started  ==============================================================
#Now perform Metric selection based on metric violation score MVS (metric violation score), vt (violation table)
#MVS of i metric=Ci/Fi,  Ci is no. of violation in ith metiric and Fi is the no. of metric values in the i-th metric(no. of instances)
#After labeling the instances now perform metric selection
hp1_wb_crt_nRowsDf<-nrow(hp1_wb_crtp)
hp1_wb_crt_nColsDf<-ncol(hp1_wb_crtp)
hp1_wb_crt_vt<-matrix(nrow=hp1_wb_crt_nRowsDf,ncol=hp1_wb_crt_nColsDf)
for(i in 1:hp1_wb_crt_nColsDf)
{
  hp1_wb_crt_cy<-0  #no. of yes in coloumn cy
  hp1_wb_crt_cn<-0
  for(j in 1:hp1_wb_crt_nRowsDf)
  {
    if( hp1_wb_crt$label_bug[j]=='YES')
    {
      r <- if(hp1_wb_crt[j,i] <= hp1_wb_crt_thr[i]) 1 else 0
      hp1_wb_crt_cy<-hp1_wb_crt_cy+r
      hp1_wb_crt_vt[j,i]<-hp1_wb_crt[j,i]-hp1_wb_crt_thr[i]     #its for clami+ to write difference
    }
    else if( hp1_wb_crt$label_bug[j]=='NO')
    {
      s<- if(hp1_wb_crt[j,i] > hp1_wb_crt_thr[i]) 1 else 0
      
      hp1_wb_crt_cn<-hp1_wb_crt_cn+s
      hp1_wb_crt_vt[j,i]<-hp1_wb_crt[j,i]-hp1_wb_crt_thr[i]  #its difference table used in clami+
      
    }
  }
  hp1_wb_crt[hp1_wb_crt_nRowsDf+1,i]<-(hp1_wb_crt_cy+hp1_wb_crt_cn)  #the last row is added as the no. of violation row for each metrics. correct
}

#Now arrange the metrics in ascending order based on the value of last row and choose first ceilling(log2(n)) metrics for next phase.
#Firstly remove the last coloumn as label_bug for ordering metrics because its factor coloumn
hp1_wb_crt2<-hp1_wb_crt[,1:ncol(hp1_wb_crt)-1]
hp1_wb_crt2=hp1_wb_crt2[,order(hp1_wb_crt2[nrow(hp1_wb_crt2),])]

# select the metrics with minimum mvs score
hp1_wb_crt_mvs_min<-min(hp1_wb_crt[hp1_wb_crt_nRowsDf+1,1:hp1_wb_crt_nColsDf])  #min of 8th row
#Select the metrics which has 8th row contain the hp1_wb_crt_mvs_min value and store it in hp1_wb_crt_mvs means dataset with select metrics with min mvs
hp1_wb_crt_nmvs<- sum(hp1_wb_crt_mvs_min == hp1_wb_crt[hp1_wb_crt_nRowsDf+1,1:hp1_wb_crt_nColsDf])  #number of minimum vs
cat("number of selected metrics=",hp1_wb_crt_nmvs)
hp1_wb_crt_mvs<-matrix()
hp1_wb_crt_mvs_row<-hp1_wb_crt[hp1_wb_crt_nRowsDf+1,]  #last row as MVS values of each metrics
hp1_wb_crt_mvs_row_order<-hp1_wb_crt_mvs_row[,order(hp1_wb_crt_mvs_row)]


hp1_wb_crt_wmvs<-hp1_wb_crt[1:hp1_wb_crt_nRowsDf,]   #Removing the last row means mvs row from hp1_wb_crt, without mvs row
hp1_wb_crt_thr_st<-c()  #st means selected threshold of selected metrics

#Select n% of metrics using hp1_wb_crt_select_nm


hp1_wb_crt_select_nm<-ceiling(log2(hp1_wb_crt_nColsDf))  #Select first ceiling(log2(n)) metrics based on mvs
hp1_wb_crt_ams<-matrix()   #hp1_wb_crt after metric selection
for(i in 1:hp1_wb_crt_select_nm)
{
  hp1_wb_crt_ams<-cbind(hp1_wb_crt_ams,hp1_wb_crt2[i])   #hp1_wb_crt2 used
}

hp1_wb_crt_ams<-cbind(hp1_wb_crt_ams[2:ncol(hp1_wb_crt_ams)],label_bug=hp1_wb_crt$label_bug)

str(hp1_wb_crt_ams)
hp1_wb_crt_ams_colnames<-colnames(hp1_wb_crt_ams[,1:hp1_wb_crt_select_nm])
hp1_wb_crt_mvs_row_order_first3<-hp1_wb_crt_mvs_row_order[,1:hp1_wb_crt_select_nm]
#write.csv(as.data.frame(hp1_wb_crt_mvs_row_order_first3),"hp1_wb_crt_ams_first3.csv")
#write.csv(hp1_wb_crt_mvs_row_order,'hp1_wb_crt_mvs_row_order.csv')
##=======================================================================================
##metric selection completed now work on instance selection on 18.00, 03/01/2021
hp1_wb_crt_ams_nColsDf<-ncol(hp1_wb_crt_ams)-1
hp1_wb_crt_ams_nRowsDf<-nrow(hp1_wb_crt_ams)-1

for(i in 1:hp1_wb_crt_ams_nRowsDf)
{
  hp1_DS_iy<-0  #no. of yes in coloumn cy
  hp1_DS_in<-0
  for(j in 1:hp1_wb_crt_ams_nColsDf)
  {
    if( hp1_wb_crt_ams$label_bug[i]=='YES')
    {
      r <- if(hp1_wb_crt_ams[i,j] <= hp1_wb_crt_ams[hp1_wb_crt_ams_nRowsDf+1,j]) 1 else 0
      hp1_DS_iy<-hp1_DS_iy+r
      #hp1_DS_vt[j,i]<-hp1_wb_crt_ams[j,i]-hp1_DS_thr_st[i]     #its for clami+ to write difference
    }
    else if( hp1_wb_crt_ams$label_bug[i]=='NO')
    {
      s<- if(hp1_wb_crt_ams[i,j] > hp1_wb_crt_ams[hp1_wb_crt_ams_nRowsDf+1,j]) 1 else 0
      
      hp1_DS_in<-hp1_DS_in+s
      #hp1_DS_ivt[j,i]<-hp1_wb_crt_ams[j,i]-hp1_DS_thr[i]
    }
  }
  hp1_wb_crt_ams[i,hp1_wb_crt_ams_nColsDf+2]<-(hp1_DS_iy+hp1_DS_in)  #the 4th coloumn is the no. of violation row for each instances correct done 13.1 6/1/2021
}
#View(hp1_wb_crt_ams)
count(hp1_wb_crt_ams[,ncol(hp1_wb_crt_ams)])  #instance violation score completed   data[,ncol(data)]
## UP to here it executed completely   its final ds after metric selection
################################################################################
#Now arrange the instances in ascending order based on the value of last coloumn and choose first n-n/12 metrics for next phase.
#Firstly remove the last coloumn as label_bug for ordering rows because its factor coloumn
hp1_wb_crt_ams_row_ord<-hp1_wb_crt_ams[order(hp1_wb_crt_ams[,ncol(hp1_wb_crt_ams)]),]
hp1_wb_crt_ams_row_ord_nrow<-nrow(hp1_wb_crt_ams)
hp1_wb_crt_ams_row_ord<-hp1_wb_crt_ams_row_ord[-hp1_wb_crt_ams_row_ord_nrow,]
str(hp1_wb_crt_ams_row_ord)

hp1_wb_crt_ams_cnv<-count(hp1_wb_crt_ams_row_ord$V5)   #count the number of violations existed in the last coloumn V5
hp1_wb_crt_ams_cnv


#Now select top n% of instances using hp1_const_sr_crtp


hp1_wb_crt_ams_row_ord<-head(hp1_wb_crt_ams_row_ord,hp1_const_sr_crtp*hp1_wb_crt_ams_row_ord_nrow)   #select no. of rows
hp1_wb_crt_ams_is_nc<-count(hp1_wb_crt_ams_row_ord$label_bug)  #count the number of yes and no
hp1_wb_crt_ams_is_nc
hp1_wb_crt_ams_is_nc_r<-hp1_wb_crt_ams_is_nc[1,2]/hp1_wb_crt_ams_is_nc[2,2]  #ratio
hp1_wb_crt_ams_is_nc_r

hp1_wb_crt_ams_row_ord <- hp1_wb_crt_ams_row_ord[sample(nrow(hp1_wb_crt_ams_row_ord)),]   #randomize instances
hp1_wb_crt_ams_ais<-hp1_wb_crt_ams_row_ord[,-ncol(hp1_wb_crt_ams_row_ord)]                  #Remove last coloumn as ivs

hp1_wb_crt_ams_ais_nrow_ratio<-cbind(nrow(hp1_wb_crt_ams_ais),hp1_wb_crt_ams_is_nc_r)
#write.csv(as.data.frame(hp1_wb_crt_ams_ais_nrow_ratio),"hp1_wb_crt_ams_ais_nrow_ratio.csv")
#View(hp1_wb_crt_ams_ais)  #This is the final input data
######################################################################

#Now apply the SMOTE: Synthetic Minority Oversampling Technique To Handle Class Imbalancy In Binary Classification
count(hp1_wb_crt_ams_ais$label_bug)
hp1_wb_crt_ams_ais_rose<-ROSE(label_bug~.,data=hp1_wb_crt_ams_ais,seed = 786)$data   #it is better than smote
count(hp1_wb_crt_ams_ais_rose$label_bug)
set.seed(200)
hp1_wb_crt_ams_ais_rose<-hp1_wb_crt_ams_ais_rose[sample(nrow(hp1_wb_crt_ams_ais_rose)),]

#=======================================================================================================================
# #data<-read.csv(file.choose())  *********change 1:#m,#m+2
hp1_wb_crt_data<-hp1_wb_crt_ams_ais_rose     #label name   "label_bug"
hp1_crtp_org<-hp1_org                      #label name "bug"
hp1_wb_crt_predictors<-colnames(hp1_wb_crt_data[,1:ncol(hp1_wb_crt_data)-1])
#write.csv(hp1_wb_crt_predictors,"hp1_wb_crt_predictors.csv")
hp1_crtp_org_predictors<-colnames(hp1_crtp_org[,1:ncol(hp1_crtp_org)-1])
hp1_wb_crt_outcomeName<-'label_bug'
hp1_crtp_org_outcomeName<-'bug'

fitControl <- trainControl(      #train control for TCLP
  method = "cv",
  number = 10,
  savePredictions = 'final',  summaryFunction = multiClassSummary,
  classProbs = T )

set.seed(3)
hp1_wb_crt_model_rf<-train(hp1_wb_crt_data[,hp1_wb_crt_predictors],hp1_wb_crt_data[,hp1_wb_crt_outcomeName],method='rf',trControl=fitControl,tuneLength=5)
hp1_wb_crt_model_rf[["resample"]]<-na.mean(hp1_wb_crt_model_rf[["resample"]])
hp1_wb_crt_model_rf
hp1_crtp_org$predrf_bug<-predict(object = hp1_wb_crt_model_rf,hp1_crtp_org[,hp1_crtp_org_predictors],type='raw')   #predict on original data using pre trained model
hp1_wb_crt_cm_rf<-caret::confusionMatrix(factor(hp1_crtp_org$bug),factor(hp1_crtp_org$predrf_bug), mode= "prec_recall")
hp1_wb_crt_cm_rf
#View(hp1_crtp_org$predrf_bug)

hp1_wb_crtp_01<-hp1_crtp_org       #store the actual data in temp data with labelin with 0 1 as bug
hp1_wb_crtp_01$predrf_bug<-ifelse(hp1_wb_crtp_01$predrf_bug == "YES", 1,0)
hp1_wb_crtp_mcc<-mccr(hp1_01b$bug,hp1_wb_crtp_01$predrf_bug)



#Silhouette Plot and AvgSilWidth
hp1_crtp_org$bnum<-ifelse(hp1_crtp_org$predrf_bug=='YES',2,1)
hp1_crtp_org_dist<-dist(hp1_crtp_org[,1:6],method = "euclidean")^2
hp1_crtp_org_sil=silhouette(hp1_crtp_org$bnum,hp1_crtp_org_dist)
set.seed(21)
pdf("hp1_crtp_org_sil_plot.pdf")
plot(hp1_crtp_org_sil)  #0.78 for eth)
dev.off()
hp1_crtp_org_sil_summary<-summary(hp1_crtp_org_sil, FUN=mean)


#Melting data according to predrf_bug
hp1_crtp_org_melt<- melt(hp1_crtp_org[,c(1,2,3,4,5,6,8)], id.var = "predrf_bug")
hp1_crtp_org_melt$predrf_bug<-as.factor(hp1_crtp_org_melt$predrf_bug)
pdf("hp1_crtp_org_bp_predrf_bug_wise.pdf")
p <- ggplot(data = hp1_crtp_org_melt, aes(x=variable, y=value)) + geom_boxplot(aes(fill=predrf_bug))
p + facet_wrap( ~ variable, scales="free")   #Boxplot as per requirment need to scale y axis
dev.off()
#######################################################################################################

#Start simple majority voting ensemble (TCL,SRT,CRT) =TSCP
hp1_comb_predictedYN<-cbind(hp1_org, TCL=factor(hp1_wb_tcl[1:hp1_org_nrow,]$label_bug),TCLP=hp1_tclp_org$predrf_bug,SRT=factor(hp1_wb_srt[1:hp1_org_nrow,]$label_bug), SRTP=hp1_srtp_org$predrf_bug, CRT=factor(hp1_wb_crt[1:hp1_org_nrow,]$label_bug),CRTP=hp1_crtp_org$predrf_bug)



hp1_comb_predictedYN$hp1_pred_TSC<-ifelse(hp1_comb_predictedYN$TCL=='YES' & hp1_comb_predictedYN$CRT=='YES','YES',ifelse(hp1_comb_predictedYN$TCL=='YES' & hp1_comb_predictedYN$SRT=='YES','YES',ifelse(hp1_comb_predictedYN$CRT=='YES' & hp1_comb_predictedYN$SRT=='YES','YES','NO')))
hp1_comb_predictedYN$hp1_pred_TSC<- factor(hp1_comb_predictedYN$hp1_pred_TSC, levels = c("NO","YES"));
hp1_cm_TSC<-caret::confusionMatrix(factor(hp1_comb_predictedYN$bug),factor(hp1_comb_predictedYN$hp1_pred_TSC), mode= "prec_recall")
hp1_cm_TSC #STK3   0.66

hp1_comb_predictedYN$hp1_pred_TSC_01<-hp1_comb_predictedYN$hp1_pred_TSC       #store the actial data in temp data with labelin with 0 1 as Bug
hp1_comb_predictedYN$hp1_pred_TSC_01<-ifelse(hp1_comb_predictedYN$hp1_pred_TSC_01 == "YES", 1,0)
hp1_mcc_TSC<-mccr(hp1_01b$bug,hp1_comb_predictedYN$hp1_pred_TSC_01)
hp1_mcc_TSC


#The majority vote (unweighted average) tclp, acl, clami
hp1_comb_predictedYN$hp1_pred_TSCP<-ifelse(hp1_comb_predictedYN$TCLP=='YES' & hp1_comb_predictedYN$CRTP=='YES','YES',ifelse(hp1_comb_predictedYN$TCLP=='YES' & hp1_comb_predictedYN$SRTP=='YES','YES',ifelse(hp1_comb_predictedYN$CRTP=='YES' & hp1_comb_predictedYN$SRTP=='YES','YES','NO')))
hp1_comb_predictedYN$hp1_pred_TSCP<- factor(hp1_comb_predictedYN$hp1_pred_TSCP, levels = c("NO","YES"));
hp1_cm_TSCP<-caret::confusionMatrix(factor(hp1_comb_predictedYN$bug),factor(hp1_comb_predictedYN$hp1_pred_TSCP), mode= "prec_recall")
hp1_cm_TSCP #STK3   0.66

hp1_comb_predictedYN$hp1_pred_TSCP_01<-hp1_comb_predictedYN$hp1_pred_TSCP       #store the actial data in temp data with labelin with 0 1 as Bug
hp1_comb_predictedYN$hp1_pred_TSCP_01<-ifelse(hp1_comb_predictedYN$hp1_pred_TSCP_01 == "YES", 1,0)
hp1_mcc_TSCP<-mccr(hp1_01b$bug,hp1_comb_predictedYN$hp1_pred_TSCP_01)
hp1_mcc_TSCP
#==========================================================================================================================

hp1_acc<-matrix(c(hp1_wb_tcl_cm[["overall"]][["Accuracy"]]*100,hp1_wb_tcl_cm_rf[["overall"]][["Accuracy"]]*100,hp1_wb_srt_cm[["overall"]][["Accuracy"]]*100,hp1_wb_srt_cm_rf[["overall"]][["Accuracy"]]*100,hp1_wb_crt_cm[["overall"]][["Accuracy"]]*100,hp1_wb_crt_cm_rf[["overall"]][["Accuracy"]]*100,hp1_cm_TSC[["overall"]][["Accuracy"]]*100,hp1_cm_TSCP[["overall"]][["Accuracy"]]*100),ncol = 8, byrow = TRUE)
colnames(hp1_acc)<-c("TCL","TCLP","SRT","SRTP","CRT","CRTP","TSC","TSCP")
rownames(hp1_acc)<-c("hp1_acc")
hp1_acc<-as.table(hp1_acc)
hp1_acc

hp1_fm<-matrix(c(hp1_wb_tcl_cm[["byClass"]][["F1"]],hp1_wb_tcl_cm_rf[["byClass"]][["F1"]],hp1_wb_srt_cm[["byClass"]][["F1"]],hp1_wb_srt_cm_rf[["byClass"]][["F1"]],hp1_wb_crt_cm[["byClass"]][["F1"]],hp1_wb_crt_cm_rf[["byClass"]][["F1"]],hp1_cm_TSC[["byClass"]][["F1"]],hp1_cm_TSCP[["byClass"]][["F1"]]),ncol = 8, byrow = TRUE)
colnames(hp1_fm)<-c("TCL","TCLP","SRT","SRTP","CRT","CRTP","TSC","TSCP")
rownames(hp1_acc)<-c("hp1_fm")
hp1_fm<-as.table(hp1_fm)
hp1_fm

hp1_mcc<-matrix(c(hp1_wb_tcl_mcc,hp1_wb_tclp_mcc,hp1_wb_srt_mcc,hp1_wb_srtp_mcc,hp1_wb_crt_mcc,hp1_wb_crtp_mcc,hp1_mcc_TSCP,hp1_mcc_TSCP),ncol = 8, byrow = TRUE)
colnames(hp1_mcc)<-c("TCL","TCLP","SRT","SRTP","CRT","CRTP","TSC","TSCP")
rownames(hp1_mcc)<-c("hp1_mcc")
hp1_mcc<-as.table(hp1_mcc)
hp1_mcc

hp1_cmb_table<-matrix(c(hp1_org_ncol,hp1_org_nrow,hp1_org_bugpr,hp1_acc,hp1_fm,hp1_mcc),ncol = 27,byrow=FALSE)
colnames(hp1_cmb_table)<-c('NCOL','NROW','bug%',"TCL","TCLP","SRT","SRTP","CRT","CRTP","TSC","TSCP","TCL","TCLP","SRT","SRTP","CRT","CRTP","TSC","TSCP","TCL","TCLP","SRT","SRTP","CRT","CRTP","TSC","TSCP")
rownames(hp1_cmb_table) <- c("hp1_DS")
hp1_cmb_table<-as.table(hp1_cmb_table)
write.csv(hp1_cmb_table,"hp1_cmb_table_TSCP.csv")




 
