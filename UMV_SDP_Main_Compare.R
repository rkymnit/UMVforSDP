
set.seed(786)
packageslist <- list("tidyverse", "cluster", "factoextra", "e1071", "caret", "gridExtra", "fBasics", "dplyr", "plyr", "ggplot2", "reshape2", "VGAM", "corrplot","mlbench", "MASS", "RColorBrewer", "RANN", "ROCR", "pROC", "MLmetrics", "imputeTS", "xgboost","mccr"
                     ,"kohonen","vegan","mclust","cclust","FCPS","DMwR","ROSE")
load_packages <- lapply(packageslist, require, character.only = T)
#Random oversampling example : its generate artificaial synthetic data on sampling methods and smoothed bootstrap approach.
getmode<-function(x){
  uniq<-unique(x)
  uniq[which.max(tabulate(match(x,uniq)))]
}

ONS=0.99     #Om Namah Shivay
etherium_const_sr_tclp<-ONS
etherium_const_sr_clami<-ONS
etherium_const_sr_aclp<-ONS
etherium_const_sr_clmp<-ONS    #Means actual clami+

etherium_org<-read.csv("etherium_defect.csv")
str(etherium_org)
#etherium_org<-etherium_org[,1:6]
 #etherium_org$bug<-as.factor(etherium_org$bug>0)
#levels(etherium_org$bug) <- c("NO", "YES")

etherium_org_ncol<-ncol(etherium_org)
etherium_org_nrow<-nrow(etherium_org)
etherium_org_count<-count(etherium_org$bug)
etherium_org_bugpr<-(etherium_org_count[2,2]/etherium_org_nrow)*100   #bug % module in dataset

etherium_01b<-etherium_org
levels(etherium_01b$bug)<-0:1  #level bug as 0 and 1 for No and YES
str(etherium_01b)
# remove the label bug from etherium_org, it become etherium_wb
etherium_wb<-etherium_org[,-ncol(etherium_org)]

etherium_wb_cla<-etherium_wb
etherium_wb_clami<-etherium_wb
etherium_wb_tcl<-etherium_wb
etherium_wb_tclp<-etherium_wb
etherium_wb_km<-etherium_wb
etherium_sl_cv<-etherium_org
# 
####################################
## Threshold calculation start and implement tcl
#Check with calculate log without adding 1
etherium_wb_tcl1<-etherium_wb_tcl+1  #Add 1 to the metrics so that easily calculated the log of the metrics. because log 0 is not defined
etherium_wb_tcl1[,c(1:ncol(etherium_wb_tcl1))] <- log(etherium_wb_tcl1[,c(1:ncol(etherium_wb_tcl1))])   #take the log of all metrics

etherium_wb_tcl1_mean<- colMeans(etherium_wb_tcl1)
etherium_wb_tcl1_median<- apply(etherium_wb_tcl1,2,median)  #2 is margin means col, sd is function
etherium_wb_tcl1_mode<- apply(etherium_wb_tcl1,2,getmode)  #2 is margin means col, sd is function
etherium_wb_tcl1_sd<- apply(etherium_wb_tcl1,2,sd)  #2 is margin means col, sd is function

#===============================================================================================================
# now calculate threshold of each  metric for the logrithmic value of metrics given metric values thr=exp(T'), T'=mean+sd
etherium_wb_tcl1_th_m<-c()
for(i in 1:ncol(etherium_wb_tcl1))
{
  etherium_wb_tcl1_th_m[i]<-etherium_wb_tcl1_mean[i]+etherium_wb_tcl1_sd[i]
}
etherium_wb_tcl1_thr_exp<-exp(etherium_wb_tcl1_th_m) #Actual threshold used for tcl methods
etherium_wb_tcl1_threshold<-etherium_wb_tcl1_thr_exp  #Actual threshold
#=================Derived threshold completed====================================================
#Now for better ment mean median mode should be aomost equal to the threshold
etherium_wb_tcl_basic_stats<-basicStats(etherium_wb_tcl)
etherium_wb_tcl1_basic_stats<-basicStats(etherium_wb_tcl1)
etherium_wb_tcl_mean<- colMeans(etherium_wb_tcl)
etherium_wb_tcl_median<- apply(etherium_wb_tcl,2,median)  #2 is margin means col, sd is function
etherium_wb_tcl_mode<- apply(etherium_wb_tcl,2,getmode)
etherium_wb_tcl_skewness<-matrix(c(etherium_wb_tcl_basic_stats[15,],etherium_wb_tcl1_basic_stats[15,],etherium_wb_tcl_mean,etherium_wb_tcl_median,etherium_wb_tcl_mode,etherium_wb_tcl1_threshold),ncol = ncol(etherium_wb_tcl1), byrow = TRUE)
rownames(etherium_wb_tcl_skewness)<- c('Skewness_before_log','Skewness_after_log','mean','median','mode','threshold')
#write.csv(etherium_wb_tcl_skewness,"etherium_wb_tcl_skewness_bfrLog_aftrLog_mean_med_mode_thrs.csv")
#=====================================================================================
#Find the value of k for each instances and and label as YES or NO
etherium_wb_tcl_nRowsDf<-nrow(etherium_wb_tcl)
etherium_wb_tcl_nColsDf<-ncol(etherium_wb_tcl)
etherium_wb_tcl_thr<- etherium_wb_tcl1_threshold  #use threshold calculated using Log transform methods

#Now calculate the values of k for each instances
for(i in 1:etherium_wb_tcl_nRowsDf)
{
  loop=0
  for(j in 1:etherium_wb_tcl_nColsDf)
  {
    if(etherium_wb_tcl[i,j]>=etherium_wb_tcl_thr[j])
    {
      loop=loop+1
    }
    etherium_wb_tcl$k[i]<-loop
  }
}

#Lable top half as YES and bottom half as NO, find max of k values then divide by 2 and get floor value.
for(i in 1:etherium_wb_tcl_nRowsDf)
{
  if(etherium_wb_tcl$k[i]>floor(max(etherium_wb_tcl$k)/2))
  {
    etherium_wb_tcl$label_bug[i]='YES'
  }
  else {
    etherium_wb_tcl$label_bug[i]='NO'
  }
}
#K value calculation and Labeling completed Now compute performance of TCL
etherium_wb_tcl_cm<-confusionMatrix(factor(etherium_org$bug), factor(etherium_wb_tcl$label_bug), mode="prec_recall")
#View(etherium_wb_tcl$label_bug)

etherium_wb_tcl_01<-etherium_wb_tcl       #store the actial data in temp data with labelin with 0 1 as bug
etherium_wb_tcl_01$label_bug<-ifelse(etherium_wb_tcl_01$label_bug == "YES", 1,0)
etherium_wb_tcl_mcc<-mccr(etherium_01b$bug,etherium_wb_tcl_01$label_bug)

#============== TCL completed   ========================================
#Start TCLP
##========   Metric selection  started  ==============================================================
#Now perform Metric selection based on metric violation score MVS (metric violation score), vt (violation table)
#MVS of i metric=Ci/Fi,  Ci is no. of violation in ith metiric and Fi is the no. of metric values in the i-th metric(no. of instances)
#After labeling the instances now perform metric selection
etherium_wb_tcl_nRowsDf<-nrow(etherium_wb_tclp)
etherium_wb_tcl_nColsDf<-ncol(etherium_wb_tclp)
etherium_wb_tcl_vt<-matrix(nrow=etherium_wb_tcl_nRowsDf,ncol=etherium_wb_tcl_nColsDf)
for(i in 1:etherium_wb_tcl_nColsDf)
{
  etherium_wb_tcl_cy<-0  #no. of yes in coloumn cy
  etherium_wb_tcl_cn<-0
  for(j in 1:etherium_wb_tcl_nRowsDf)
  {
    if( etherium_wb_tcl$label_bug[j]=='YES')
    {
      r <- if(etherium_wb_tcl[j,i] <= etherium_wb_tcl_thr[i]) 1 else 0
      etherium_wb_tcl_cy<-etherium_wb_tcl_cy+r
      etherium_wb_tcl_vt[j,i]<-etherium_wb_tcl[j,i]-etherium_wb_tcl_thr[i]     #its for clami+ to write difference
    }
    else if( etherium_wb_tcl$label_bug[j]=='NO')
    {
      s<- if(etherium_wb_tcl[j,i] > etherium_wb_tcl_thr[i]) 1 else 0
      
      etherium_wb_tcl_cn<-etherium_wb_tcl_cn+s
      etherium_wb_tcl_vt[j,i]<-etherium_wb_tcl[j,i]-etherium_wb_tcl_thr[i]  #its difference table used in clami+
      
    }
  }
  etherium_wb_tcl[etherium_wb_tcl_nRowsDf+1,i]<-(etherium_wb_tcl_cy+etherium_wb_tcl_cn)  #the last row is added as the no. of violation row for each metrics. correct
}

#Now arrange the metrics in ascending order based on the value of last row and choose first ceilling(log2(n)) metrics for next phase.
#Firstly remove the last coloumn as label_bug for ordering metrics because its factor coloumn
etherium_wb_tcl2<-etherium_wb_tcl[,1:ncol(etherium_wb_tcl)-1]
etherium_wb_tcl2=etherium_wb_tcl2[,order(etherium_wb_tcl2[nrow(etherium_wb_tcl2),])]

# select the metrics with minimum mvs score
etherium_wb_tcl_mvs_min<-min(etherium_wb_tcl[etherium_wb_tcl_nRowsDf+1,1:etherium_wb_tcl_nColsDf])  #min of 8th row
#Select the metrics which has 8th row contain the etherium_wb_tcl_mvs_min value and store it in etherium_wb_tcl_mvs means dataset with select metrics with min mvs
etherium_wb_tcl_nmvs<- sum(etherium_wb_tcl_mvs_min == etherium_wb_tcl[etherium_wb_tcl_nRowsDf+1,1:etherium_wb_tcl_nColsDf])  #number of minimum vs
cat("number of selected metrics=",etherium_wb_tcl_nmvs)
etherium_wb_tcl_mvs<-matrix()
etherium_wb_tcl_mvs_row<-etherium_wb_tcl[etherium_wb_tcl_nRowsDf+1,]  #last row as MVS values of each metrics
etherium_wb_tcl_mvs_row_order<-etherium_wb_tcl_mvs_row[,order(etherium_wb_tcl_mvs_row)]


etherium_wb_tcl_wmvs<-etherium_wb_tcl[1:etherium_wb_tcl_nRowsDf,]   #Removing the last row means mvs row from etherium_wb_tcl, without mvs row
etherium_wb_tcl_thr_st<-c()  #st means selected threshold of selected metrics

etherium_wb_tcl_select_nm<-ceiling(log2(etherium_wb_tcl_nColsDf))  #Select first ceiling(log2(n)) metrics based on mvs
etherium_wb_tcl_ams<-matrix()   #etherium_wb_tcl after metric selection
for(i in 1:etherium_wb_tcl_select_nm)
{
  etherium_wb_tcl_ams<-cbind(etherium_wb_tcl_ams,etherium_wb_tcl2[i])   #etherium_wb_tcl2 used
}

etherium_wb_tcl_ams<-cbind(etherium_wb_tcl_ams[2:ncol(etherium_wb_tcl_ams)],label_bug=etherium_wb_tcl$label_bug)

str(etherium_wb_tcl_ams)
etherium_wb_tcl_ams_colnames<-colnames(etherium_wb_tcl_ams[,1:etherium_wb_tcl_select_nm])
etherium_wb_tcl_mvs_row_order_first5<-etherium_wb_tcl_mvs_row_order[,1:etherium_wb_tcl_select_nm]
write.csv(as.data.frame(etherium_wb_tcl_mvs_row_order_first5),"etherium_wb_tcl_ams_first5.csv")
write.csv(etherium_wb_tcl_mvs_row_order,'etherium_wb_tcl_mvs_row_order.csv')
##=======================================================================================
##metric selection completed now work on instance selection on 18.00, 03/01/2021
etherium_wb_tcl_ams_nColsDf<-ncol(etherium_wb_tcl_ams)-1
etherium_wb_tcl_ams_nRowsDf<-nrow(etherium_wb_tcl_ams)-1

for(i in 1:etherium_wb_tcl_ams_nRowsDf)
{
  etherium_DS_iy<-0  #no. of yes in coloumn cy
  etherium_DS_in<-0
  for(j in 1:etherium_wb_tcl_ams_nColsDf)
  {
    if( etherium_wb_tcl_ams$label_bug[i]=='YES')
    {
      r <- if(etherium_wb_tcl_ams[i,j] <= etherium_wb_tcl_ams[etherium_wb_tcl_ams_nRowsDf+1,j]) 1 else 0
      etherium_DS_iy<-etherium_DS_iy+r
      #etherium_DS_vt[j,i]<-etherium_wb_tcl_ams[j,i]-etherium_DS_thr_st[i]     #its for clami+ to write difference
    }
    else if( etherium_wb_tcl_ams$label_bug[i]=='NO')
    {
      s<- if(etherium_wb_tcl_ams[i,j] > etherium_wb_tcl_ams[etherium_wb_tcl_ams_nRowsDf+1,j]) 1 else 0
      
      etherium_DS_in<-etherium_DS_in+s
      #etherium_DS_ivt[j,i]<-etherium_wb_tcl_ams[j,i]-etherium_DS_thr[i]
    }
  }
  etherium_wb_tcl_ams[i,etherium_wb_tcl_ams_nColsDf+2]<-(etherium_DS_iy+etherium_DS_in)  #the 4th coloumn is the no. of violation row for each instances correct done 13.1 6/1/2021
}
#View(etherium_wb_tcl_ams)
count(etherium_wb_tcl_ams[,ncol(etherium_wb_tcl_ams)])  #instance violation score completed   data[,ncol(data)]
## UP to here it executed completely 20.15  03/01/2021  its final ds after metric selection
################################################################################
#Now arrange the instances in ascending order based on the value of last coloumn and choose first n-n/12 metrics for next phase.
#Firstly remove the last coloumn as label_bug for ordering rows because its factor coloumn
etherium_wb_tcl_ams_row_ord<-etherium_wb_tcl_ams[order(etherium_wb_tcl_ams[,ncol(etherium_wb_tcl_ams)]),]
etherium_wb_tcl_ams_row_ord_nrow<-nrow(etherium_wb_tcl_ams)
etherium_wb_tcl_ams_row_ord<-etherium_wb_tcl_ams_row_ord[-etherium_wb_tcl_ams_row_ord_nrow,]
str(etherium_wb_tcl_ams_row_ord)

etherium_wb_tcl_ams_cnv<-count(etherium_wb_tcl_ams_row_ord$V7)   #count the number of violations existed in the last coloumn V7
etherium_wb_tcl_ams_cnv
etherium_wb_tcl_ams_row_ord<-head(etherium_wb_tcl_ams_row_ord,etherium_const_sr_tclp*etherium_wb_tcl_ams_row_ord_nrow)   #select no. of rows
etherium_wb_tcl_ams_is_nc<-count(etherium_wb_tcl_ams_row_ord$label_bug)  #count the number of yes and no
etherium_wb_tcl_ams_is_nc
etherium_wb_tcl_ams_is_nc_r<-etherium_wb_tcl_ams_is_nc[1,2]/etherium_wb_tcl_ams_is_nc[2,2]  #ratio
etherium_wb_tcl_ams_is_nc_r

etherium_wb_tcl_ams_row_ord <- etherium_wb_tcl_ams_row_ord[sample(nrow(etherium_wb_tcl_ams_row_ord)),]   #randomize instances
etherium_wb_tcl_ams_ais<-etherium_wb_tcl_ams_row_ord[,-ncol(etherium_wb_tcl_ams_row_ord)]                  #Remove last coloumn as ivs

etherium_wb_tcl_ams_ais_nrow_ratio<-cbind(nrow(etherium_wb_tcl_ams_ais),etherium_wb_tcl_ams_is_nc_r)
write.csv(as.data.frame(etherium_wb_tcl_ams_ais_nrow_ratio),"etherium_wb_tcl_ams_ais_nrow_ratio.csv")
#View(etherium_wb_tcl_ams_ais)  #This is the final input data
######################################################################

#Now apply the SMOTE: Synthetic Minority Oversampling Technique To Handle Class Imbalancy In Binary Classification
count(etherium_wb_tcl_ams_ais$label_bug)
etherium_wb_tcl_ams_ais_rose<-ROSE(label_bug~.,data=etherium_wb_tcl_ams_ais,seed = 786)$data   #it is better than smote
count(etherium_wb_tcl_ams_ais_rose$label_bug)
set.seed(200)
etherium_wb_tcl_ams_ais_rose<-etherium_wb_tcl_ams_ais_rose[sample(nrow(etherium_wb_tcl_ams_ais_rose)),]

#=======================================================================================================================
# #data<-read.csv(file.choose())  *********change 1:#m,#m+2
etherium_wb_tcl_data<-etherium_wb_tcl_ams_ais_rose     #label name   "label_bug"
etherium_tclp_org<-etherium_org                      #label name "bug"
etherium_wb_tcl_predictors<-colnames(etherium_wb_tcl_data[,1:ncol(etherium_wb_tcl_data)-1])
#write.csv(etherium_wb_tcl_predictors,"etherium_wb_tcl_predictors.csv")
etherium_tclp_org_predictors<-colnames(etherium_tclp_org[,1:ncol(etherium_tclp_org)-1])
etherium_wb_tcl_outcomeName<-'label_bug'
etherium_tclp_org_outcomeName<-'bug'

fitControl <- trainControl(      #train control for TCLP
  method = "cv",
  number = 10,
  savePredictions = 'final',  summaryFunction = multiClassSummary,
  classProbs = T )

set.seed(3)
etherium_wb_tcl_model_rf<-train(etherium_wb_tcl_data[,etherium_wb_tcl_predictors],etherium_wb_tcl_data[,etherium_wb_tcl_outcomeName],method='rf',trControl=fitControl,tuneLength=5)
etherium_wb_tcl_model_rf[["resample"]]<-na.mean(etherium_wb_tcl_model_rf[["resample"]])
etherium_wb_tcl_model_rf
etherium_tclp_org$predrf_bug<-predict(object = etherium_wb_tcl_model_rf,etherium_tclp_org[,etherium_tclp_org_predictors],type='raw')   #predict on original data using pre trained model
etherium_wb_tcl_cm_rf<-caret::confusionMatrix(factor(etherium_tclp_org$bug),factor(etherium_tclp_org$predrf_bug), mode= "prec_recall")
etherium_wb_tcl_cm_rf
#View(etherium_tclp_org$predrf_bug)

etherium_wb_tclp_01<-etherium_tclp_org       #store the actual data in temp data with labelin with 0 1 as bug
etherium_wb_tclp_01$predrf_bug<-ifelse(etherium_wb_tclp_01$predrf_bug == "YES", 1,0)
etherium_wb_tclp_mcc<-mccr(etherium_01b$bug,etherium_wb_tclp_01$predrf_bug)

############################################################################################################
#============================================================================================
#ACL Started
#step1: calculate half average of each metric (HAF)

# remove the label bug from etherium_org, it become etherium_wb
#etherium_wb<-etherium_org[,-ncol(etherium_org)]
etherium_wb_acl<-etherium_wb
etherium_wb_acl_mean<-colMeans(etherium_wb_acl)
etherium_wb_acl_mean_half<-etherium_wb_acl_mean/2     #(HAM) Half average of metric
etherium_wb_acl_nRowsDf<-nrow(etherium_wb_acl)
etherium_wb_acl_nColsDf<-ncol(etherium_wb_acl)
#etherium_wb_acl_thr<-etherium_wb_acl_mean_half
etherium_wb_acl_thr<-etherium_wb_tcl_thr
#Step2: Create metric violation matrix based on HAF, means for a metric, if the value of instace is greater than the HAF then set 1 else 0.
#Now calculate the values of k or rsum or MIVS for each instances
for(i in 1:etherium_wb_acl_nRowsDf)
{
  loop=0
  for(j in 1:etherium_wb_acl_nColsDf)
  {
    if(etherium_wb_acl[i,j]>=etherium_wb_acl_thr[j])
    {
      loop=loop+1
    }
    etherium_wb_acl$rsum[i]<-loop
  }
}

#MIVS (metric of instance violation score) coloumn calculated of ACL method. (rsum)
#HAM: Half average of metrics (etherium_wb_acl_mean_half or etherium_wb_acl_thr)
#MVM: metric violation metrics  (etherium_wb_acl_vt)
#Algo 1 completed exceot 4th point
#etherium_wb_acl_vt$rsum<-apply(etherium_wb_acl_vt,1,sum)/etherium_wb_acl_nColsDf
#Step4: find the cutoff of rsum according to ACL to decide which instances have higher mivs or buggy
#step 4 is calculating the cutoff and labeling instances
etherium_wb_acl_hms<-etherium_wb_acl_nColsDf/2   #half of metrics set (hms)
etherium_wb_acl_ni=0     # #number of instances which mivs(rsum) is greater than hms (pd:)

for(i in 1:etherium_wb_acl_nRowsDf)
{
  if(etherium_wb_acl$rsum[i]>=etherium_wb_acl_hms)   #use >=
  {
    etherium_wb_acl_ni<-etherium_wb_acl_ni+1
    #etherium_wb_acl_vt$label_bug[i]='YES'
  }
  # else {
  #   etherium_wb_acl_vt$label_bug[i]='NO'
  # }
}

#Now  calculate possible rate of defect (PDr) based on pd or etherium_wb_acl_ni
etherium_wb_acl_PDr<-etherium_wb_acl_ni/etherium_wb_acl_nRowsDf   #PDr=PD/NO. of instances

#AMIVS+:Average of Metrics of Instance Violation Score (using only distinct MIVS)
#AMIVS+:average of distinct MIVS 
etherium_wb_acl_dmivs<-unique(etherium_wb_acl$rsum)   #distinct values of mivs
etherium_wb_acl_admivs<-mean(etherium_wb_acl_dmivs)    #AMIVS+
#AMIVS: is the average of MIVS
etherium_wb_acl_amivs<-mean(etherium_wb_acl$rsum)
#MVM:metric violation matrix                                                    
#MMIVS: median of MIVS(rsum) 
etherium_wb_acl_mmivs<-median(etherium_wb_acl$rsum)   

#HMIVS:Harmonic mean of AMIVS+ (admivs) and MMIVS
etherium_wb_acl_hmivs<-(2*etherium_wb_acl_admivs*etherium_wb_acl_mmivs)/(etherium_wb_acl_admivs+etherium_wb_acl_mmivs)

#If the possible rate of defect is less than 0.5, it is meaning the real rate of defect is low. And the labeling cutoff 
#should be high, then the labeling cutoff is calculated by equation 2. I
#if PDr is less than 0.5 then eqn2 else eqn1 (from ppr algo 2, 3rd point)
{
  if(etherium_wb_acl_PDr<0.5)
  {
    etherium_wb_acl_cutoff<-(2*etherium_wb_acl_nColsDf*(1-etherium_wb_acl_PDr)*etherium_wb_acl_hmivs)/((etherium_wb_acl_nColsDf*(1-etherium_wb_acl_PDr))+etherium_wb_acl_hmivs)
  }
  else {
    etherium_wb_acl_cutoff<-(etherium_wb_acl_hmivs*etherium_wb_acl_PDr)+((etherium_wb_acl_nColsDf-etherium_wb_acl_amivs)*(1-etherium_wb_acl_PDr))    #Eqn1
  }
}


#Now if the MiVS of instances is greater than labeling cutoff, the instances is labeled as faulty, otherwise non-faulty. rsum >etherium_wb_acl_cutoff
for(i in 1:etherium_wb_acl_nRowsDf)
{
  if(etherium_wb_acl$rsum[i] > etherium_wb_acl_cutoff)
  {
    etherium_wb_acl$label_bug[i]='YES'
  }
  else {
    etherium_wb_acl$label_bug[i]='NO'
  }
}

etherium_wb_acl_cm<-confusionMatrix(factor(etherium_org$bug), factor(etherium_wb_acl$label_bug), mode="prec_recall")
#View(etherium_wb_acl$label_bug)
etherium_wb_acl_01<-etherium_wb_acl       #store the actial data in temp data with labelin with 0 1 as bug
etherium_wb_acl_01$label_bug<-ifelse(etherium_wb_acl_01$label_bug == "YES", 1,0)
etherium_wb_acl_mcc<-mccr(etherium_01b$bug,etherium_wb_acl_01$label_bug)

#Successful first phase of ACL as clustering and labeling


#============================================================================================
#Now we can perform ACLP means acl plus metric selection, instance selection then ML
########################   Metric selection  started  ###############################################################################################
#Now perform Metric selection based on metric violation score MVS (metric violation score), vt (violation table)
#MVS of i metric=Ci/Fi,  Ci is no. of violation in ith metiric and Fi is the no. of metric values in the i-th metric(no. of instances)
#After labeling the instances now perform metric selection
#etherium_wb_aclp<-etherium_org[,-ncol(etherium_org)]
etherium_wb_acl_nRowsDf<-nrow(etherium_wb)
etherium_wb_acl_nColsDf<-ncol(etherium_wb)
etherium_wb_acl_vt<-matrix(nrow=etherium_wb_acl_nRowsDf,ncol=etherium_wb_acl_nColsDf)
for(i in 1:etherium_wb_acl_nColsDf)
{
  etherium_wb_acl_cy<-0  #no. of yes in coloumn cy
  etherium_wb_acl_cn<-0
  for(j in 1:etherium_wb_acl_nRowsDf)
  {
    if( etherium_wb_acl$label_bug[j]=='YES')
    {
      r <- if(etherium_wb_acl[j,i] <= etherium_wb_acl_thr[i]) 1 else 0
      etherium_wb_acl_cy<-etherium_wb_acl_cy+r
      etherium_wb_acl_vt[j,i]<-etherium_wb_acl[j,i]-etherium_wb_acl_thr[i]     #its for clami+ to write difference
    }
    else if( etherium_wb_acl$label_bug[j]=='NO')
    {
      s<- if(etherium_wb_acl[j,i] > etherium_wb_acl_thr[i]) 1 else 0

      etherium_wb_acl_cn<-etherium_wb_acl_cn+s
      etherium_wb_acl_vt[j,i]<-etherium_wb_acl[j,i]-etherium_wb_acl_thr[i]  #its difference table used in clami+

    }
  }
  etherium_wb_acl[etherium_wb_acl_nRowsDf+1,i]<-(etherium_wb_acl_cy+etherium_wb_acl_cn)  #the last row is added as the no. of violation row for each metrics. correct
}

#Now arrange the metrics in ascending order based on the value of last row and choose first ceilling(log2(n)) metrics for next phase.
#Firstly remove the last coloumn as label_bug for ordering metrics because its factor coloumn
etherium_wb_acl2<-etherium_wb_acl[,1:ncol(etherium_wb_acl)-1]
etherium_wb_acl2=etherium_wb_acl2[,order(etherium_wb_acl2[nrow(etherium_wb_acl2),])]

# select the metrics with minimum mvs score
etherium_wb_acl_mvs_min<-min(etherium_wb_acl[etherium_wb_acl_nRowsDf+1,1:etherium_wb_acl_nColsDf])  #min of 8th row
#Select the metrics which has 8th row contain the etherium_wb_acl_mvs_min value and store it in etherium_wb_acl_mvs means dataset with select metrics with min mvs
etherium_wb_acl_nmvs<- sum(etherium_wb_acl_mvs_min == etherium_wb_acl[etherium_wb_acl_nRowsDf+1,1:etherium_wb_acl_nColsDf])  #number of minimum vs
cat("number of selected metrics=",etherium_wb_acl_nmvs)
etherium_wb_acl_mvs<-matrix()
etherium_wb_acl_wmvs<-etherium_wb_acl[1:etherium_wb_acl_nRowsDf,]   #Removing the last row means mvs row from etherium_wb_acl, without mvs row
etherium_wb_acl_thr_st<-c()  #st means selected threshold of selected metrics

etherium_wb_acl_select_nm<-ceiling(log2(etherium_wb_acl_nColsDf))  #Select first ceiling(log2(n)) metrics based on mvs
etherium_wb_acl_ams<-matrix()   #etherium_wb_acl after metric selection
for(i in 1:etherium_wb_acl_select_nm)
{
  etherium_wb_acl_ams<-cbind(etherium_wb_acl_ams,etherium_wb_acl2[i])   #etherium_wb_acl2 used
}

etherium_wb_acl_ams<-cbind(etherium_wb_acl_ams[2:ncol(etherium_wb_acl_ams)],label_bug=etherium_wb_acl$label_bug)

str(etherium_wb_acl_ams)

#####################################################################################################
##metric selection completed now work on instance selection on 18.00, 03/01/2021
etherium_wb_acl_ams_nColsDf<-ncol(etherium_wb_acl_ams)-1
etherium_wb_acl_ams_nRowsDf<-nrow(etherium_wb_acl_ams)-1

for(i in 1:etherium_wb_acl_ams_nRowsDf)
{
  etherium_DS_iy<-0  #no. of yes in coloumn cy
  etherium_DS_in<-0
  for(j in 1:etherium_wb_acl_ams_nColsDf)
  {
    if( etherium_wb_acl_ams$label_bug[i]=='YES')
    {
      r <- if(etherium_wb_acl_ams[i,j] <= etherium_wb_acl_ams[etherium_wb_acl_ams_nRowsDf+1,j]) 1 else 0
      etherium_DS_iy<-etherium_DS_iy+r
      #etherium_DS_vt[j,i]<-etherium_wb_acl_ams[j,i]-etherium_DS_thr_st[i]     #its for clami+ to write difference
    }
    else if( etherium_wb_acl_ams$label_bug[i]=='NO')
    {
      s<- if(etherium_wb_acl_ams[i,j] > etherium_wb_acl_ams[etherium_wb_acl_ams_nRowsDf+1,j]) 1 else 0

      etherium_DS_in<-etherium_DS_in+s
      #etherium_DS_ivt[j,i]<-etherium_wb_acl_ams[j,i]-etherium_DS_thr[i]
    }
  }
  etherium_wb_acl_ams[i,etherium_wb_acl_ams_nColsDf+2]<-(etherium_DS_iy+etherium_DS_in)  #the 4th coloumn is the no. of violation row for each instances correct done 13.1 6/1/2021
}
#View(etherium_wb_acl_ams)
count(etherium_wb_acl_ams[,ncol(etherium_wb_acl_ams)])  #instance violation score completed   data[,ncol(data)]
## UP to here it executed completely  its final ds after metric selection
################################################################################
#Now arrange the instances in ascending order based on the value of last coloumn and choose first n-n/12 metrics for next phase.
#Firstly remove the last coloumn as label_bug for ordering rows because its factor coloumn
etherium_wb_acl_ams_row_ord<-etherium_wb_acl_ams[order(etherium_wb_acl_ams[,ncol(etherium_wb_acl_ams)]),]
etherium_wb_acl_ams_row_ord_nrow<-nrow(etherium_wb_acl_ams)
etherium_wb_acl_ams_row_ord<-etherium_wb_acl_ams_row_ord[-etherium_wb_acl_ams_row_ord_nrow,]
str(etherium_wb_acl_ams_row_ord)

etherium_wb_acl_ams_cnv<-count(etherium_wb_acl_ams_row_ord$V7)   #count the number of violations existed in the last coloumn V7
etherium_wb_acl_ams_cnv
etherium_wb_acl_ams_row_ord<-head(etherium_wb_acl_ams_row_ord,etherium_const_sr_aclp*etherium_wb_acl_ams_row_ord_nrow)
etherium_wb_acl_ams_is_nc<-count(etherium_wb_acl_ams_row_ord$label_bug)  #count the number of yes and no
etherium_wb_acl_ams_is_nc
etherium_wb_acl_ams_is_nc_r<-etherium_wb_acl_ams_is_nc[1,2]/etherium_wb_acl_ams_is_nc[2,2]  #ratio
etherium_wb_acl_ams_is_nc_r

etherium_wb_acl_ams_row_ord <- etherium_wb_acl_ams_row_ord[sample(nrow(etherium_wb_acl_ams_row_ord)),]   #randomize instances
etherium_wb_acl_ams_ais<-etherium_wb_acl_ams_row_ord[,-ncol(etherium_wb_acl_ams_row_ord)]                  #Remove last coloumn as ivs

#View(etherium_wb_acl_ams_ais)  #This is the final input data
######################################################################

#Now apply the SMOTE: Synthetic Minority Oversampling Technique To Handle Class Imbalancy In Binary Classification
count(etherium_wb_acl_ams_ais$label_bug)
etherium_wb_acl_ams_ais_rose<-ROSE(label_bug~.,data=etherium_wb_acl_ams_ais,seed = 786)$data   #it is better than smote
count(etherium_wb_acl_ams_ais_rose$label_bug)
set.seed(200)
etherium_wb_acl_ams_ais_rose<-etherium_wb_acl_ams_ais_rose[sample(nrow(etherium_wb_acl_ams_ais_rose)),]



##############################################################################
# #data<-read.csv(file.choose())  *********change 1:#m,#m+2
etherium_wb_acl_data<-etherium_wb_acl_ams_ais_rose     #label name   "label_bug"
etherium_aclp_org<-etherium_org                      #label name "bug"
etherium_wb_acl_predictors<-colnames(etherium_wb_acl_data[,1:ncol(etherium_wb_acl_data)-1])
#write.csv(etherium_wb_acl_predictors,"etherium_wb_acl_predictors.csv")
etherium_aclp_org_predictors<-colnames(etherium_aclp_org[,1:ncol(etherium_aclp_org)-1])
etherium_wb_acl_outcomeName<-'label_bug'
etherium_aclp_org_outcomeName<-'bug'

fitControl <- trainControl(      #train control for TCLP
  method = "cv",
  number = 10,
  savePredictions = 'final',  summaryFunction = multiClassSummary,
  classProbs = T )

set.seed(3)
etherium_wb_acl_model_rf<-train(etherium_wb_acl_data[,etherium_wb_acl_predictors],etherium_wb_acl_data[,etherium_wb_acl_outcomeName],method='rf',trControl=fitControl,tuneLength=5)
etherium_wb_acl_model_rf[["resample"]]<-na.mean(etherium_wb_acl_model_rf[["resample"]])
etherium_wb_acl_model_rf
etherium_aclp_org$predrf_bug<-predict(object = etherium_wb_acl_model_rf,etherium_aclp_org[,etherium_aclp_org_predictors])   #predict on original data using pre trained model
etherium_wb_acl_cm_rf<-caret::confusionMatrix(factor(etherium_aclp_org$bug),factor(etherium_aclp_org$predrf_bug), mode= "prec_recall")
etherium_wb_acl_cm_rf

#View(etherium_aclp_org$predrf_bug)

etherium_wb_aclp_01<-etherium_aclp_org       #store the actual data in temp data with labelin with 0 1 as bug
etherium_wb_aclp_01$predrf_bug<-ifelse(etherium_wb_aclp_01$predrf_bug == "YES", 1,0)
etherium_wb_aclp_mcc<-mccr(etherium_01b$bug,etherium_wb_aclp_01$predrf_bug)

#==============================================================================================================================

############################################################################################################
#CLA Methods started
etherium_wb_cla_median<- apply(etherium_wb_cla,2,median)
etherium_wb_cla_nRowsDf<-nrow(etherium_wb_cla)
etherium_wb_cla_nColsDf<-ncol(etherium_wb_cla)
etherium_wb_cla_thr<- etherium_wb_cla_median  #use threshold calculated using Log transform methods

#Now calculate the values of k for each instances
for(i in 1:etherium_wb_cla_nRowsDf)
{
  loop=0
  for(j in 1:etherium_wb_cla_nColsDf)
  {
    if(etherium_wb_cla[i,j]>=etherium_wb_cla_thr[j])
    {
      loop=loop+1
    }
    etherium_wb_cla$k[i]<-loop
  }
}

#Lable top half as YES and bottom half as NO, find max of k values then divide by 2 and get floor value.
for(i in 1:etherium_wb_cla_nRowsDf)
{
  if(etherium_wb_cla$k[i]>floor(max(etherium_wb_cla$k)/2))
  {
    etherium_wb_cla$label_bug[i]='YES'
  }
  else {
    etherium_wb_cla$label_bug[i]='NO'
  }
}
#View(etherium_wb_tcl_01)  #K value calculation and Labeling completed Now compute performance of TCL
etherium_wb_cla_cm<-confusionMatrix(factor(etherium_org$bug), factor(etherium_wb_cla$label_bug), mode="prec_recall")

#View(etherium_wb_cla$label_bug)

etherium_wb_cla_01<-etherium_wb_cla       #store the actial data in temp data with labelin with 0 1 as bug
etherium_wb_cla_01$label_bug<-ifelse(etherium_wb_cla_01$label_bug == "YES", 1,0)
etherium_wb_cla_mcc<-mccr(etherium_01b$bug,etherium_wb_cla_01$label_bug)

##############################################################################################################
#Start CLAMI
########################   Metric selection  started  ###############################################################################################
#Now perform Metric selection based on metric violation score MVS (metric violation score), vt (violation table)
#MVS of i metric=Ci/Fi,  Ci is no. of violation in ith metiric and Fi is the no. of metric values in the i-th metric(no. of instances)
#After labeling the instances now perform metric selection
etherium_wb_cla_nRowsDf<-nrow(etherium_wb_clami)
etherium_wb_cla_nColsDf<-ncol(etherium_wb_clami)
etherium_wb_cla_thr<-etherium_wb_cla_median
etherium_wb_cla_vt<-matrix(nrow=etherium_wb_cla_nRowsDf,ncol=etherium_wb_cla_nColsDf)
for(i in 1:etherium_wb_cla_nColsDf)
{
  etherium_wb_cla_cy<-0  #no. of yes in coloumn cy
  etherium_wb_cla_cn<-0
  for(j in 1:etherium_wb_cla_nRowsDf)
  {
    if( etherium_wb_cla$label_bug[j]=='YES')
    { 
      r <- if(etherium_wb_cla[j,i] <= etherium_wb_cla_thr[i]) 1 else 0
      etherium_wb_cla_cy<-etherium_wb_cla_cy+r
      etherium_wb_cla_vt[j,i]<-etherium_wb_cla[j,i]-etherium_wb_cla_thr[i]     #its for clami+ to write difference  
    }
    else if( etherium_wb_cla$label_bug[j]=='NO')
    {
      s<- if(etherium_wb_cla[j,i] > etherium_wb_cla_thr[i]) 1 else 0
      
      etherium_wb_cla_cn<-etherium_wb_cla_cn+s
      etherium_wb_cla_vt[j,i]<-etherium_wb_cla[j,i]-etherium_wb_cla_thr[i]  #its difference table used in clami+
      
    }
  }
  etherium_wb_cla[etherium_wb_cla_nRowsDf+1,i]<-(etherium_wb_cla_cy+etherium_wb_cla_cn)  #the last row is added as the no. of violation row for each metrics. correct
}

#Now arrange the metrics in ascending order based on the value of last row and choose first ceilling(log2(n)) metrics for next phase.
#Firstly remove the last coloumn as label_bug for ordering metrics because its factor coloumn
etherium_wb_cla2<-etherium_wb_cla[,1:ncol(etherium_wb_cla)-1]
etherium_wb_cla2=etherium_wb_cla2[,order(etherium_wb_cla2[nrow(etherium_wb_cla2),])]

# select the metrics with minimum mvs score
etherium_wb_cla_mvs_min<-min(etherium_wb_cla[etherium_wb_cla_nRowsDf+1,1:etherium_wb_cla_nColsDf])  #min of 8th row
#Select the metrics which has 8th row contain the etherium_wb_cla_mvs_min value and store it in etherium_wb_cla_mvs means dataset with select metrics with min mvs
etherium_wb_cla_nmvs<- sum(etherium_wb_cla_mvs_min == etherium_wb_cla[etherium_wb_cla_nRowsDf+1,1:etherium_wb_cla_nColsDf])  #number of minimum vs
cat("number of selected metrics=",etherium_wb_cla_nmvs)
etherium_wb_cla_mvs<-matrix()
etherium_wb_cla_wmvs<-etherium_wb_cla[1:etherium_wb_cla_nRowsDf,]   #Removing the last row means mvs row from etherium_wb_cla, without mvs row
etherium_wb_cla_thr_st<-c()  #st means selected threshold of selected metrics

etherium_wb_cla_select_nm<-ceiling(log2(etherium_wb_cla_nColsDf))  #Select first ceiling(log2(n)) metrics based on mvs
etherium_wb_cla_ams<-matrix()   #etherium_wb_cla after metric selection
for(i in 1:etherium_wb_cla_select_nm)
{
  etherium_wb_cla_ams<-cbind(etherium_wb_cla_ams,etherium_wb_cla2[i])   #etherium_wb_cla2 used
}

etherium_wb_cla_ams<-cbind(etherium_wb_cla_ams[2:ncol(etherium_wb_cla_ams)],label_bug=etherium_wb_cla$label_bug)  

#str(etherium_wb_cla_ams)
#View(etherium_wb_cla_ams)
#####################################################################################################
##metric selection completed now work on instance selection on 18.00, 03/01/2021
etherium_wb_cla_ams_nColsDf<-ncol(etherium_wb_cla_ams)-1
etherium_wb_cla_ams_nRowsDf<-nrow(etherium_wb_cla_ams)-1

for(i in 1:etherium_wb_cla_ams_nRowsDf)
{
  etherium_DS_iy<-0  #no. of yes in coloumn cy
  etherium_DS_in<-0
  for(j in 1:etherium_wb_cla_ams_nColsDf)
  {
    if( etherium_wb_cla_ams$label_bug[i]=='YES')
    { 
      r <- if(etherium_wb_cla_ams[i,j] <= etherium_wb_cla_ams[etherium_wb_cla_ams_nRowsDf+1,j]) 1 else 0
      etherium_DS_iy<-etherium_DS_iy+r
      #etherium_DS_vt[j,i]<-etherium_wb_cla_ams[j,i]-etherium_DS_thr_st[i]     #its for clami+ to write difference  
    }
    else if( etherium_wb_cla_ams$label_bug[i]=='NO')
    {
      s<- if(etherium_wb_cla_ams[i,j] > etherium_wb_cla_ams[etherium_wb_cla_ams_nRowsDf+1,j]) 1 else 0
      
      etherium_DS_in<-etherium_DS_in+s
      #etherium_DS_ivt[j,i]<-etherium_wb_cla_ams[j,i]-etherium_DS_thr[i]
    }
  }
  etherium_wb_cla_ams[i,etherium_wb_cla_ams_nColsDf+2]<-(etherium_DS_iy+etherium_DS_in)  #the 4th coloumn is the no. of violation row for each instances correct done 13.1 6/1/2021 
}
#View(etherium_wb_cla_ams)
count(etherium_wb_cla_ams[,ncol(etherium_wb_cla_ams)])   #instance violation score completed [,ncol(etherium_wb_tcl_ams)]
## UP to here it executed completely 20.15  03/01/2021  its final ds after metric selection
################################################################################
#Now arrange the instances in ascending order based on the value of last coloumn and choose first n-n/12 metrics for next phase.
#Firstly remove the last coloumn as label_bug for ordering rows because its factor coloumn
etherium_wb_cla_ams_row_ord<-etherium_wb_cla_ams[order(etherium_wb_cla_ams[,ncol(etherium_wb_cla_ams)]),]
etherium_wb_cla_ams_row_ord_nrow<-nrow(etherium_wb_cla_ams)
etherium_wb_cla_ams_row_ord<-etherium_wb_cla_ams_row_ord[-etherium_wb_cla_ams_row_ord_nrow,]

etherium_wb_cla_ams_row_ord<-head(etherium_wb_cla_ams_row_ord,etherium_const_sr_clami*etherium_wb_cla_ams_row_ord_nrow)
count(etherium_wb_cla_ams_row_ord$label_bug)
count(etherium_wb_cla_ams_row_ord$V7)
etherium_wb_cla_ams_row_ord <- etherium_wb_cla_ams_row_ord[sample(nrow(etherium_wb_cla_ams_row_ord)),]   #randomize instances
etherium_wb_cla_ams_ais<-etherium_wb_cla_ams_row_ord[,-ncol(etherium_wb_cla_ams_row_ord)]                  #Remove last coloumn as ivs

etherium_wb_cla_ams_is_nc<-count(etherium_wb_cla_ams_row_ord$label_bug)  #count the number of yes and no
etherium_wb_cla_ams_is_nc
etherium_wb_cla_ams_is_nc_r<-etherium_wb_cla_ams_is_nc[1,2]/etherium_wb_cla_ams_is_nc[2,2]  #ratio
etherium_wb_cla_ams_is_nc_r
#View(etherium_wb_cla_ams_ais)  #This is the final input data for clami ml
######################################################################


# #data<-read.csv(file.choose())  *********change 1:#m,#m+2
etherium_wb_cla_data<-etherium_wb_cla_ams_ais     #label name   "label_bug"
etherium_clami_org<-etherium_org                      #label name "bug"
etherium_wb_cla_predictors<-colnames(etherium_wb_cla_data[,1:ncol(etherium_wb_cla_data)-1])

etherium_clami_org_predictors<-colnames(etherium_clami_org[,1:ncol(etherium_clami_org)-1])
etherium_wb_cla_outcomeName<-'label_bug'
etherium_clami_org_outcomeName<-'bug'

fitControl <- trainControl(
  method = "cv", 
  number = 10,  
  savePredictions = 'final',  summaryFunction = multiClassSummary,
  classProbs = T )

set.seed(3)
etherium_wb_cla_model_rf<-train(etherium_wb_cla_data[,etherium_wb_cla_predictors],etherium_wb_cla_data[,etherium_wb_cla_outcomeName],method='rf',trControl=fitControl,tuneLength=5)
etherium_wb_cla_model_rf[["resample"]]<-na.mean(etherium_wb_cla_model_rf[["resample"]])
etherium_wb_cla_model_rf

etherium_clami_org$predrf_bug<-predict(object = etherium_wb_cla_model_rf,etherium_clami_org[,etherium_clami_org_predictors],type='raw')   #predict on original data using pre trained model
etherium_wb_cla_cm_rf<-caret::confusionMatrix(factor(etherium_clami_org$bug),factor(etherium_clami_org$predrf_bug), mode= "prec_recall")              
etherium_wb_cla_cm_rf 

#View(etherium_clami_org$predrf_bug)

etherium_wb_clami_01<-etherium_clami_org       #store the actual data in temp data with labelin with 0 1 as bug
etherium_wb_clami_01$predrf_bug<-ifelse(etherium_wb_clami_01$predrf_bug == "YES", 1,0)
etherium_wb_clami_mcc<-mccr(etherium_01b$bug,etherium_wb_clami_01$predrf_bug)
etherium_wb_clami_mcc
#=====================================================================================================================================================

#Clami+ methods started

etherium_wb_clm<-etherium_org[,-ncol(etherium_org)]
etherium_wb_clm_median<- apply(etherium_wb_clm,2,median)
etherium_wb_clm_nRowsDf<-nrow(etherium_wb_clm)
etherium_wb_clm_nColsDf<-ncol(etherium_wb_clm)
etherium_wb_clm_thr<- etherium_wb_clm_median    #use threshold

etherium_wb_clm_vt<-matrix(nrow=nrow(etherium_wb_clm),ncol=ncol(etherium_wb_clm)) #empty violation table
for(i in 1:ncol(etherium_wb_clm))
{
  for(j in 1:nrow(etherium_wb_clm))
  {
    #its for clami+ to write difference  
    etherium_wb_clm_vt[j,i]<- 1/(1+exp(-(etherium_wb_clm[j,i]-etherium_wb_clm_thr[i])))
  }
}
colnames(etherium_wb_clm_vt)<-colnames(etherium_wb_clm)
#########################################################################################################
#View(etherium_wb_clm_vt)   #Correct
#etherium_wb_clm_vt$rsum<-rowSums(etherium_wb_clm_vt[,col_no])
etherium_wb_clm_vt<-as.data.frame(etherium_wb_clm_vt)
etherium_wb_clm_vt$rsum<-apply(etherium_wb_clm_vt,1,sum)/etherium_wb_clm_nColsDf

#Lable top half as YES and bottom half as NO, based on rsum >=0.5
for(i in 1:etherium_wb_clm_nRowsDf)
{
  if(etherium_wb_clm_vt$rsum[i]>=0.5)
  {
    etherium_wb_clm_vt$label_bug[i]='YES'
  }
  else {
    etherium_wb_clm_vt$label_bug[i]='NO'
  }
}
str(etherium_wb_clm_vt)
etherium_wb_clm<-etherium_wb_clm_vt
#Successful first phase of CLAmi+ as clustering and labeling
#View(etherium_wb_clm_vt)  #K value calculation and Labeling completed Now compute performance of clm
etherium_wb_clm_cm<-caret::confusionMatrix(factor(etherium_org$bug),factor(etherium_wb_clm$label_bug), mode="prec_recall")
#View(etherium_wb_clm$label_bug)

count(etherium_wb_clm$label_bug)
count(etherium_org$bug)
etherium_wb_clm_01<-etherium_wb_clm      #store the actial data in temp data with labelin with 0 1 as bug
etherium_wb_clm_01$label_bug<-ifelse(etherium_wb_clm_01$label_bug == "YES", 1,0)
etherium_wb_clm_mcc<-mccr(etherium_01b$bug,etherium_wb_clm_01$label_bug)
##Successfully above but result is so less
###################################################################################
#Now perform metric selection then instance selection then machine learning
#Start clmp
########################   Metric selection  started  ###############################################################################################
#Now perform Metric selection based on metric violation score MVS (metric violation score), vt (violation table)
#MVS of i metric=Ci/Fi,  Ci is no. of violation in ith metiric and Fi is the no. of metric values in the i-th metric(no. of instances)
#After labeling the instances now perform metric selection
#etherium_wb_clmp<-cbind(etherium_wb_clmp,
etherium_wb_clmp<-etherium_org[,-ncol(etherium_org)]
etherium_wb_clm_nRowsDf<-nrow(etherium_wb_clmp)
etherium_wb_clm_nColsDf<-ncol(etherium_wb_clmp)
etherium_wb_clm_vt<-matrix(nrow=etherium_wb_clm_nRowsDf,ncol=etherium_wb_clm_nColsDf)
for(i in 1:etherium_wb_clm_nColsDf)
{
  etherium_wb_clm_cy<-0  #no. of yes in coloumn cy
  etherium_wb_clm_cn<-0
  for(j in 1:etherium_wb_clm_nRowsDf)
  {
    if( etherium_wb_clm$label_bug[j]=='YES')
    {
      r <- if(etherium_wb_clm[j,i] <= etherium_wb_clm_thr[i]) 1 else 0
      etherium_wb_clm_cy<-etherium_wb_clm_cy+r
      etherium_wb_clm_vt[j,i]<-etherium_wb_clm[j,i]-etherium_wb_clm_thr[i]     #its for clami+ to write difference
    }
    else if( etherium_wb_clm$label_bug[j]=='NO')
    {
      s<- if(etherium_wb_clm[j,i] > etherium_wb_clm_thr[i]) 1 else 0
      
      etherium_wb_clm_cn<-etherium_wb_clm_cn+s
      etherium_wb_clm_vt[j,i]<-etherium_wb_clm[j,i]-etherium_wb_clm_thr[i]  #its difference table used in clami+
      
    }
  }
  etherium_wb_clm[etherium_wb_clm_nRowsDf+1,i]<-(etherium_wb_clm_cy+etherium_wb_clm_cn)  #the last row is added as the no. of violation row for each metrics. correct
}

#Now arrange the metrics in ascending order based on the value of last row and choose first ceilling(log2(n)) metrics for next phase.
#Firstly remove the last coloumn as label_bug for ordering metrics because its factor coloumn
etherium_wb_clm2<-etherium_wb_clm[,1:ncol(etherium_wb_clm)-1]
etherium_wb_clm2=etherium_wb_clm2[,order(etherium_wb_clm2[nrow(etherium_wb_clm2),])]

# select the metrics with minimum mvs score
etherium_wb_clm_mvs_min<-min(etherium_wb_clm[etherium_wb_clm_nRowsDf+1,1:etherium_wb_clm_nColsDf])  #min of 8th row
#Select the metrics which has 8th row contain the etherium_wb_clm_mvs_min value and store it in etherium_wb_clm_mvs means dataset with select metrics with min mvs
etherium_wb_clm_nmvs<- sum(etherium_wb_clm_mvs_min == etherium_wb_clm[etherium_wb_clm_nRowsDf+1,1:etherium_wb_clm_nColsDf])  #number of minimum vs
cat("number of selected metrics=",etherium_wb_clm_nmvs)
etherium_wb_clm_mvs<-matrix()
etherium_wb_clm_wmvs<-etherium_wb_clm[1:etherium_wb_clm_nRowsDf,]   #Removing the last row means mvs row from etherium_wb_clm, without mvs row
etherium_wb_clm_thr_st<-c()  #st means selected threshold of selected metrics

etherium_wb_clm_select_nm<-ceiling(log2(etherium_wb_clm_nColsDf))  #Select first ceiling(log2(n)) metrics based on mvs
etherium_wb_clm_ams<-matrix()   #etherium_wb_clm after metric selection
for(i in 1:etherium_wb_clm_select_nm)
{
  etherium_wb_clm_ams<-cbind(etherium_wb_clm_ams,etherium_wb_clm2[i])   #etherium_wb_clm2 used
}

etherium_wb_clm_ams<-cbind(etherium_wb_clm_ams[2:ncol(etherium_wb_clm_ams)],label_bug=etherium_wb_clm$label_bug)

str(etherium_wb_clm_ams)

##=======================================================================================
##metric selection completed now work on instance selection 
etherium_wb_clm_ams_nColsDf<-ncol(etherium_wb_clm_ams)-1
etherium_wb_clm_ams_nRowsDf<-nrow(etherium_wb_clm_ams)-1

for(i in 1:etherium_wb_clm_ams_nRowsDf)
{
  etherium_DS_iy<-0  #no. of yes in coloumn cy
  etherium_DS_in<-0
  for(j in 1:etherium_wb_clm_ams_nColsDf)
  {
    if( etherium_wb_clm_ams$label_bug[i]=='YES')
    {
      r <- if(etherium_wb_clm_ams[i,j] <= etherium_wb_clm_ams[etherium_wb_clm_ams_nRowsDf+1,j]) 1 else 0
      etherium_DS_iy<-etherium_DS_iy+r
      #etherium_DS_vt[j,i]<-etherium_wb_clm_ams[j,i]-etherium_DS_thr_st[i]     #its for clami+ to write difference
    }
    else if( etherium_wb_clm_ams$label_bug[i]=='NO')
    {
      s<- if(etherium_wb_clm_ams[i,j] > etherium_wb_clm_ams[etherium_wb_clm_ams_nRowsDf+1,j]) 1 else 0
      
      etherium_DS_in<-etherium_DS_in+s
      #etherium_DS_ivt[j,i]<-etherium_wb_clm_ams[j,i]-etherium_DS_thr[i]
    }
  }
  etherium_wb_clm_ams[i,etherium_wb_clm_ams_nColsDf+2]<-(etherium_DS_iy+etherium_DS_in)  #the 4th coloumn is the no. of violation row for each instances correct done 13.1 6/1/2021
}
#View(etherium_wb_clm_ams)
count(etherium_wb_clm_ams[,ncol(etherium_wb_clm_ams)])  #instance violation score completed   data[,ncol(data)]
## UP to here it executed completely   its final ds after metric selection
################################################################################
#Now arrange the instances in ascending order based on the value of last coloumn and choose first n-n/12 metrics for next phase.
#Firstly remove the last coloumn as label_bug for ordering rows because its factor coloumn
etherium_wb_clm_ams_row_ord<-etherium_wb_clm_ams[order(etherium_wb_clm_ams[,ncol(etherium_wb_clm_ams)]),]
etherium_wb_clm_ams_row_ord_nrow<-nrow(etherium_wb_clm_ams)
etherium_wb_clm_ams_row_ord<-etherium_wb_clm_ams_row_ord[-etherium_wb_clm_ams_row_ord_nrow,]
str(etherium_wb_clm_ams_row_ord)

etherium_wb_clm_ams_cnv<-count(etherium_wb_clm_ams_row_ord$V7)   #count the number of violations existed in the last coloumn V7
etherium_wb_clm_ams_cnv
etherium_wb_clm_ams_row_ord<-head(etherium_wb_clm_ams_row_ord,etherium_const_sr_clmp*etherium_wb_clm_ams_row_ord_nrow)
etherium_wb_clm_ams_is_nc<-count(etherium_wb_clm_ams_row_ord$label_bug)  #count the number of yes and no
etherium_wb_clm_ams_is_nc
etherium_wb_clm_ams_is_nc_r<-etherium_wb_clm_ams_is_nc[1,2]/etherium_wb_clm_ams_is_nc[2,2]  #ratio
etherium_wb_clm_ams_is_nc_r

etherium_wb_clm_ams_row_ord <- etherium_wb_clm_ams_row_ord[sample(nrow(etherium_wb_clm_ams_row_ord)),]   #randomize instances
etherium_wb_clm_ams_ais<-etherium_wb_clm_ams_row_ord[,-ncol(etherium_wb_clm_ams_row_ord)]                  #Remove last coloumn as ivs

#View(etherium_wb_clm_ams_ais)  #This is the final input data
######################################################################

#Now apply the SMOTE: Synthetic Minority Oversampling Technique To Handle Class Imbalancy In Binary Classification
count(etherium_wb_clm_ams_ais$label_bug)
etherium_wb_clm_ams_ais_rose<-ROSE(label_bug~.,data=etherium_wb_clm_ams_ais,seed = 786)$data   #it is better than smote
count(etherium_wb_clm_ams_ais_rose$label_bug)
set.seed(200)
etherium_wb_clm_ams_ais_rose<-etherium_wb_clm_ams_ais_rose[sample(nrow(etherium_wb_clm_ams_ais_rose)),]

#=======================================================================================================================
# #data<-read.csv(file.choose())  *********change 1:#m,#m+2
etherium_wb_clm_data<-etherium_wb_clm_ams_ais_rose     #label name   "label_bug"
etherium_clmp_org<-etherium_org                      #label name "bug"
etherium_wb_clm_predictors<-colnames(etherium_wb_clm_data[,1:ncol(etherium_wb_clm_data)-1])
#write.csv(etherium_wb_clm_predictors,"etherium_wb_clm_predictors.csv")
etherium_clmp_org_predictors<-colnames(etherium_clmp_org[,1:ncol(etherium_clmp_org)-1])
etherium_wb_clm_outcomeName<-'label_bug'
etherium_clmp_org_outcomeName<-'bug'

fitControl <- trainControl(      #train control for TCLP
  method = "cv",
  number = 10,
  savePredictions = 'final',  summaryFunction = multiClassSummary,
  classProbs = T )

set.seed(3)
etherium_wb_clm_model_rf<-train(etherium_wb_clm_data[,etherium_wb_clm_predictors],etherium_wb_clm_data[,etherium_wb_clm_outcomeName],method='naive_bayes',trControl=fitControl,tuneLength=5)
etherium_wb_clm_model_rf[["resample"]]<-na.mean(etherium_wb_clm_model_rf[["resample"]])
etherium_wb_clm_model_rf
etherium_clmp_org$predrf_bug<-predict(object = etherium_wb_clm_model_rf,etherium_clmp_org[,etherium_clmp_org_predictors],type='raw')   #predict on original data using pre trained model
etherium_wb_clm_cm_rf<-caret::confusionMatrix(factor(etherium_clmp_org$bug),factor(etherium_clmp_org$predrf_bug), mode= "prec_recall")
etherium_wb_clm_cm_rf
#View(etherium_clmp_org$predrf_bug)

etherium_wb_clmp_01<-etherium_clmp_org       #store the actual data in temp data with labelin with 0 1 as bug
etherium_wb_clmp_01$predrf_bug<-ifelse(etherium_wb_clmp_01$predrf_bug == "YES", 1,0)
etherium_wb_clmp_mcc<-mccr(etherium_01b$bug,etherium_wb_clmp_01$predrf_bug)
etherium_wb_clmp_mcc


#=============================================================================================================================
#Unsupervised Learning
#_wp means without preprocessing
etherium_ipkm_wp<- etherium_wb
set.seed(13)
etherium_ipkm_wp_kmean_obj<-kmeans(etherium_ipkm_wp,2,nstart = 50)
etherium_ipkm_wp_dist<-dist(etherium_ipkm_wp,method = "euclidean")^2

etherium_ipkm_wp_label<-cbind(etherium_ipkm_wp,prd_bug=etherium_ipkm_wp_kmean_obj$cluster)

etherium_ipkm_wp_label_count<-count(etherium_ipkm_wp_label$prd_bug)  
#Now change levels according cluster majority means if first large and second small then first else II

{
  if (etherium_ipkm_wp_label_count[1,2]>etherium_ipkm_wp_label_count[2,2])
  {
    etherium_ipkm_wp_label$prd_bug<-ifelse(etherium_ipkm_wp_label$prd_bug>1,'YES','NO')
    etherium_ipkm_wp_label$prd_bug<- factor(etherium_ipkm_wp_label$prd_bug, levels = c('NO','YES'))
    count(etherium_ipkm_wp_label$prd_bug)
  }
  else
  {
    etherium_ipkm_wp_label$prd_bug<-ifelse(etherium_ipkm_wp_label$prd_bug>1,'NO','YES')
    etherium_ipkm_wp_label$prd_bug<- factor(etherium_ipkm_wp_label$prd_bug, levels = c('YES','NO'))
    count(etherium_ipkm_wp_label$prd_bug)
  }
}

etherium_ipkm_wp_cm<-confusionMatrix(factor(etherium_org$bug),factor(etherium_ipkm_wp_label$prd_bug), mode="prec_recall")
etherium_ipkm_wp_cm
#View(etherium_ipkm_wp_label$prd_bug)

etherium_ipkm_wp_label_01<-etherium_ipkm_wp_label       
etherium_ipkm_wp_label_01$prd_bug<-ifelse(etherium_ipkm_wp_label_01$prd_bug == "YES", 1,0)
etherium_ipkm_wp_label_mcc<-mccr(etherium_01b$bug,etherium_ipkm_wp_label_01$prd_bug)
etherium_ipkm_wp_label_mcc
#===================================================================================================

#A Neural gas clusgering in cclust package
etherium_ipngc_wp<- etherium_wb
set.seed(13)

etherium_ipngc_wp_ngc_obj<-NeuralGasClustering(as.matrix(etherium_ipngc_wp), 2,PlotIt=FALSE)
etherium_ipngc_wp_label<-cbind(etherium_ipngc_wp,prd_bug=etherium_ipngc_wp_ngc_obj[["Object"]][["cluster"]])

etherium_ipngc_wp_label_count<-count(etherium_ipngc_wp_label$prd_bug)  
#Now change levels according cluster majority means if first large and second small then first else II

{
  if (etherium_ipngc_wp_label_count[1,2]>etherium_ipngc_wp_label_count[2,2])
  {
    etherium_ipngc_wp_label$prd_bug<-ifelse(etherium_ipngc_wp_label$prd_bug>1,'YES','NO')
    etherium_ipngc_wp_label$prd_bug<- factor(etherium_ipngc_wp_label$prd_bug, levels = c('NO','YES'))
    count(etherium_ipngc_wp_label$prd_bug)
  }
  else
  {
    etherium_ipngc_wp_label$prd_bug<-ifelse(etherium_ipngc_wp_label$prd_bug>1,'NO','YES')
    etherium_ipngc_wp_label$prd_bug<- factor(etherium_ipngc_wp_label$prd_bug, levels = c('YES','NO'))
    count(etherium_ipngc_wp_label$prd_bug)
  }
}

etherium_ipngc_wp_cm<-confusionMatrix(factor(etherium_org$bug),factor(etherium_ipngc_wp_label$prd_bug), mode="prec_recall")
etherium_ipngc_wp_cm

#View(etherium_ipngc_wp_label$prd_bug)

etherium_ipngc_wp_label_01<-etherium_ipngc_wp_label       
etherium_ipngc_wp_label_01$prd_bug<-ifelse(etherium_ipngc_wp_label_01$prd_bug == "YES", 1,0)
etherium_ipngc_wp_label_mcc<-mccr(etherium_01b$bug,etherium_ipngc_wp_label_01$prd_bug)
etherium_ipngc_wp_label_mcc
#============================================================================================================
# Model Based Clustering

etherium_ipmbc_wp<- etherium_wb
set.seed(13)
etherium_ipmbc_wp_mclust_obj<-Mclust(etherium_ipmbc_wp,2,nstart = 50)
etherium_ipmbc_wp_dist<-dist(etherium_ipmbc_wp,method = "euclidean")^2

etherium_ipmbc_wp_label<-cbind(etherium_ipmbc_wp,prd_bug=etherium_ipmbc_wp_mclust_obj[["classification"]])

etherium_ipmbc_wp_label_count<-count(etherium_ipmbc_wp_label$prd_bug)  
#Now change levels according cluster majority means if first large and second small then first else II

{
  if (etherium_ipmbc_wp_label_count[1,2]>etherium_ipmbc_wp_label_count[2,2])
  {
    etherium_ipmbc_wp_label$prd_bug<-ifelse(etherium_ipmbc_wp_label$prd_bug>1,'YES','NO')
    etherium_ipmbc_wp_label$prd_bug<- factor(etherium_ipmbc_wp_label$prd_bug, levels = c('NO','YES'))
    count(etherium_ipmbc_wp_label$prd_bug)
  }
  else
  {
    etherium_ipmbc_wp_label$prd_bug<-ifelse(etherium_ipmbc_wp_label$prd_bug>1,'NO','YES')
    etherium_ipmbc_wp_label$prd_bug<- factor(etherium_ipmbc_wp_label$prd_bug, levels = c('YES','NO'))
    count(etherium_ipmbc_wp_label$prd_bug)
  }
}

etherium_ipmbc_wp_cm<-confusionMatrix(factor(etherium_org$bug),factor(etherium_ipmbc_wp_label$prd_bug), mode="prec_recall")
etherium_ipmbc_wp_cm
#View(etherium_ipmbc_wp_label$prd_bug)

etherium_ipmbc_wp_label_01<-etherium_ipmbc_wp_label       
etherium_ipmbc_wp_label_01$prd_bug<-ifelse(etherium_ipmbc_wp_label_01$prd_bug == "YES", 1,0)
etherium_ipmbc_wp_label_mcc<-mccr(etherium_01b$bug,etherium_ipmbc_wp_label_01$prd_bug)
etherium_ipmbc_wp_label_mcc
#==========================================================================================================
# Hierarchical clustering Clustering
# Ward Hierarchical Clustering or complete
etherium_iphc_wp<- etherium_wb
set.seed(13)
etherium_iphc_wp_dist<-dist(etherium_iphc_wp,method = "euclidean")

etherium_iphc_wp_hclust_obj<-hclust(etherium_iphc_wp_dist,method="complete")   #method="complete"  ,ward.D2
#plot(etherium_iphc_wp_hclust_obj) # display dendogram
etherium_iphc_groups <- cutree(etherium_iphc_wp_hclust_obj, k=2) # cut tree into 2 clusters
# draw dendogram with red borders around the 5 clusters
#rect.hclust(etherium_iphc_wp_hclust_obj, k=2, border="red")

etherium_iphc_wp_label<-cbind(etherium_iphc_wp,prd_bug=etherium_iphc_groups)

etherium_iphc_wp_label_count<-count(etherium_iphc_wp_label$prd_bug)  
#Now change levels according cluster majority means if first large and second small then first else II

{
  if (etherium_iphc_wp_label_count[1,2]>etherium_iphc_wp_label_count[2,2])
  {
    etherium_iphc_wp_label$prd_bug<-ifelse(etherium_iphc_wp_label$prd_bug>1,'YES','NO')
    etherium_iphc_wp_label$prd_bug<- factor(etherium_iphc_wp_label$prd_bug, levels = c('NO','YES'))
    count(etherium_iphc_wp_label$prd_bug)
  }
  else
  {
    etherium_iphc_wp_label$prd_bug<-ifelse(etherium_iphc_wp_label$prd_bug>1,'NO','YES')
    etherium_iphc_wp_label$prd_bug<- factor(etherium_iphc_wp_label$prd_bug, levels = c('YES','NO'))
    count(etherium_iphc_wp_label$prd_bug)
  }
}

etherium_iphc_wp_cm<-confusionMatrix(factor(etherium_org$bug),factor(etherium_iphc_wp_label$prd_bug), mode="prec_recall")
etherium_iphc_wp_cm
#View(etherium_iphc_wp_label$prd_bug)

etherium_iphc_wp_label_01<-etherium_iphc_wp_label       
etherium_iphc_wp_label_01$prd_bug<-ifelse(etherium_iphc_wp_label_01$prd_bug == "YES", 1,0)
etherium_iphc_wp_label_mcc<-mccr(etherium_01b$bug,etherium_iphc_wp_label_01$prd_bug)
etherium_iphc_wp_label_mcc
#===============================================================================================

##############################################################################################
#Supervised Learning With 10 fold cross validation 
#############################################################################################

etherium_slcv_org<-etherium_org                      #label name "bug"
etherium_slcv_org_predictors<-colnames(etherium_slcv_org[,1:ncol(etherium_slcv_org)-1])
etherium_slcv_org_outcomeName<-'bug'

fitControl_2cv <- trainControl(
  method = "cv", 
  number = 10,  
  savePredictions = 'final',  summaryFunction = multiClassSummary,
  classProbs = T )

set.seed(44)
etherium_slcv_model_nb<-train(etherium_slcv_org[,etherium_slcv_org_predictors],etherium_slcv_org[,etherium_slcv_org_outcomeName],method='naive_bayes',trControl=fitControl_2cv,tuneLength=5)
etherium_slcv_model_nb[["resample"]]<-na.mean(etherium_slcv_model_nb[["resample"]])
etherium_slcv_model_nb
etherium_slcv_org$etherium_OOF_predcv_nb<-etherium_slcv_model_nb$pred$Y[order(etherium_slcv_model_nb$pred$rowIndex)]
etherium_slcv_org$etherium_OOF_predcv_nb<-na.mean(etherium_slcv_org$etherium_OOF_predcv_nb)


etherium_slcv_org$etherium_OOF_predcv_nbl<-as.factor(ifelse(etherium_slcv_org$etherium_OOF_predcv_nb>0.5,'YES','NO'))
etherium_slcv_org$etherium_OOF_predcv_nbl_01<-etherium_slcv_org$etherium_OOF_predcv_nbl
etherium_slcv_org$etherium_OOF_predcv_nbl_01<-ifelse(etherium_slcv_org$etherium_OOF_predcv_nbl_01 == "YES", 1,0)
etherium_slcv_mcc_nb<-mccr(etherium_01b$bug,etherium_slcv_org$etherium_OOF_predcv_nbl_01)

set.seed(45)
etherium_slcv_model_svm<-train(etherium_slcv_org[,etherium_slcv_org_predictors],etherium_slcv_org[,etherium_slcv_org_outcomeName],method='svmRadial',trControl=fitControl_2cv,tuneLength=5)
etherium_slcv_model_svm[["resample"]]<-na.mean(etherium_slcv_model_svm[["resample"]])
etherium_slcv_model_svm
etherium_slcv_org$etherium_OOF_predcv_svm<-etherium_slcv_model_svm$pred$Y[order(etherium_slcv_model_svm$pred$rowIndex)]
etherium_slcv_org$etherium_OOF_predcv_svm<-na.mean(etherium_slcv_org$etherium_OOF_predcv_svm)


etherium_slcv_org$etherium_OOF_predcv_svml<-as.factor(ifelse(etherium_slcv_org$etherium_OOF_predcv_svm>0.5,'YES','NO'))
etherium_slcv_org$etherium_OOF_predcv_svetherium_01<-etherium_slcv_org$etherium_OOF_predcv_svml
etherium_slcv_org$etherium_OOF_predcv_svetherium_01<-ifelse(etherium_slcv_org$etherium_OOF_predcv_svetherium_01 == "YES", 1,0)
etherium_slcv_mcc_svm<-mccr(etherium_01b$bug,etherium_slcv_org$etherium_OOF_predcv_svetherium_01)

set.seed(46)
etherium_slcv_model_knn<-train(etherium_slcv_org[,etherium_slcv_org_predictors],etherium_slcv_org[,etherium_slcv_org_outcomeName],method='knn',trControl=fitControl_2cv,tuneLength=5)
etherium_slcv_model_knn[["resample"]]<-na.mean(etherium_slcv_model_knn[["resample"]])
etherium_slcv_model_knn
etherium_slcv_org$etherium_OOF_predcv_knn<-etherium_slcv_model_knn$pred$Y[order(etherium_slcv_model_knn$pred$rowIndex)]
etherium_slcv_org$etherium_OOF_predcv_knn<-na.mean(etherium_slcv_org$etherium_OOF_predcv_knn)


etherium_slcv_org$etherium_OOF_predcv_knnl<-as.factor(ifelse(etherium_slcv_org$etherium_OOF_predcv_knn>0.5,'YES','NO'))
etherium_slcv_org$etherium_OOF_predcv_knnl_01<-etherium_slcv_org$etherium_OOF_predcv_knnl
etherium_slcv_org$etherium_OOF_predcv_knnl_01<-ifelse(etherium_slcv_org$etherium_OOF_predcv_knnl_01 == "YES", 1,0)
etherium_slcv_mcc_knn<-mccr(etherium_01b$bug,etherium_slcv_org$etherium_OOF_predcv_knnl_01)

set.seed(47)
etherium_slcv_model_rf<-train(etherium_slcv_org[,etherium_slcv_org_predictors],etherium_slcv_org[,etherium_slcv_org_outcomeName],method='rf',trControl=fitControl_2cv,tuneLength=5)
etherium_slcv_model_rf[["resample"]]<-na.mean(etherium_slcv_model_rf[["resample"]])
etherium_slcv_model_rf
etherium_slcv_org$etherium_OOF_predcv_rf<-etherium_slcv_model_rf$pred$Y[order(etherium_slcv_model_rf$pred$rowIndex)]
etherium_slcv_org$etherium_OOF_predcv_rf<-na.mean(etherium_slcv_org$etherium_OOF_predcv_rf)


etherium_slcv_org$etherium_OOF_predcv_rfl<-as.factor(ifelse(etherium_slcv_org$etherium_OOF_predcv_rf>0.5,'YES','NO'))
etherium_slcv_org$etherium_OOF_predcv_rfl_01<-etherium_slcv_org$etherium_OOF_predcv_rfl
etherium_slcv_org$etherium_OOF_predcv_rfl_01<-ifelse(etherium_slcv_org$etherium_OOF_predcv_rfl_01 == "YES", 1,0)
etherium_slcv_mcc_rf<-mccr(etherium_01b$bug,etherium_slcv_org$etherium_OOF_predcv_rfl_01)

set.seed(48)
etherium_slcv_model_C50<-train(etherium_slcv_org[,etherium_slcv_org_predictors],etherium_slcv_org[,etherium_slcv_org_outcomeName],method='C5.0',trControl=fitControl_2cv,tuneLength=5)
etherium_slcv_model_C50[["resample"]]<-na.mean(etherium_slcv_model_C50[["resample"]])   #naive_bayes, nb
etherium_slcv_model_C50

etherium_slcv_org$etherium_OOF_predcv_C50<-etherium_slcv_model_C50$pred$Y[order(etherium_slcv_model_C50$pred$rowIndex)]
etherium_slcv_org$etherium_OOF_predcv_C50<-na.mean(etherium_slcv_org$etherium_OOF_predcv_C50)


etherium_slcv_org$etherium_OOF_predcv_C50l<-as.factor(ifelse(etherium_slcv_org$etherium_OOF_predcv_C50>0.5,'YES','NO'))
etherium_slcv_org$etherium_OOF_predcv_C50l_01<-etherium_slcv_org$etherium_OOF_predcv_C50l
etherium_slcv_org$etherium_OOF_predcv_C50l_01<-ifelse(etherium_slcv_org$etherium_OOF_predcv_C50l_01 == "YES", 1,0)
etherium_slcv_mcc_C50<-mccr(etherium_01b$bug,etherium_slcv_org$etherium_OOF_predcv_C50l_01)



# #======================================================================================================================

#combine accuracy nd f-measure and mcc auc kappa
etherium_accuracy<-matrix(c(etherium_wb_cla_cm_rf[["overall"]][["Accuracy"]]*100,etherium_wb_tcl_cm[["overall"]][["Accuracy"]]*100,etherium_wb_tcl_cm_rf[["overall"]][["Accuracy"]]*100,etherium_wb_acl_cm[["overall"]][["Accuracy"]]*100,etherium_wb_acl_cm_rf[["overall"]][["Accuracy"]]*100,etherium_wb_clm_cm[["overall"]][["Accuracy"]]*100,etherium_wb_clm_cm_rf[["overall"]][["Accuracy"]]*100,
                            
                            etherium_ipkm_wp_cm[["overall"]][["Accuracy"]]*100,etherium_ipngc_wp_cm[["overall"]][["Accuracy"]]*100,etherium_ipmbc_wp_cm[["overall"]][["Accuracy"]]*100,etherium_iphc_wp_cm[["overall"]][["Accuracy"]]*100,
                            mean(etherium_slcv_model_nb[["results"]][["Accuracy"]])*100,mean(etherium_slcv_model_svm[["results"]][["Accuracy"]])*100,mean(etherium_slcv_model_knn[["results"]][["Accuracy"]])*100,
                            mean(etherium_slcv_model_rf[["results"]][["Accuracy"]])*100,mean(etherium_slcv_model_C50[["results"]][["Accuracy"]])*100
),ncol=16,byrow=TRUE)
colnames(etherium_accuracy)<-c('CLAMI','TCL','TCLP','ACL','ACLP','CLM','CLMP','KMS','NGC','MBC','HC','NB','SVM','KNN','RF','C50')
rownames(etherium_accuracy) <- c("etherium_Accuracy")
etherium_accuracy <- as.table(etherium_accuracy)
etherium_accuracy

etherium_fm<-matrix(c(etherium_wb_cla_cm_rf[["byClass"]][["F1"]],etherium_wb_tcl_cm[["byClass"]][["F1"]],etherium_wb_tcl_cm_rf[["byClass"]][["F1"]],etherium_wb_acl_cm[["byClass"]][["F1"]],etherium_wb_acl_cm_rf[["byClass"]][["F1"]],etherium_wb_clm_cm[["byClass"]][["F1"]],etherium_wb_clm_cm_rf[["byClass"]][["F1"]],
                      
                      etherium_ipkm_wp_cm[["byClass"]][["F1"]],etherium_ipngc_wp_cm[["byClass"]][["F1"]],etherium_ipmbc_wp_cm[["byClass"]][["F1"]],etherium_iphc_wp_cm[["byClass"]][["F1"]],
                      mean(etherium_slcv_model_nb[["results"]][["F1"]]),mean(etherium_slcv_model_svm[["results"]][["F1"]]),mean(etherium_slcv_model_knn[["results"]][["F1"]]),
                      mean(etherium_slcv_model_rf[["results"]][["F1"]]),mean(etherium_slcv_model_C50[["results"]][["F1"]])
),ncol=16,byrow=TRUE)
colnames(etherium_fm)<-c('CLAMI','TCL','TCLP','ACL','ACLP','CLM','CLMP','KMS','NGC','MBC','HC','NB','SVM','KNN','RF','C50')
rownames(etherium_fm) <- c("etherium_fm")
etherium_fm <- as.table(etherium_fm)
etherium_fm
#etherium_cm_mv3[["byClass"]][["F1"]],max(etherium_slcv_model_C50[["results"]][["F1"]])
#View(etherium_slcv_model_rf)

etherium_mcc<-matrix(c(etherium_wb_clami_mcc,etherium_wb_tcl_mcc,etherium_wb_tclp_mcc,etherium_wb_acl_mcc,etherium_wb_aclp_mcc,etherium_wb_clm_mcc,etherium_wb_clmp_mcc,
                       etherium_ipkm_wp_label_mcc,etherium_ipngc_wp_label_mcc,etherium_ipmbc_wp_label_mcc,etherium_iphc_wp_label_mcc
                       ,etherium_slcv_mcc_nb,etherium_slcv_mcc_svm,etherium_slcv_mcc_knn,etherium_slcv_mcc_rf,etherium_slcv_mcc_C50
),ncol=16,byrow=TRUE)
colnames(etherium_mcc)<-c('CLAMI','TCL','TCLP','ACL','ACLP','CLM','CLMP','KMS','NGC','MBC','HC','NB','SVM','KNN','RF','C50')
rownames(etherium_mcc) <- c("etherium_mcc")
etherium_mcc <- as.table(etherium_mcc)
etherium_mcc

etherium_cmb_table<-matrix(c(etherium_org_ncol,etherium_org_nrow,etherium_org_bugpr,etherium_accuracy,etherium_fm,etherium_mcc),ncol = 51,byrow=FALSE)
colnames(etherium_cmb_table)<-c('NCOL','NROW','bug%','CLAMI','TCL','TCLP','ACL','ACLP','CLM','CLMP','KMS','NGC','MBC','HC','NB','SVM','KNN','RF','C50','CLAMI','TCL','TCLP','ACL','ACLP','CLM','CLMP','KMS','NGC','MBC','HC','NB','SVM','KNN','RF','C50','CLAMI','TCL','TCLP','ACL','ACLP','CLM','CLMP','KMS','NGC','MBC','HC','NB','SVM','KNN','RF','C50')
rownames(etherium_cmb_table) <- c("etherium_DS")
etherium_cmb_table<-as.table(etherium_cmb_table)
write.csv(etherium_cmb_table,"etherium_cmb_table_TCLP.csv")
