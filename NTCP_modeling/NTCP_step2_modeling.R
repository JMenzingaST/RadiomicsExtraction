rm(list =ls())
#Sys.setenv(JAVA_HOME="");library(rJava)
library(rms);library(plyr);#library(xlsx)
library(pracma);library(caret);library(ResourceSelection)
library(pROC)

outputdir = r"(\Modeling\)"

normaliseer <- function(b) {
  db <- dim(b)
  g=b
  del=matrix();meansd=c()
  b_n<-matrix(nrow=db[1], ncol=db[2]) 
  for (f in 1:db[2]) {
    g[,f] <- (b[,f]-mean(b[,f]))/sd(b[,f])
  }
  
  return(g)
  
}

###############################
# Test on whole clin data set #
###############################

load(r"(\preprocessed_clinical_variables.R)")
clin_data = DATA

# fit a simple univariable NTCP model
dd <- datadist(clin_data); options(datadist='dd')
Model_ref0 <- lrm(y ~  LungMINGTV_doseGEM + Age_StartRT + smoking_3months, x = T, y = T, data = clin_data)               #MODIFY to a specific/random variable in your set
Model_ref0

# load Dataset
DATA <- read.csv(r"(\merged_clinvar_radiomics_all.csv)")

##### Doing the Train Test Split #####
missing_v=apply(DATA,2, function(x) sum(is.na(x)))
print(sort(missing_v[missing_v>5]))

# delete buggy variables

DATA=DATA[,!colnames(DATA) %in% c("Patient_ID","train","Index", 
                                  grep('Sphericity', colnames(DATA), value=T), 
                                  grep('ClusterShade', colnames(DATA), value=T),
                                  grep('ClusterProminence', colnames(DATA), value=T),
                                  grep('Complexity', colnames(DATA), value=T))]
DATA = na.omit(DATA)
DATA$smoking_3months = as.factor(DATA$smoking_3months)
DATA = na.omit(DATA)

print(paste0("event rate: ",round(100*sum(DATA$y)/nrow(DATA),0),"%;  ",sum(DATA$y),"pts out of total ",nrow(DATA) ))

set.seed(28)
trainperc=0.7; all_index=1:nrow(DATA)
train_sample=sort(sample(all_index,round(trainperc*nrow(DATA))))
test_sample=all_index[! all_index %in% train_sample]

DATA=DATA[train_sample,]
DATA_test=DATA[test_sample,]

print(paste0("event rate: ",round(100*sum(DATA$y)/nrow(DATA),0),"%;  ",sum(DATA$y),"pts out of total ",nrow(DATA) ))

# Only leave in Niezink variables, Y and radiomics:
avg_radiomics <- grep('avg', colnames(DATA), value=T)
bp60_radiomics <- grep('bp60', colnames(DATA), value=T)
bp0_radiomics <- grep('bp0', colnames(DATA), value=T)
delta_radiomics <- grep('delta', colnames(DATA), value=T)

niezink_cols = c("Patient_ID", "y" ,"LungMINGTV_doseGEM", 
                 "Age_StartRT", "smoking_3months")

random_cols = grep('Random', colnames(DATA), value = T)

###############################################################################
####### Change this to vary the radiomics you want to make models with! #######
DATA = DATA[, colnames(DATA) %in% c(niezink_cols, delta_radiomics, random_cols)]
DATA_test = DATA_test[, colnames(DATA_test) %in% c(niezink_cols, delta_radiomics, random_cols)]
###############################################################################

ynam="y"

feats_of_interest=colnames(DATA)[!colnames(DATA) %in% c(ynam)]
# t(t(colnames(DATA))) # This is just a 'print' statement to look at all the features, one feature per line
radiomics_cols = grep('original', colnames(DATA), value=T)
DATA[,colnames(DATA) %in% radiomics_cols] <- normaliseer(DATA[,colnames(DATA) %in% radiomics_cols])

DATA_num=DATA
for (io in feats_of_interest){
  DATA_num[,io]=as.numeric(DATA[,io]) -1
}

dd <- datadist(DATA_num);options(datadist='dd')

###########################
##### REFERENCE MODEL #####
###########################

##### Refit of Niezink Model #####
niezink_refit_model <- lrm(y ~ MLD + Age + Smoking, x=T, y=T, data=DATA_num)
niezink_refit_model

# Evaluation
NTCP_niezink_refit_train <- predict(niezink_refit_model, DATA, type="fitted")
roc_niezink_refit_train <- roc(DATA$y, NTCP_niezink_refit_train, ci = T, plot = T)
print(paste0("Train AUC: ", roc_niezink_refit_train$auc)); print(roc_niezink_refit_train$ci)

NTCP_niezink_refit_test <- predict(niezink_refit_model, DATA_test, type="fitted")
roc_niezink_refit_test <- roc(DATA_test$y, NTCP_niezink_refit_test,ci = T, plot = T)
print(paste0("Test AUC: ", roc_niezink_refit_test$auc)); print(roc_niezink_refit_test$ci)

######################
###### ANALYSES ######
######################

########## univariable analyses ##########
source(r"(\univariable_logistic_regression_vMay2023.R)")

oo <- univariable_logistic_regression(DATA_num, ynam)
print(round(oo,3)); cat("\n","Univariable significant: ", nrow(oo[oo[,'p_LRT']<0.05,])," of ", nrow(oo))

# Save the Univariate Logistic regression outcomes
write.csv(oo, paste0(outputdir,"univariable_regression_pneumonitis.csv"))

# Some plots #
plot(DATA$LungMINGTV_doseGEM,DATA$original_glrlm_GrayLevelNonUniformityNormalized,col=c("black","red")[1+DATA$y],type= 'p',pch = c(16,16))

#############################
##### Forward selection #####
#############################

## pre-selection (optional)
source("/preselection_correlation_v2.R")
Preslection <- preselection_correlation(DATA_num, ynam)#
DATA_sel <- DATA_num[,colnames(DATA_num) %in% c("y",rownames(Preslection))]# Alter base data frame

DATA_sel = na.omit(DATA_sel)
dd <- datadist(DATA_sel)
options(datadist='dd')
# # Nu door met deze nieuwe pre-selected dataset (DATA_sel)
source(r"(\forward_sel_logistic.R)")

run_i <- 0;
excel_export <- 0
output_fs <- F
boot <- 0
go <- forward_sel_logistic(DATA_sel,ynam,boot,excel_export,output_fs, 0.10)
round(go,2)
write.table(go, file="clipboard-16384", sep="\t", row.names=T, col.names=T)

############################################
###### Bootstrapped Forward selection ######
############################################
presel=0
iterations=1000
inTrain <- createResample (DATA$y, iterations, list = FALSE)
totaal_alles=1:size(inTrain,1);
for (run_i in 1:iterations){
  test_set=totaal_alles[!totaal_alles %in% inTrain[,run_i]]
  # print(length(test_set))
}

statistiek=matrix(NA,nrow = iterations, ncol = 9);coef_nam=matrix(NA,nrow = iterations, ncol = 20);coef_num=matrix(NA,nrow = iterations, ncol = 20);coef_pval=matrix(NA,nrow = iterations, ncol = 20);
cal_slopen=matrix(NA,nrow = iterations, ncol = 2);dis_slopen=matrix(NA,nrow = iterations, ncol = 2)
for (run_i in 1:iterations){
  if (run_i>1){rm(DATA_boot)}
  
  DATA_boot=DATA_num[inTrain[,run_i],]
  print(paste0('N.o. events: ',sum(DATA_boot$y)))
  
  ###pre-selection
  print("pre-selection")
  dd <- datadist(DATA_boot); options(datadist='dd')
  if (presel>0){
    Preslection=preselection_correlation(DATA_boot,ynam)
  }
  
  ###forward selection
  print('Forward Selection')
  go=try_default(forward_sel_logistic(DATA_boot,ynam,0,0,"",0.05), 0)
  int=grep("Intercept",rownames(go))
  if (length(int)>1){
    
    final_model_go= go[int[length(int)-1]:(int[length(int)]-1),]
    
    print(run_i)
    
    lang=length(final_model_go[,1])
    final_vars=names(final_model_go[,1])[2:lang]
    coef_nam[run_i,1:(lang-1)]=final_vars
    coef_num[run_i,1:lang]=final_model_go[c(2:(lang),1),1]
    coef_pval[run_i,1:lang]=final_model_go[c(2:(lang),1),3]
    
    ### final model performance in bootstrap sample
    print('final model performance in bs')
    frm <- as.formula(paste("y ~ ",paste(final_vars, collapse = ' + ')) )
    finalmodel = lrm(formula = frm, x = T, y = T,  data = DATA_boot)
    predict.fit2 <- predict(finalmodel,DATA_boot, type="fitted")
    q.fit<-  t.test(predict.fit2~DATA_boot$y)
    dis_slopen[run_i,1] <- q.fit$estimate[2]-q.fit$estimate[1]
    
    ### final model performance in entire set
    print('fin mod perf in ent set')
    NTCP <- predict(finalmodel,DATA_num, type="fitted")
    lp=  log(NTCP/(1-NTCP))
    
    if (sum(is.infinite(lp))==0){
      fit.enter<- lrm(y ~ lp, x = T, data = DATA_num)
      validatie_stat=fit.enter$stats[c(5:6,10)];names(validatie_stat)=paste(names(validatie_stat),'_val',sep='')
      cal_slopen[run_i,1:2]=fit.enter$coefficients
      q.fit<-  t.test(NTCP~DATA_num$y)
      dis_slopen[run_i,2] <-  q.fit$estimate[2]-q.fit$estimate[1]
    } else {validatie_stat=c(NA,NA,NA);names(validatie_stat)=paste(names(validatie_stat),'_val',sep='')}
    
    
    statistiek[run_i,]=c(final_model_go[1,c(4:6,11:13)],validatie_stat)
    colnames(statistiek)=c("P"  ,    "C"     , "R2" , "AIC","BIC",   "p_LRT",   "P_val",  "C_val" , "R2_val")
  }
}

### post processing stappen
colnames(cal_slopen)=c("intercept_lp","Slope_lp");colnames(dis_slopen)=c("DS_boot","DS_val")

C_diff=cbind(statistiek[,"P"]-statistiek[,"P_val"],statistiek[,"C"]-statistiek[,"C_val"],statistiek[,"R2"]-statistiek[,"R2_val"],dis_slopen[,2]-dis_slopen[,1])
colnames(C_diff)=c("P_diff","C_diff","R2_diff","DS_diff")
statistiek=cbind(statistiek,dis_slopen,C_diff,cal_slopen)

coef_nam_num_pval=c()
for (t in 1:20){coef_nam_num_pval=cbind(coef_nam_num_pval, cbind(coef_nam[,t],coef_num[,t],coef_pval[,t]))}

print(mean(statistiek[!is.na(statistiek[,"C_diff"]),"C_diff"]))

### calculate frequencies of selection
feats_of_interest=colnames(DATA)[!colnames(DATA) %in% ynam]
frequencies=matrix(0,nrow = length(feats_of_interest),ncol = 5);rownames(frequencies) = feats_of_interest
for (t in 1:length(feats_of_interest)){
  frequencies[t,1]=sum(coef_nam %in% feats_of_interest[t])
  frequencies[t,2]=sum(coef_nam[,1] %in% feats_of_interest[t])
  frequencies[t,3]=sum(coef_nam[,2] %in% feats_of_interest[t])
  frequencies[t,4]=sum(coef_nam[,3] %in% feats_of_interest[t])
  frequencies[t,5]=sum(coef_nam[,4] %in% feats_of_interest[t])
}
frequencies=frequencies[order(-frequencies[,1]),]
frequencies

filename="BootstrappedFwdSel"

write.csv(frequencies, paste(outputdir,paste(filename,"_FREQUENCIES__",iterations,".csv",sep=""),sep="/"))
write.csv(cbind(statistiek,coef_nam_num_pval),    paste(outputdir,paste(filename,"_stat_coef_nam__"   ,iterations,"csv",sep=""),sep="/"))


#########################################
##### Model Creation and Validation #####
#########################################
bp0_lung_minus_tumor
# Creation
final_model<- lrm(y ~ MLD + Age + Smoking + delta_originalshape_Maximum2DDiameterSlice + delta_originalngtdm_Strength,
                  x = T, y = T, data = DATA_num)
final_model

# Validation
NTCP_new_train <- predict(final_model, DATA, type="fitted")
roc_new_train <- roc(DATA$y, NTCP_new_train, ci = T, plot = T)
print(paste0("Train AUC: ", roc_new_train$auc)); print(roc_new_train$ci)

NTCP_new_test <- predict(final_model, DATA_test, type="fitted")
roc_new_test <- roc(DATA_test$y, NTCP_new_test, ci = T, plot = T)
print(paste0("Test AUC: ", roc_new_test$auc)); print(roc_new_test$ci)

# Correlation between radiomics features
cor(DATA_num$delta_originalshape_Maximum2DDiameterSlice, 
    DATA_num$delta_originalngtdm_Strength,
    method = 'spearman')

##################
##### REFITS ##### 
##################

##### Refit of Niezink Model #####
niezink_refit_model <- lrm(y ~ MLD + Age + Smoking, x=T, y=T, data=DATA_num)
niezink_refit_model

# Evaluation
NTCP_niezink_refit_train <- predict(niezink_refit_model, DATA, type="fitted")
roc_niezink_refit_train <- roc(DATA$y, NTCP_niezink_refit_train, ci = T, plot = T)
print(paste0("Train AUC: ", roc_niezink_refit_train$auc)); print(roc_niezink_refit_train$ci)

NTCP_niezink_refit_test <- predict(niezink_refit_model, DATA_test, type="fitted")
roc_niezink_refit_test <- roc(DATA_test$y, NTCP_niezink_refit_test,ci = T, plot = T)
print(paste0("Test AUC: ", roc_niezink_refit_test$auc)); print(roc_niezink_refit_test$ci)