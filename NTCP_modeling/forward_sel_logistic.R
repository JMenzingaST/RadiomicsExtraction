forward_sel_logistic <- function(DATA,ynam,boot,excel_export,outputdir,p_val_thresh) {

#force in canindate variables
sel=c('Age_StartRT', 'LungMINGTV_doseGEM', 'smoking_3months');
TEST_pvalue=0;#p_val_thresh=0.05

run_i0=0
feats_of_interest=colnames(DATA)[!colnames(DATA) %in% ynam]
Steps_together=c()
while (TEST_pvalue<p_val_thresh & length(sel)<length(feats_of_interest) ) {
  run_i0=run_i0+1
  
  ####-- delete variables already in the model --
  selfor=sel[1];if (length(sel)>1){ for (z in 2:length(sel)) {selfor=paste(selfor,sel[z],sep = " + ")}};selfor
  feats_of_interest=feats_of_interest[!feats_of_interest %in%  sel]
  
  
  ####-- make reference model
  if (length(sel)>0){
    frm <- as.formula(paste("y ~ ",selfor) )
    Model_ref = lrm(formula = frm, x = T, y = T,  data = DATA) }  
  
  
  Model_predict=list(type = any)
  stats.s   =matrix(NA,nrow = 18, ncol = length(feats_of_interest));
  stats.val =matrix(NA,nrow = 11, ncol = length(feats_of_interest));
  stats.LLRT=matrix(NA,nrow = 1,  ncol = length(feats_of_interest))
  AIC_BIC   =matrix(NA,nrow = 2,  ncol = length(feats_of_interest))
  beta      =matrix(NA,nrow = 2,  ncol = length(feats_of_interest))
  for (i in 1:length(feats_of_interest)) {
    
    frm <- as.formula(paste("y ~ ",selfor," + ",feats_of_interest[i])) 
    Model_predict[[i]] = try(
      (lrm(formula = frm, x = T, y = T,  data = DATA) ))
    
    if(inherits(Model_predict[[i]], "try-error")){
      stats.LLRT[i]=NA;
      }else if(length( Model_predict[[i]])<3  ) {
         h="k"
    } else {
      ###-- basic performance measures
      stats.s[,i]=Model_predict[[i]]$stats
      AIC_BIC[,i]=c(AIC(Model_predict[[i]]),BIC( Model_predict[[i]]))
      beta[,i] =c(Model_predict[[i]]$coefficients[[feats_of_interest[i]]], exp(Model_predict[[i]]$coefficients[[feats_of_interest[i]]]))
      ###-- bootstrap per model
      if (boot>0){   
        p=validate(Model_predict[[i]],method = "boot", B=boot); #print(p)
        stats.val[,i]=p[,5];
      } else { stats.val[,i]=NA} 
      
      ###-- likelihood ratio test
      if (length(sel)>0) {LLRT=lrtest(Model_ref,Model_predict[[i]]);stats.LLRT[i]=LLRT$stats[3]
      }else {LLRT=anova(Model_predict[[i]]);stats.LLRT[i]=LLRT["TOTAL","P"]}
      
    }
  }
  
  
  colnames(stats.s)  = feats_of_interest;
  rownames(stats.s)  = names(Model_predict[[i]]$stats);rownames(stats.LLRT)="p_LRT"
  rownames(stats.val)= paste0("val_",c( "Dxy","R2","Intercept","Slope","Emax","D","U","Q","B","g","gp"))
  rownames(AIC_BIC)  = c("AIC","BIC")
  rownames(beta)     = c("beta","OR")
  
  stats.val[1,]=(stats.val[1,]+1)/2
  # rownames(stats.val)<-paste("val_",rownames(p),sep = '');
  
  vallalles=t(rbind(stats.s[c(5:6,10),],stats.val[1:2,],beta,AIC_BIC,stats.LLRT));
  
  vgl=vallalles[order(vallalles[,"p_LRT"]),]
  # is0=which(vgl[,"p_LRT"]==0);vgl_is0=vgl[is0,];
  # vgl[is0,]=vgl_is0[order(-vgl_is0[,"R2"]),]
  
  round(vgl,3)
  
  TEST_pvalue=vgl[1,"p_LRT"]; if (is.na(as.numeric(TEST_pvalue))){TEST_pvalue=1}
  
  selected_variable=rownames(vgl)[1]
  sel=c(sel,selected_variable)
  
  print(paste0(round(TEST_pvalue,4),"-->  ",selfor," + ",selected_variable ))
  # 

  # print(run_i)
  
  if (excel_export==1){
    write.xlsx(vgl, paste(outputdir,paste("Preselection_step_",run_i0,".xlsx",sep=""),sep="\\"), row.names=TRUE)
  }
  
  
  frm <- as.formula(paste("y ~ " ,selfor," + ",selected_variable) )
  step_model<-lrm(formula = frm, x = T, y = T,  data = DATA) 
  ss=summary(step_model)
  step_model_coef=cbind(step_model$coefficients, exp(step_model$coefficients),c(NA,anova( step_model)[1:(length(step_model$coefficients)-1),3]))
  
  step_model_params=matrix(NA, size(step_model_coef)[1],3+length(vgl[1,]));
  step_model_params[,1:3]=step_model_coef
  step_model_params[1,4:(3+length(vgl[1,]))]=vgl[1,]
  if (length(sel)>1){
    rownames(step_model_params)=rownames(step_model_coef)
    colnames(step_model_params)=c(c("beta","OR","p_val"),  names(vgl[1,]))
  }else{
    rownames(step_model_params)=rownames(step_model_coef)
    colnames(step_model_params)=c(c("beta","OR","p_val"),  names(vgl[1,]))
  }
  
  Steps_together=rbind(Steps_together,step_model_params)
}
round(Steps_together,2)

return(Steps_together)
}