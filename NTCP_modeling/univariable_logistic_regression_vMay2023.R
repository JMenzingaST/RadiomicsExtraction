univariable_logistic_regression <- function(DATA,ynam) {
 
  result <- tryCatch({
    
  feats_of_interest=colnames(DATA)[!colnames(DATA) %in% ynam]

  Model_predict=list(type = any)
  stats.s=matrix()
  OR_95=matrix(NA,nrow = length(feats_of_interest), ncol=4)
  AIC_BIC=matrix(NA,nrow = length(feats_of_interest), ncol=2)
  stats.s=matrix(NA,nrow = length(feats_of_interest), ncol=19)
  nums=matrix(NA,nrow = length(feats_of_interest), ncol=2)
  for (i in 1:length(feats_of_interest)) {
    frm <- as.formula(paste( ynam, " ~ ",feats_of_interest[i]))# LGRE_GLRLM_max
    Model_predict[[i]] = lrm(formula = frm, x = T, y = T,  data = DATA) # , x = T, y = T, data = pred.frame
    
    StanER=sqrt(diag(vcov( Model_predict[[i]]))) 
    OR_95[i,]=cbind(Model_predict[[i]]$coefficients[[feats_of_interest[i]]], exp(Model_predict[[i]]$coefficients[[feats_of_interest[i]]]),
                    exp(Model_predict[[i]]$coefficients[[feats_of_interest[i]]]-StanER[[feats_of_interest[i]]]*1.96),
                    exp(Model_predict[[i]]$coefficients[[feats_of_interest[i]]]+StanER[[feats_of_interest[i]]]*1.96))
    stats.s[i,]=c(Model_predict[[i]]$stats,anova(Model_predict[[i]])[1,3])
    AIC_BIC[i,]=c(AIC( Model_predict[[i]]),BIC( Model_predict[[i]]))
    nums[i,]=c(Model_predict[[i]] $freq[2],sum(Model_predict[[i]] $freq))

  }
  
  colnames(stats.s)=names(c(Model_predict[[i]]$stats,anova(Model_predict[[i]])[1,3]))
  colnames(OR_95)=c("beta","OR","OR95-","OR95+")
  colnames(AIC_BIC)=c("AIC","BIC")
  colnames(nums)=c("event","n")
  rownames(stats.s)<- feats_of_interest;colnames(stats.s)[19]="p_LRT"
  
  index_p=order(stats.s[,"p_LRT"])
  # index_p=order(AIC_BIC[,"AIC"])
  univar_result=cbind(stats.s[index_p,c(19,6,10)],  OR_95[index_p,],AIC_BIC[index_p,],nums[index_p,])
  
  
  # print(univar_result)#[univar_uitk2[,'uP']<0.05,])
  # cat("\n","Univariable significant: ", nrow(univar_uitk2[univar_uitk2[,'p_LRT']<0.05,])," of ", nrow(univar_uitk2))
  
  return(univar_result)
 
  
  
  
  
  }, error = function(e) {
    print(frm)
  })
  
  
}