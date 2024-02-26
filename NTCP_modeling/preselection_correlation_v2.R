#Building block: preselection_correlation

preselection_correlation <- function(DATA,ynam) {
  
feats_of_interest=colnames(DATA)#[,c(2:133,135:274,276:(dim(DATA)[2]))]);
 feats_of_interest <- feats_of_interest[!feats_of_interest %in% c(ynam)];t(t(feats_of_interest))
dd <- datadist(DATA)
options(datadist='dd')
Model_predict=list(type = any)
stats.s=matrix()
OR_95=matrix(NA,nrow = length(feats_of_interest), ncol=4)

for (i in 1:length(feats_of_interest)) {
  frm <- as.formula(paste( ynam, " ~ ",feats_of_interest[i]))# LGRE_GLRLM_max
  Model_predict[[i]] = lrm(formula = frm, x = T, y = T,  data = DATA) # , x = T, y = T, data = pred.frame
  
          StanER=sqrt(diag(vcov( Model_predict[[i]]))) 
        OR_95[i,]=cbind(Model_predict[[i]]$coefficients[[feats_of_interest[i]]], exp(Model_predict[[i]]$coefficients[[feats_of_interest[i]]]),
                 exp(Model_predict[[i]]$coefficients[[feats_of_interest[i]]]-StanER[[feats_of_interest[i]]]*1.96),
                 exp(Model_predict[[i]]$coefficients[[feats_of_interest[i]]]+StanER[[feats_of_interest[i]]]*1.96))
  
  if (i>1) {
    stats.s=cbind(stats.s,c(Model_predict[[i]]$stats,anova(Model_predict[[i]])[1,3]))
  } else {stats.s=c(Model_predict[[i]]$stats,anova(Model_predict[[i]])[1,3])}
}
colnames(stats.s)<- feats_of_interest;rownames(stats.s)[19]="uP"
univar_uitk2=cbind(t(stats.s[c(19,6,10),order(stats.s[19,])]),OR_95[order(stats.s[19,]),]); colnames(univar_uitk2)[4:7]=c("beta","OR","OR95-","OR95+")
# print(univar_uitk2)#[univar_uitk2[,'uP']<0.05,])
cat("\n","Univariable significant: ", nrow(univar_uitk2[univar_uitk2[,'uP']<0.05,])," of ", nrow(univar_uitk2))
p005=univar_uitk2[univar_uitk2[,'uP']<0.05,]

# ms=ms[order(stats.s[15,]),];uni=cbind(univar_uitk2,ms)

var.sel1.names=rownames(univar_uitk2[univar_uitk2[,'uP']<1,])
var.sel1.stats=univar_uitk2[univar_uitk2[,'uP']<1,]
z=0
correlatiens=c()#matrix(nrow = length(var.sel1.names),ncol = length(var.sel1.names))
for (i in 1:length(var.sel1.names)) {
  for (j in 1:length(var.sel1.names)) {
    if (i!=j & i<j){
      
      correlat=cor(DATA[,var.sel1.names[i]],DATA[,var.sel1.names[j]])
      if (abs(correlat)>0.8){
        z=z+1
        correlatiens=rbind(correlatiens,c(var.sel1.names[i],var.sel1.stats[i],var.sel1.names[j],var.sel1.stats[j],correlat))
      }
    }
  }
}

correlatiens_dyn=correlatiens
i1=0;preselected_vars=c()

while (size(correlatiens_dyn,1)>1){
i1=i1+1
superior_feature=correlatiens_dyn[1,1]
exluderen=correlatiens_dyn[correlatiens_dyn[,1] %in%  superior_feature,3]
correlatiens_dyn=correlatiens_dyn[!correlatiens_dyn[,1] %in%  exluderen,]
if (size(correlatiens_dyn,1)>1){
correlatiens_dyn=correlatiens_dyn[!correlatiens_dyn[,1] %in%  superior_feature,]}
preselected_vars[i1]=superior_feature
print(size(correlatiens_dyn,1));#
# print(t(t(exluderen)))
}
if (size(correlatiens_dyn,1)==1){
preselected_vars[i1+1]=correlatiens_dyn[1]}
no_corr_at_all=var.sel1.names[!var.sel1.names %in% c(correlatiens[,1],correlatiens[,3])]

# print(preselected_vars)

features_of_interest0=var.sel1.names[!var.sel1.names %in% correlatiens[,3]]
features_of_interest=var.sel1.names[var.sel1.names %in% c(no_corr_at_all,preselected_vars)]


Preselection=var.sel1.stats[var.sel1.names %in% features_of_interest,]
cat("\n","Selected: ", nrow(Preselection)," of ", nrow(univar_uitk2))
cat("\n")
print(round(Preselection[Preselection[,1]<0.05,],3))

return(Preselection)
}

# outputdir="//zkh/dfs/Gebruikers26/DijkLV/Data/AA_PROJECTS_DATA/A3_PET_IBMS/univariable_analysis"
# write.xlsx(Preselection, paste(outputdir,"TEX_preSelection_XER_OR_Univariable_intensity.xlsx",sep="/"))
# write.xlsx(correlatiens, paste(outputdir,"TEX_correlatiens.xlsx",sep="/"))
# write.xlsx(univar_uitk2, paste(outputdir,"TEX_univar_intensity.xlsx",sep="/"))
