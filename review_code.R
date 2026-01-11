### Two-dimensional data examination by exploratory functional data analysis 
### to improve detection of scattered soil contamination by Cu-bearing pesticides

### TMG, ST, MH, KR, IP, SS, JM

### supplementary code
################################################################################

#########
### initial setting
#########
setwd("C:/Users/Ivana/Desktop/OneDrive - Univerzita Palackého v Olomouci/p h d/Karel/Densities_new2/Densities/Standa_05_2025")

library(class)
library(cluster)
library(fda)
library(gplots)
library(robustbase)
library(splines2)


# abbreviations for all districts
abbr77=c("BN", "BE", "BK", "BV", "BM", "BI", "BR", "CL", "CB", "CK", "DC", "DO", "FM", "HB", "HO",
         "HK", "CH", "CV", "CR", "JN", "JE", "JC", "JI", "JH", "KV", "KI", "KD", "KT", "KO", "KM",
         "KH", "LI", "LT", "LN", "ME", "MB", "MO", "NA", "NJ", "NB", "OC", "OP", "OV", "PU", "PE",
         "PI", "PJ", "PM", "PS", "A",  "PY", "PZ", "PT", "PR", "PB", "PV", "RA", "RO", "RK", "SM",
         "SO", "ST", "SU", "SY", "TA", "TC", "TP", "TR", "TU", "UH", "UL", "UO", "VS", "VY", "ZR", 
         "ZL", "ZN")

###############
### FIGURE 4

b_coef_Cu_ar=as.matrix(read.csv2("Cu_arithmetic_new.csv",header=F))
range_Cu=c(0.588,4.593)

basis_Cu=create.bspline.basis(rangeval = range_Cu,breaks=seq(range_Cu[1],range_Cu[2],length=5+2),norder=3)
fd_data_Cu_ar=fd(t(b_coef_Cu_ar),basisobj = basis_Cu)
plot(basis_Cu)
plot.fd(fd_data_Cu_ar)

t.fine=seq(range_Cu[1],range_Cu[2],length=1000)

#
Cu_ar=t(eval.fd(t.fine,fd_data_Cu_ar))
rownames(Cu_ar)=t(abbr77)

rgb.palette <- colorRampPalette(c("blue4","turquoise","white","orange","red4"),
                                space = "rgb")
quant=quantile(c(Cu_ar),probs = seq(0,1,0.01)) #color breaks according to quantiles

#Fig. 4
HM_fig4=heatmap.2(Cu_ar,col = rgb.palette(100), key=FALSE,keysize=0.8,density.info="none", symkey=FALSE, 
          breaks=quant, trace="none",cexCol=0.8,cexRow=0.4, scale = "none",
          Colv=FALSE,main="Cu - arithmetic marginals", dendrogram="row",
          labCol=round(exp(t.fine),1),xlab="Cu(mg kg-1)",labRow = abbr77)

rm(list = setdiff(ls(), c("abbr77","range_Cu","t.fine","basis_Cu")))
dev.off()
###############
### FIGURE 5

setwd("C:/Users/Ivana/Desktop/OneDrive - Univerzita Palackého v Olomouci/p h d/Karel/Densities_new2/Densities/Standa_05_2025")

overall_coeffs=list()
{
  overall_coeffs[[1]]=usti_nad_orlici=as.matrix(read.csv("B_overall_usti_nad_orlici.csv",header=F,sep=";",dec=","))
  overall_coeffs[[2]]=rakovnik=as.matrix(read.csv("B_overall_rakovnik.csv",header=F,sep=";",dec=","))
  overall_coeffs[[3]]=kutna_hora=as.matrix(read.csv("B_overall_kutna_hora.csv",header=F,sep=";",dec=","))
  overall_coeffs[[4]]=plzen_sever=as.matrix(read.csv("B_overall_plzen_sever.csv",header=F,sep=";",dec=","))
}

abbr4=c("UO","RA","KH","PS")

# parameters
alpha=0.5
k=l=2

g=h=5

p=q=1

range_Cu=c(0.588,4.593)
range_Zn=c(1.457,5.665)

knots_Cu=seq(range_Cu[1],range_Cu[2],length=7)
knots_Zn=seq(range_Zn[1],range_Zn[2],length=7)

grid_Cu=seq(range_Cu[1],range_Cu[2],length=60)
grid_Zn=seq(range_Zn[1],range_Zn[2],length=60)

B_Cu <- bSpline(grid_Cu, knots = knots_Cu[-c(1,7)], degree = 2, intercept = TRUE)
B_Zn <- bSpline(grid_Zn, knots = knots_Zn[-c(1,7)], degree = 2, intercept = TRUE)

rgb.palette <- colorRampPalette(c("blue4","turquoise","white","orange","red4"),
                                space = "rgb")
par(pty="s",mfrow=c(2,2))
for (ind in 1:4){
  
  Z<- outer(1:length(grid_Cu), 1:length(grid_Zn), Vectorize(function(i, j) {
    sum(B_Cu[i, ] %*% as.matrix(overall_coeffs[[ind]]) %*% as.matrix(B_Zn[j, ]))
  }))
  
  image(exp(grid_Cu),exp(grid_Zn),Z,col=rgb.palette(100),main=paste(abbr4[ind],": overall"),
        xlab="Cu",ylab="Zn",log="xy",cex.axis=0.9)
}
dev.off()

###############
### FIGURE 8
setwd("C:/Users/Ivana/Desktop/OneDrive - Univerzita Palackého v Olomouci/p h d/Karel/Densities_new2/Densities/regrese_06_2025")

data=read.csv2("Initially_proc_data.csv",header=T,dec=",",sep=";")
data[,1]=as.factor(data[,1])
districts=levels(as.factor(data[,1]))

###
reg=lmrob(CU_HNO3~ZN_HNO3,data=data,method="MM")
reg

range(reg$residuals)

hist1=hist(reg$residuals,breaks=50,probability = T)
breaks=seq(-1.9,2.7,by=0.1)
counts <- table(cut(reg$residuals, breaks = breaks, right = FALSE))

#loop for all districts - proportions
districts_residuals_prop=matrix(NA,nrow=length(districts),ncol=length(breaks)-1)
for(i in 1:length(districts)){
  districts_residuals_prop[i,]=c(table(cut(reg$residuals[data$NAZEV==districts[i]], breaks = breaks, right = FALSE)))/length(reg$residuals[data$NAZEV==districts[i]])
}
rownames(districts_residuals_prop)=abbr77

rgb.palette <- colorRampPalette(c("blue4","turquoise","white","orange","red4"),
                                space = "rgb")
quant=quantile(c(districts_residuals_prop)[districts_residuals_prop>0],probs = seq(0,1,0.01))

#Fig.8
heatmap.2(districts_residuals_prop,col = rgb.palette(100), key=T,keysize=1,density.info="none", symkey=FALSE,key.title=NA,key.xlab=NA,key.ylab = NA,
          breaks=quant, trace="none",cexCol=0.6,cexRow=0.4, scale = "none",
          Colv=FALSE,main="", dendrogram="row",
          labCol=round(breaks,1),xlab="RRes ln(Cu)",labRow = abbr77,cex.main=0.5,cex.lab=0.5)

rm(list = setdiff(ls(), c("abbr77","range_Cu","t.fine","basis_Cu")))
dev.off()

###############
### FIGURE 9
setwd("C:/Users/Ivana/Desktop/OneDrive - Univerzita Palackého v Olomouci/p h d/Karel/Densities_new2/Densities/Standa_05_2025")
b_coef_Cu=as.matrix(read.csv2("B_Cu_geometric.csv",header=F))

rownames(b_coef_Cu)=t(abbr77)

basis_Cu=create.bspline.basis(rangeval = range_Cu,breaks=seq(range_Cu[1],range_Cu[2],length=5+2),norder=3)
fd_data_Cu_ge=fd(t(b_coef_Cu),basisobj = basis_Cu)

plot(basis_Cu)
plot.fd(fd_data_Cu_ge)

#
Cu_ge=t(eval.fd(t.fine,fd_data_Cu_ge))
rgb.palette <- colorRampPalette(c("blue4","turquoise","white","orange","red4"),
                                space = "rgb")
quant=quantile(c(Cu_ge),probs = seq(0,1,0.01)) #color breaks according to quantiles
HM_fig9=heatmap.2(Cu_ge,col = rgb.palette(100), key=T,keysize=1,density.info="none", symkey=FALSE,key.title=NA,key.xlab=NA,key.ylab = NA,
          breaks=quant, trace="none",cexCol=0.8,cexRow=0.4, scale = "none",
          Colv=FALSE,main="Cu - geometric marginals", dendrogram="row",
          labCol=round(exp(t.fine),1),xlab="Cu(mg kg-1)",labRow = abbr77)

