#---------------Clean work area, from the beginning------------------
rm(list = ls()) #Remove all objects
graphics.off()  #Remove all graphics
cat("\014")     #Remove script in windows console

## Load relevant libraries
library(chron)
library(grDevices)
library(optimx)
library(plantecophys)
library(DEoptim)
library(weathermetrics)
library(minpack.lm)
library(readxl)


#------------------------------------------ 
#Change the directory where you want to save your data
#setwd("/Users/maquellegarcia/Documents/GitHub/DAT_Proj/Results")

# Give a name to the files (table and pdf) that will receive the results
arquivo <-"A_ci_fit_DAT_Tapajos_20230106"
# creates the pdf file that will receive the graphs
  pdf(file=paste(arquivo, ".pdf", sep=""),height=10,width=20)
  
# Creating a matrix that stores the column names of the text file that will receive the results of this script
  output.names     <- matrix (c("Tree_id","DAT","person","Mean.Tleaf","Max.Tleaf",
                                "Asat", "Amax","Amin", "Ci.at.Asat","E.at.Asat","Ca","A/Ci","Ci/Ca",
                                "gs.at.Asat","gs.index","Mean.gs","mean.RH","Flow",
                                "mean.Press","Max.VPD","Ciindex",
                                
                                "Vcmax_1ptR1",
                                "Vcmax_1ptR2",
                                "Vcmax_1ptR3",
                                
                                "vcmax_Best_Model","Jmax_Best","TPU_Best","Rd_Best","p5","sequencia","erros",
                                #"Vcmax_PlantEcophys", "Jmax_PlantEcophys", "TPU_PlantEcophys", "Rd_PlantEcophys", "Error_PlantEcophys",
                                
                                "Vcmax_One_Point_Method_25C_1",
                                "Vcmax_One_Point_Method_25C_2",
                                "Vcmax_One_Point_Method_25C_3",
                                "Best_Vcmax_25C",
                                "Best Jmax_25C",
                                
                                #"Vcmax_PlantEcophys_25C",  
                                #"Jmax_PlantEcophys_25C",
                                #"Vcmax_PlantEcophys_25C_R", 
                                #"Jmax_PlantEcophys_25C_R",
                                "Best_Vcmax_25C_Bernac",
                                "Best Jmax_25C_Bernac",
                                
                                "GammaStar",
                                "Kmi", 
                                "TleafK.Asat",
                                "Gma2",
                                "Km2"),nrow=1)

  colnames (output.names)	<- output.names
  write.table (output.names, paste(arquivo, ".csv", sep=""), append=TRUE, sep=",",
             row.names=FALSE, col.names=FALSE)
  

## Function for fitting Farquhar's biochemical model based on A-Ci curves 
## Error estimated as the sum of squared differences between modeled and observed assimilation rates
## there is a 'penalty function' to keep the switch from Av to Aj close to Ci = 0.7 x Ca (or 28 Pa)
A_Ci_Fit1  <- function(x){
  Vcmax <- x[1]    # variable 1
  Jmax  <- x[2]    # variable 2
  TPU   <- x[3]		 # variable 3
  Rd    <- x[4]		 # variable 4
  
  J    <- (Qabs + Jmax - sqrt((Qabs + Jmax)^2 - 4 * curv * Qabs * Jmax))/(2*curv) # electron transport rate when light is not saturating
  aj   <- ifelse(Curve$Pci>=10, J * (Curve$Pci - gstari_Tleaf) / (4*Curve$Pci + 8 * gstari_Tleaf)-Rd,999)
  av   <- ifelse(Curve$Pci<50, Vcmax * (Curve$Pci - gstari_Tleaf) / (Kmi + Curve$Pci)-Rd, 999)# carboxylation rate
  atpu <- ifelse (Curve$Pci == tail(Curve$Pci,1),3 * TPU -Rd , 999)	          # triose phosphate utilization
  if (tail(Curve$A,1) <= Curve$A[length(Curve$A)-1]){
    aj   <- ifelse(Curve$Pci>=10 & Curve$Pci< max(Curve$Pci),J * (Curve$Pci - gstari_Tleaf) / (4*Curve$Pci + 8 * gstari_Tleaf)-Rd,999)
                  }  
  a     <- pmin(av,aj)
  adiff <- sum((Curve$A - a)^2)				 # least-square
  bdiff <- sum((Curve$A[Curve$Pci == tail(Curve$Pci,1)] - atpu[Curve$Pci == tail(Curve$Pci,1)])^2)
  #diffsum <- adiff + bdiff +  abs(28-Cseq[max(head(sort((av-aj)^2, index.return=TRUE)$ix,2))]) + (Rd*(-1))#essa parte adiciona a penalizacao
}

## Function for fitting the A-Ci model $$ Error estimated as mean squared differences between modeled and observed assimilation rates
## It uses calculated absorbed PAR, instead of Incident PAR
A_Ci_Fit2  <- function(x){
  Vcmax <- x[1]    # variable 1
  Jmax  <- x[2]    # variable 2
  TPU   <- x[3]		 # variable 3
  Rd    <- x[4]		 # variable 4
  
  J    <- (Qabs + Jmax - sqrt((Qabs + Jmax)^2 - 4 * curv * Qabs * Jmax))/(2*curv) # electron transport rate when light is not saturating
  aj   <- ifelse(Curve$Pci>=10, J * (Curve$Pci - gstari_Tleaf) / (4*Curve$Pci + 8 * gstari_Tleaf)-Rd,999)
  av   <- ifelse(Curve$Pci<50, Vcmax * (Curve$Pci - gstari_Tleaf) / (Kmi + Curve$Pci)-Rd, 999)# carboxylation rate
  atpu <- ifelse (Curve$Pci == tail(Curve$Pci,1),3 * TPU -Rd , 999)	          # triose phosphate utilization
  if (tail(Curve$A,1) <= Curve$A[length(Curve$A)-1]){
    aj   <- ifelse(Curve$Pci>=10 & Curve$Pci< max(Curve$Pci),J * (Curve$Pci - gstari_Tleaf) / (4*Curve$Pci + 8 * gstari_Tleaf)-Rd,999)
                }  
  a <- pmin(av,aj)
  adiff <- sum((Curve$A - a)^2)	# least-square
  bdiff <- sum((Curve$A[Curve$Pci == tail(Curve$Pci,1)] - atpu[Curve$Pci == tail(Curve$Pci,1)])^2)
  #diffsum <- adiff + bdiff + abs((Jmax/Vcmax) - (-0.07* mean(Curve$Tleaf) + 3.8)) +    (Rd*(-1)) #
  # correao pela temperatura +    (Rd*(-1)) # + abs((aj/av) - (-0.07* mean(Curve$Tleaf) + 3.8))
}


## Same as A-Ci_Fit1, but using incident light (Qin) instead of absorbed light (Qin).
A_Ci_Fit3  <- function(x){
  Vcmax <- x[1]    # variable 1
  Jmax  <- x[2]    # variable 2
  TPU   <- x[3]		 # variable 3
  Rd    <- x[4]		 # variable 4
 curv2  <- x[5]
 
  J    <- (Qabs + Jmax - sqrt((Qabs + Jmax)^2 - 4 * curv2 * Qabs * Jmax))/(2*curv2) # electron transport rate when light is not saturating
  aj   <- ifelse(Curve$Pci>=10, J * (Curve$Pci - gstari_Tleaf) / (4*Curve$Pci + 8 * gstari_Tleaf)-Rd,999)
  av   <- ifelse(Curve$Pci<50, Vcmax * (Curve$Pci - gstari_Tleaf) / (Kmi + Curve$Pci)-Rd, 999)				            # carboxylation rate
  atpu <- ifelse (Curve$Pci == tail(Curve$Pci,1),3 * TPU - Rd , 999)	          # triose phosphate utilization
  if (tail(Curve$A,1) <= Curve$A[length(Curve$A)-1]){
   aj   <- ifelse(Curve$Pci>=10 & Curve$Pci< max(Curve$Pci),J * (Curve$Pci - gstari_Tleaf) / (4*Curve$Pci + 8 * gstari_Tleaf)- Rd,999)
            }
  a     <- pmin(av,aj)
  adiff <- sum((Curve$A - a)^2)# least-square
  bdiff <- sum((Curve$A[Curve$Pci == tail(Curve$Pci,1)] - atpu[Curve$Pci == tail(Curve$Pci,1)])^2)
  #diffsum <- adiff + bdiff +  abs(28-Cseq[max(head(sort((av-aj)^2, index.return=TRUE)$ix,2))]) + (Rd*(-1))#penalizacao
}

#### Função para conversão parâmetros para temperatura padrão 25°C, de acordo com Bernacchi
Temp_Correction <- function (p1,p2,p3,p4) {
  a <- p1/(exp(26.355 - (65.33 /(R*meanTleafK))))
  b <- p2/(exp( 17.71 - (43.9  /(R*meanTleafK))))
  c <- p3/((exp(21.46 - 53.1   /(R*meanTleafK)))/(1+exp((0.65*meanTleafK-201.8)/(R*meanTleafK))))
  d <- p4/(exp(18.715 - (46.39 /(R*meanTleafK))))
  return(c(a,b,c,d))
}
### Peaked function Medlin et al. 2002PCE
### Ea = Activation Energy, DeltaS = entropy factor = 200, Hd = Deactivation Energy)
### updated to Kumarathunge et al 2019, New Phytologist.(222: 768-784) doi: 10.1111/nph.15668
Ea_V     <- 82620.87 #58550 or 82620.87 or
Delta_V  <- 645.1013   #629.26 or 645.1013 or
Ed       <- 200000  
b_V <- 1+exp((298.15 * Delta_V - Ed) / (298.15  *R2))

Ea_J     <- 39676.89 # or 39676.89 
Delta_J  <- 641.3615 # or 641.3615
#EdVJ    <- 200000
b_J <- 1 + exp((298.15 * Delta_J - Ed) / (298.15 * R2))

Vc_peaked_25C <- function (Vc_LeafTemp,Ea_V,Delta_V,Ed) {
  a <- exp(Ea_V * (meanTleafK - 298.15) / (298.15 * R2 * meanTleafK))
  c <- 1+exp((meanTleafK * Delta_V - Ed) / (meanTleafK * R2))
  Vc_25 <- Vc_LeafTemp / (a * (b_V / c))
}

J_peaked_25C <- function (J_LeafTemp,Ea_J,Delta_J,Ed) {
  a <- exp(Ea_J * (meanTleafK - 298.15)/(298.15 * R2 * meanTleafK))
  c <- 1 + exp((meanTleafK * Delta_J - Ed) / (meanTleafK * R2))
  Jmax_25 <- J_LeafTemp / (a * (b_J / c))
}    

### Functions for fitted model
modeled_lines <- function (p1,p2,p3,p4){
  av   <- p1 * (Cseq - mean(gstari_Tleaf))/(mean(Kmi) + Cseq ) - p4 
  J    <- mean((Qabs + p2 - sqrt((Qabs + p2)^2 - 4 * curv * Qabs * p2))/(2*curv) )
  aj   <- J * (Cseq - mean(gstari_Tleaf))/(4*Cseq + 8*mean(gstari_Tleaf)) - p4 
  atpu <- p3 * 3 - p4 + (Cseq /10^10) 
  amin <- pmin(av,aj,atpu)
  return(cbind(av,aj,atpu,amin))
}
modeled_points <- function (p1,p2,p3,p4){
  avb   <- p1 * (Curve$Pci - mean(gstari_Tleaf))/(mean(Kmi) + Curve$Pci ) - p4
  J    <- mean((Qabs + p2 - sqrt((Qabs + p2)^2 - 4 * curv * Qabs * p2))/(2*curv) )
  ajb   <- J * (Curve$Pci - mean(gstari_Tleaf)) / (4*Curve$Pci + 8 * mean(gstari_Tleaf)) - p4
  atpub <- p3 * 3 - p4 + (Curve$Pci /10^10)
  aminb <- pmin(avb,ajb,atpub)
  return(cbind(avb,ajb,atpub,aminb))
}  

### Functions for fitted model
modeled_lines2 <- function (p1,p2,p3,p4,p5){
  av   <- p1 * (Cseq - mean(gstari_Tleaf))/(mean(Kmi) + Cseq ) - p4 
  J    <- mean((Qabs + p2 - sqrt((Qabs + p2)^2 - 4 * p5 * Qabs * p2))/(2*p5) )
  aj   <- J * (Cseq - mean(gstari_Tleaf))/(4*Cseq + 8*mean(gstari_Tleaf)) - p4 
  atpu <- p3 * 3 - p4 + (Cseq /10^10) 
  amin <- pmin(av,aj,atpu)
  return(cbind(av,aj,atpu,amin))
}

modeled_points2 <- function (p1,p2,p3,p4,p5){
  avb   <- p1 * (Curve$Pci - mean(gstari_Tleaf))/(mean(Kmi) + Curve$Pci ) - p4
  J    <- mean((Qabs + p2 - sqrt((Qabs + p2)^2 - 4 * p5 * Qabs * p2))/(2*p5) )
  ajb   <- J * (Curve$Pci - mean(gstari_Tleaf)) / (4*Curve$Pci + 8 * mean(gstari_Tleaf)) - p4
  atpub <- p3 * 3 - p4 + (Curve$Pci /10^10)
  aminb <- pmin(avb,ajb,atpub)
  return(cbind(avb,ajb,atpub,aminb))
}  

modeled_lines3 <- function (par1,par2,par3,par4){
    av   <- par1 * (Cseq - mean(gstari_Tleaf))/(mean(Kmi) + Cseq ) - par4 
    J    <- mean((Qabs + par2 - sqrt((Qabs + par2)^2 - 4 * curv * Qabs * par2))/(2*curv) )
    aj   <- J * (Cseq - mean(gstari_Tleaf))/(4*Cseq + 8*mean(gstari_Tleaf)) - par4 
    atpu <- par3 * 3 - par4 + (Cseq /10^10) 
    amin <- pmin(av,aj,atpu)
    return(cbind(av,aj,atpu,amin))
}


modeled_points3 <- function (par1,par2,par3,par4){
    avb   <- par1 * (Curve$Pci - mean(gstari_Tleaf))/(mean(Kmi) + Curve$Pci ) - par4
    J    <- mean((Qabs + par2 - sqrt((Qabs + par2)^2 - 4 * curv * Qabs * par2))/(2*curv) )
    ajb   <- J * (Curve$Pci - mean(gstari_Tleaf)) / (4*Curve$Pci + 8 * mean(gstari_Tleaf)) - par4
    atpub <- par3 * 3 - par4 + (Curve$Pci /10^10)
    aminb <- pmin(avb,ajb,atpub)
    return(cbind(avb,ajb,atpub,aminb))
}  

#-----------------------------------------------
setwd("/Users/maquellegarcia/Documents/GitHub/DAT_Proj/Inputs")
dir()
curvas<-read_excel("Aci_no_out.xlsx")

curvas2<-subset(curvas, Data_QC=="OK")#to exclude weird points 
curvas1<-subset(curvas, Ci>0)#to avoid negative values 

names(curvas)
  
  sp<-as.data.frame(unique(curvas1[,"unique_id"]))
  colnames(sp)<-"sp"
  sp
  

  #curve_names<-NULL
  
  for (i in 1:length(sp[,1])) {
  Curve<- subset(curvas1, unique_id==sp[i,1])
  #Curve<-subset(Curve,excluir<1)
  

  names(Curve)
  Curve$Tleaf <- as.numeric(Curve$Tleaf)
  
  ###Constants used in the Farquhar´s model ***NOT*** considering mesophyll conductance
  R             <- 0.008314    # Gas constant
  R2            <- R * 1000
  Kc   		      <- 40.49		   # Michaelis-Menten constant for CO2 (Pa) (Bernacchi et al 2001,2002) or 404 microbar von Caemmerer et al. (1994)
  delta_Kc 		  <- 79430 	     # (J mol-1) from Medlyn et al 2002
  Ko   		      <- 27.84		   # Michaelis-Menten constant for CO2 (kPa)(Bernacchi et al 2001,2002) or 248 mbar von Caemmerer et al. (1994)
  delta_Ko 		  <- 36380	     # (J mol-1) from Medlyn et al 2002
  gstar  		    <- 4.275		   # CO2 compensation point (Pa) (Bernacchi et al 2001,2002) or 36.9 microbar von Caemmerer et al. (1994)
  delta_gstar 	<- 37830	     # (J mol-1)
  
  #Other Constants
  curv          <- 0.85        # Curvature for calculating J from Evans (1989)
  Ambient.CO2   <- 400
  O2  	      	<- 21 		     # Estimated Oxygen concentration at chloroplast - (kPa)
  K0            <- 0.388       # for accounting for diffusion of CO2 through the gasket (Licor 6400 only)
  
  ## CO2 sequence (will be needed for building the figures)
  Cseq <- seq(0,200,0.5)
  
  # converts Ci from ppm to Pa
  Curve$Pci      <- Curve$Ci * Curve$Pa * 0.001


  Tleafk <- Curve$Tleaf +273.15
  Data <-format(as.POSIXct(strptime(as.character(Curve$date), format = "%Y%m%e %H:%M:%S", tz = "GMT")), "%Y%m%d")
  Data2<-format(as.POSIXct(strptime(as.character(Curve$date), format = "%Y%m%e %H:%M:%S", tz = "GMT")), "%M:%S")
  
# A standardized measure of how much stomatal conductance changed during the measurements
  gsindex <- (max(Curve$gsw) - min(Curve$gsw)) / max(Curve$gsw)
  Ciindex <- (max(Curve$Ci)- min(Curve$Ci))/max(Curve$Ci)

# Enzymatic kinects at actual leaf temperature
  Kci_Tleaf	  	<- Kc    * exp((Tleafk-298.15)*delta_Kc/   (298.15*8.314*Tleafk))
  Koi_Tleaf	   	<- Ko    * exp((Tleafk-298.15)*delta_Ko/   (298.15*8.314*Tleafk))
  gstari_Tleaf	<- gstar * exp((Tleafk-298.15)*delta_gstar/(298.15*8.314*Tleafk))
  Kmi			      <- Kci_Tleaf*(1+O2/Koi_Tleaf)
 # VON CAEMMERER, S. (2013), Steady state models of photosynthesis. Plant Cell Environ, 36: 1617-1630. doi:10.1111/pce.12098
 # Light absorbed = Incident light * 0.85(1-0.15)/2
  # PlantEcophys uses 0.24
 Qabs <- Curve$Qin * 0.85*(1 - 0.15)/2
 Curve$Qabs <- Curve$Qin * 0.85*(1 - 0.15)/2
  
#Some variable to be exported to the final results file
#During a A-Ci curve, there are usually two measurements at 400 ppm CO2 (eg. 400,300,250,200,150,100,400,600,800,1200,1500,2000)
#sometimes the first one is better, sometimes the second one is better
#this chooses the best one, based on the expected Ci/Ca ratio of 0.7
  Curve$Ci_Ca <- Curve$Ci/Curve$CO2_s
  Curve$WUE   <- Curve$A/Curve$gsw
  Curve$CUE   <- Curve$A/Curve$Ci
  Curve$VUE   <- Curve$gsw/Curve$Ci

  # CRef <- if(isTRUE(length(Curve$Ci_Ca[Curve$CO2_r<450 & Curve$CO2_r>350]) >= 1)) which(abs(Curve$Ci_Ca-0.7) == min(abs(Curve$Ci_Ca[Curve$CO2_r<450 & Curve$CO2_r>350]-0.7)))  else NA
  # CRef <- if(isTRUE(length(Curve$VUE[Curve$CO2_r<450 & Curve$CO2_r>350]) >= 1)) which(Curve$VUE == max(Curve$VUE[Curve$CO2_r<450 & Curve$CO2_r>350]))  else NA


  find.best <- try(nlsLM(VUE ~ a*Ci^-b,control=nls.control(maxiter = 100,warnOnly=T)  ,data = Curve,start= c(a=.2,b=1.5)))
  #Res_M1 <- Curve$VUE - (coef(find.best1)[1]*Curve$Ci^-coef(find.best1)[2])
  if (is(find.best, "try-error")) {
    CRef <- if(isTRUE(length(Curve$VUE[Curve$CO2_r<450 & Curve$CO2_r>350]) >= 1)) which(Curve$VUE == max(Curve$VUE[Curve$CO2_r<450 & Curve$CO2_r>350]))  else NA
  } else {
    CRef <- if(isTRUE(length(Curve$VUE[Curve$CO2_r<450 & Curve$CO2_r>350]) >= 1)) which(abs(resid(find.best)) == min(abs(resid(find.best)[Curve$CO2_r<450 & Curve$CO2_r>350])))  else NA
  }
    
  ASat              <- Curve$A[CRef][1]
  Amax              <- max(Curve$A, na.rm=TRUE)  
  Ci.Ca             <- Curve$Ci[CRef][1]/Curve$CO2_s[CRef][1]
  Ci.Asat           <- Curve$Ci[CRef][1]
  Ca.Asat           <- Curve$CO2_s[CRef][1]
  gs.Asat           <- Curve$gsw[CRef][1]
  gstari_Tleaf.Asat <- gstari_Tleaf[CRef][1]
  Kmi.Asat          <- Kmi[CRef][1]
  TleafK.Asat       <- Tleafk[CRef][1] 
  E.Asat            <- Curve$E[CRef][1]

  #
  # Starting values for the "optim" function, based on relationships between Asat or Amax  and Vcmax, Jmax, TPU, and Rd
  # these are just coefficients from fitted linear regressions derived from previous fittings of TROBIT curves
  MeanA <- mean(Curve$A)
  if (is.na(ASat)){
  start.values1 <- c(MeanA * 7.4318+16.927, MeanA * 6.0787+58.343, MeanA * 0.4158+4.2421, 0.8)
  start.values3 <- c(MeanA * 7.4318+16.927, MeanA * 6.0787+58.343, MeanA * 0.4158+4.2421, 0.8,0.85)
  } else{ 
  start.values1 <- c(ASat * 7.4318+16.927, ASat * 6.0787+58.343, ASat * 0.4158+4.2421, 0.8)
  start.values3 <- c(ASat * 7.4318+16.927, ASat * 6.0787+58.343, ASat * 0.4158+4.2421, 0.8,0.85)
  }
 MaxA <- max(Curve$A) 
  if (is.infinite(Amax)){
  start.values2 <- c(MaxA * 3.9821+2,      MaxA * 4.6045+ 7,     MaxA * 0.3281+0.1826, 0.8)
  start.values4 <- c(MaxA * 3.9821+2,      MaxA * 4.6045+ 7,     MaxA * 0.3281+0.1826, 0.8,0.85)
  } else{
    start.values2 <- c(Amax * 3.9821+2,      Amax * 4.6045+ 7,     Amax * 0.3281+0.1826, 0.8)
    start.values4 <- c(Amax * 3.9821+2,      Amax * 4.6045+ 7,     Amax * 0.3281+0.1826, 0.8,0.85)
  }

 # parameters (Vcmax, Jmax, TPU and Rd) minimum and maximum values for the "DEoptim" function
 lower.bound  <- c(-15, -15, -15,-15)     
 upper.bound  <- c(500, 500, 200, 15)
 lower.bound2 <- c(-15, -15, -15,-15, 0.80)     
 upper.bound2 <- c(500, 500, 200, 15, 0.90) 
 

   
# Fitting the A-ci model  using the "optimx" minimization algorithms
  first    <-  subset(optimx(start.values1, A_Ci_Fit1,control=list(all.methods=TRUE)), convcode < 1)
  second   <-  subset(optimx(start.values1, A_Ci_Fit2,control=list(all.methods=TRUE)), convcode < 1)
  third    <-  subset(optimx(start.values3, lower = lower.bound2, upper = upper.bound2,
                             A_Ci_Fit3,control=list(all.methods=TRUE)), convcode <= 1)
  
 fourth   <-  subset(optimx(start.values2, A_Ci_Fit1,control=list(all.methods=TRUE)), convcode < 1)
 fifth    <-  subset(optimx(start.values2, A_Ci_Fit2,control=list(all.methods=TRUE)), convcode < 1)
 sixth    <-  subset(optimx(start.values4, lower = lower.bound2, upper = upper.bound2,
                             A_Ci_Fit3,control=list(all.methods=TRUE)), convcode <= 1)
  

# organizing results from the curve fitting
  fits1     <- summary(if(sum( first[4]<0) < length( first$p1)) subset( first, p4 > 0 ) else  first, order= c(value))[1,c(1:5)]
  fits2     <- summary(if(sum(second[4]<0) < length(second$p1)) subset(second, p4 > 0 ) else second, order= c(value))[1,c(1:5)]
  fits3     <- summary(if(sum( third[4]<0) < length (third$p1)) subset( third, p4 > 0 ) else  third, order= c(value))[1,c(1:5)]
  fits4     <- summary(if(sum(fourth[4]<0) < length(fourth$p1)) subset(fourth, p4 > 0 ) else fourth, order= c(value))[1,c(1:5)]
  fits5     <- summary(if(sum( fifth[4]<0) < length( fifth$p1)) subset( fifth, p4 > 0 ) else  fifth, order= c(value))[1,c(1:5)]
  fits6     <- summary(if(sum( sixth[4]<0) < length( sixth$p1)) subset (sixth, p4 > 0 ) else  sixth, order= c(value))[1,c(1:5)]
  #fits10    <- c(coef(tenth )[-3],coef(tenth)[3])
  #fits10b   <- c(coef(tenthb )[-3],coef(tenthb)[3])
  
  #colnames(df) <- paste(colnames(df),"new",sep="_")
  fits1_algorithim  <- row.names(fits1)
  fits2_algorithim  <- row.names(fits2)
  fits3_algorithim  <- row.names(fits3)
  fits4_algorithim  <- row.names(fits4)
  fits5_algorithim  <- row.names(fits5)
  fits6_algorithim  <- row.names(fits6)
  #fits10_algorithim <- "PlantEcophys"

  fits1 <- unlist(fits1)
  fits2 <- unlist(fits2)
  fits3 <- unlist(fits3)
  fits4 <- unlist(fits4)
  fits5 <- unlist(fits5)
  fits6 <- unlist(fits6)
  
  #################### calculating rates from the curve fitting results in order to build figures
  modeled_1 <- do.call(modeled_points,as.list(fits1[1:4]))
  fits1_erro <- sqrt(sum((Curve$A - modeled_1[,4])^2))/length(Curve$A)
  result_lines_1 <- do.call(modeled_lines,as.list(fits1[1:4]))
  
  modeled_2  <- do.call(modeled_points,as.list(fits2[1:4]))
  fits2_erro <- sqrt(sum((Curve$A - modeled_2[,4])^2))/length(Curve$A)
  result_lines_2 <- do.call(modeled_lines,as.list(fits2[1:4]))
  
  modeled_3  <- do.call(modeled_points2,as.list(fits3[1:5]))
  fits3_erro <- sqrt(sum((Curve$A - modeled_3[,4])^2))/length(Curve$A)
  result_lines_3 <- do.call(modeled_lines2,as.list(fits3[1:5]))
  
  modeled_4  <- do.call(modeled_points,as.list(fits4[1:4]))
  fits4_erro <- sqrt(sum((Curve$A - modeled_4[,4])^2))/length(Curve$A)
  result_lines_4 <- do.call(modeled_lines,as.list(fits4[1:4]))
  
  modeled_5  <- do.call(modeled_points,as.list(fits5[1:4]))
  fits5_erro <- sqrt(sum((Curve$A - modeled_5[,4])^2))/length(Curve$A)
  result_lines_5 <- do.call(modeled_lines,as.list(fits5[1:4]))
  
  modeled_6  <- do.call(modeled_points,as.list(fits6[1:4]))
  fits6_erro <- sqrt(sum((Curve$A - modeled_6[,4])^2))/length(Curve$A)
  result_lines_6 <- do.call(modeled_lines,as.list(fits6[1:4]))
  

  ###################################
  Trans_fit1 <- max(head(sort((result_lines_1[,1]-result_lines_1[,2])^2, index.return=TRUE)$ix,2))
  Trans_fit2 <- max(head(sort((result_lines_2[,1]-result_lines_2[,2])^2, index.return=TRUE)$ix,2))
  Trans_fit3 <- max(head(sort((result_lines_3[,1]-result_lines_3[,2])^2, index.return=TRUE)$ix,2))
  Trans_fit4 <- max(head(sort((result_lines_4[,1]-result_lines_4[,2])^2, index.return=TRUE)$ix,2))
  Trans_fit5 <- max(head(sort((result_lines_5[,1]-result_lines_5[,2])^2, index.return=TRUE)$ix,2))
  Trans_fit6 <- max(head(sort((result_lines_6[,1]-result_lines_6[,2])^2, index.return=TRUE)$ix,2))
  #################### comparing fit errors in order to find the "Best fit". **** Unfortunately, often the best fit will not return the smallest error.
  erros <- round(c(fits1_erro,fits2_erro,fits3_erro,fits4_erro,
                   fits5_erro,fits6_erro),2)
  df2<- data.frame(rbind(fits1,fits2,fits3,fits4,fits5,fits6))
  #row.names(df2)[9] <- "fits9" 
  colnames(df2)[5] <- "p5"
  df2$sequencia <- c(1,2,3,4,5,6)
  df2$erros <- erros
  df2$Fit <- c("Model_1","Model_2","Model_3","Model_1","Model_2","Model_3")
  df2$StartValue <- c("Start_1","Start_1","Start_1","Start_2","Start_2","Start_2")
  df2$Algorithim <- c(fits1_algorithim,fits2_algorithim,fits3_algorithim,fits4_algorithim,fits5_algorithim,
                      fits6_algorithim)
  df2$Transition <- c(Trans_fit1,Trans_fit2,Trans_fit3,Trans_fit4,Trans_fit5,Trans_fit6)
  
  indice <- order(df2$erros)[!duplicated(sort(df2$erros))]
  
  um   <- which(which.min(erros)==df2$sequencia)
  dois <- ifelse (length(indice)>=2,indice[[2]], "NA")
  tres <- ifelse (length(indice)>=3,indice[[3]], "NA")
  
  best        <- df2[um  ,1:7]
  second_best <- df2[dois,1:4]
  third_best  <- df2[tres,1:4]
  
  
  # Vcmax from One Point Method by De kauwe 2016, Corrigendum (https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.14172).
  Kc2  <- 404.9 * exp(79430*(Tleafk-298.15)/(298.15*8.314*Tleafk))
  Ko2  <- 278.4 * exp(36380*(Tleafk-298.15)/(298.15*8.314*Tleafk))
  Km2  <- Kc2 * (1+210/Ko2)
  Gma2 <- 42.75 * exp(37830*(Tleafk-298.15)/(298.15*8.314*Tleafk))
  
  Vcmax_1pt     <- Curve$A / ((Curve$Ci - Gma2)/(Curve$Ci + Km2) - 0.015) 
  
  ## scaling leaf respiration to leaf temperature from Heskel et al. 2016, PNAS (www.pnas.org/cgi/doi/10.1073/pnas.1520282113) 
  RT <- 0.015 * exp(0.1012 *(Curve$Tleaf - 25) - 0.0005 *(Curve$Tleaf^2 - 25^2))
  Vcmax_1ptR     <- (Curve$A+(0.015*df2[indice[1],4])) * (Curve$Ci + Km2)/(Curve$Ci - Gma2) 
  
  Vcmax_1ptR2     <- (Curve$A+(0.015*df2[indice[1],1])) * (Curve$Ci + Km2)/(Curve$Ci - Gma2) 
  
  Ratio_Vc_Rday <- ((Curve$Ci - Gma2)/(Curve$Ci + Km2)) - (Curve$A/df2[indice[1],1])
  
  
  Ra <- exp(Ea_V * (Tleafk - 298.15) / (298.15 * R2 * Tleafk))
  Rb <- 1+exp((298.15 * Delta_V - Ed) / (298.15  *R2))
  Rc <- 1+exp((Tleafk * Delta_V - Ed) / (Tleafk * R2))
  RV <- 0.015*(Ra*Rb/Rc)
  
  Vcmax_1ptR3     <- Curve$A / ((Curve$Ci - Gma2)/(Curve$Ci + Km2) - (0.015*RV/RT)) ### not sure why it is RV/RT instead of RT/RV!!!
  
  
  ###########################
  Tleaf_k <- Curve$Tleaf +273.15
  meanTleafK    <- mean(Tleafk)
  
  Data <-format(as.POSIXct(strptime(as.character(Curve$date), format = "%Y%m%e %H:%M:%S", tz = "GMT")), "%Y%m%d")
  Data2<-format(as.POSIXct(strptime(as.character(Curve$date), format = "%Y%m%e %H:%M:%S", tz = "GMT")), "%M:%S")
  
  
  # Vcmax from One Point Method by De kauwe 2016, Corrigendum ().
  Kc        <- 404.9*exp(79430*(Tleaf_k-298.15)/(298.15*R2*Tleaf_k))
  Ko        <- 278.4*exp(36380*(Tleaf_k-298.15)/(298.15*R2*Tleaf_k))
  Km        <- Kc*(1+210/Ko)
  Gama_Star <- 42.75*exp(37830*(Tleaf_k-298.15)/(298.15*R2*Tleaf_k))
  
  #PARAMETROS VINDO DAS CURVAS DE CO2 ver as diferenças do Ecophys
  c <-30.4
  Delta_Ha <-75.3
  Delta_S <-0.79 
  Delta_Hd <- 251.3
  
  
  RT <- 0.015 * exp(0.1012 *(Curve$Tleaf - 25) - 0.0005 *(Curve$Tleaf^2 - 25^2))
  RV<-((exp(c - Delta_Ha/(R*Tleaf_k)))/(1+exp((Delta_S*Tleaf_k-Delta_Hd)/(R*Tleaf_k))))
  
  Vcmax_1point     <- Curve$A / ((Curve$Ci -  Gama_Star)/(Curve$Ci + Km) - (0.015*RT/RV)) ### not sure why it is RV/RT instead of RT/RV!!!
  
  Curve$Vcmax_1point <- Vcmax_1point 
  
  
  onept.vcmax1<-Vcmax_1pt[CRef][1]
  onept.vcmax2<-Vcmax_1ptR3[CRef][1]
  onept.vcmax3<-Vcmax_1point[CRef][1]
  
  
  Ra_Tleaf <- exp(Ea_V * (Tleafk[CRef][1] - 298.15) / (298.15 * R2 * Tleafk[CRef][1]))
  Rc_Tleaf <- 1+exp((Tleafk[CRef][1] * Delta_V - Ed) / (Tleafk[CRef][1] * R2))
  
    
  OnePoint_Vcmax_25C_1     <- Vcmax_1pt[CRef][1]/(Ra_Tleaf*Rb/Rc_Tleaf)#C02 400
  OnePoint_Vcmax_25C_2     <- Vcmax_1point[CRef][1]/(Ra_Tleaf*Rb/Rc_Tleaf)
  OnePoint_Vcmax_25C_ptR3  <- Vcmax_1ptR3[CRef][1]/(Ra_Tleaf*Rb/Rc_Tleaf)  

  
  
##### temperature conversions to a reference temperature    

  best_Vcmax_25C        <- Vc_peaked_25C(       best[[1]],Ea_V, Delta_V, Ed)
  second_best_Vcmax_25C <- Vc_peaked_25C(second_best[[1]],Ea_V, Delta_V, Ed)
  third_best_Vcmax_25C  <- Vc_peaked_25C( third_best[[1]],Ea_V, Delta_V, Ed)
  #Ecophys_Vcmax_25C     <- Vc_peaked_25C(     fits10[[1]],Ea_V, Delta_V, Ed)
  #Ecophys_Vcmax_25C_R   <- fits10b[[1]]
  
  best_Jmax_25C         <- J_peaked_25C(       best[[2]], Ea_J, Delta_J, Ed)
  second_best_Jmax_25C  <- J_peaked_25C(second_best[[2]], Ea_J, Delta_J, Ed)
  third_best_Jmax_25C   <- J_peaked_25C( third_best[[2]], Ea_J, Delta_J, Ed) 
  #Ecophys_Jmax_25C      <- J_peaked_25C(     fits10[[2]], Ea_J, Delta_J, Ed)
  #Ecophys_Jmax_25C_R    <- fits10b[[2]]
  ###################################

  best_Vcmax_25C_B        <- do.call(Temp_Correction, as.list(get(row.names(df2)[indice[1]])[1:4]))


  result_lines <- if (um[1] %in% c(3,6,9)){
                         as.data.frame(do.call(modeled_lines2, as.list(get(row.names(df2)[um[1]])[1:5])))
                         }else{
                         as.data.frame(do.call(modeled_lines, as.list(get(row.names(df2)[um[1]])[1:4])))}
                      
  result_points <- if (um[1] %in% c(3,6,9)){
                          as.data.frame(do.call(modeled_points2, as.list(get(row.names(df2)[um[1]])[1:5])))
                          }else{
                          as.data.frame(do.call(modeled_points, as.list(get(row.names(df2)[um[1]])[1:4])))}

################### making the figures
  inter <- (max(Curve$A)-min(Curve$A))/6
  f1 <-(min(Curve$A))
  f2 <-(min(Curve$A)+ inter)
  f3 <-(min(Curve$A)+ 2 * inter)
  f4 <-(min(Curve$A)+ 3 * inter)
  f5 <-(min(Curve$A)+ 4 * inter)
  f6 <-(min(Curve$A)+ 5 * inter)
  f7 <-(min(Curve$A)+ 6 * inter)
  f8 <-(min(Curve$A)+ 7 * inter)
  
  #par(mfrow=c(3,4),mar=c(4,4,3,4)+0.1,oma=c(0.5,0.5,0.5,1.5))
  par(mfrow=c(1,2),mar=c(4,4,3,4)+0.1,oma=c(0.5,0.5,0.5,1.5))

  ### Plot 1  
  plot (Curve$Pci,Curve$A,col ="blue",pch=19,xlim=c(0,200), cex = 1.5,
        ylim=c(min(Curve$A)-0.5,max(Curve$A)+max(Curve$A)/8), ylab="A", xlab="Ci, Pa")
  title(cex.main = 1.5,main = c (as.character(unique(Curve[,"unique_id"])),paste(df2$Algorithim[um[1]], df2$StartValue[um[1]], df2$Fit[um[1]])))
  lines (Cseq, result_lines$av,   col = "red",   lwd = "4")
  lines (Cseq, result_lines$aj,   col = "sienna2", lwd = "4")
  lines (Cseq, result_lines$atpu, col = "plum", lwd = "4")
  lines (c(27,27),c(-5,45), col = "black", lty = 3)
  lines (c(45,45),c(-5,45), col = "black", lty = 3)
  lines (Cseq, result_lines$amin,   col = "black", 	lwd = "4")
  text (190,f5,paste("Vcmax =",round(df2[indice[1],1],digits=1)),pos=2,col = "red",cex=1.5)
  text (190,f4,paste("Jmax =" ,round(df2[indice[1],2],digits=1)),pos=2,col = "sienna2",cex=1.5)
  text (190,f3,paste("TPU ="  ,round(df2[indice[1],3],digits=1)),pos=2,col = "plum",cex=1.5)
  text (190,f2,paste("Rd ="   ,round(df2[indice[1],4],digits=2)),pos=2,cex=1.5)
  text ( 50,f3,paste("error =",round(df2[indice[1],7],digits=3)),pos=4,col = "red", cex = 1.5)
  text ( 50,f2,paste("gs ="   ,round(mean(Curve$gsw),digits=2)),pos=4,cex=1.5)
  text ( 50,f1,paste("Tleaf =",round(mean(Curve$Tleaf),digits=1)),pos=4,cex=1.5)
  text ( 50,f4,paste("Ci/Ca =",round(Ci.Ca,digits=2)),pos=4,cex=1.5)
  if (df2$Fit[indice[1]] == "Model_3") text (190,f1,paste("Curv =",round(df2[indice[1],5],digits=3)),pos=2,cex=1.5)
  points(Curve$Pci[result_points$avb==result_points$aminb],Curve$A[result_points$avb==result_points$aminb],
         col ="red",pch=19, cex = 2)
  points(Curve$Pci[result_points$ajb==result_points$aminb],Curve$A[result_points$ajb==result_points$aminb],
         col ="sienna2",pch=19, cex = 2)
  points(Curve$Pci[result_points$atpub==result_points$aminb],Curve$A[result_points$atpub==result_points$aminb],
         col ="plum",pch=19, cex = 2)
  #points(Curve$Pci[Curve$Pci > 27 & Curve$Pci < 45],Curve$A[Curve$Pci > 27 & Curve$Pci < 45],
   #      col ="grey", pch=19, cex = 2.5)
  points(Cseq[df2$Transition[indice[1]]],result_lines$av[df2$Transition[indice[1]]],col ="green",pch=9, cex = 2)
  
  colores<- rep(3,length(Curve$A))
  
  
  names(Curve)


  #####################  export analyses to the results table
  output <- matrix (c (unique(Curve[,"unique_id"]),unique(Curve[,"Data_point"]),unique(Curve[,"Scientist"]),mean(Curve$Tleaf),max(Curve$Tleaf),
                       ASat, Amax,min(Curve$A),
                       
                       Ci.Asat,E.Asat,Ca.Asat,(ASat/Ci.Asat),Ci.Ca,gs.Asat,
                       
                       gsindex,mean(Curve$gsw), 
                       
                       mean(Curve$RHcham),mean(Curve$Flow),mean(Curve$Pa),
                       max(Curve$VPDleaf),Ciindex,
                       
                       onept.vcmax1,
                       onept.vcmax2,
                       onept.vcmax3,
                       
                       best,
                  
                       
                       OnePoint_Vcmax_25C_1,
                       OnePoint_Vcmax_25C_2,
                       OnePoint_Vcmax_25C_ptR3,
                       best_Vcmax_25C,
                       best_Jmax_25C,
                       
                       #Ecophys_Vcmax_25C,
                       #Ecophys_Jmax_25C,
                       #Ecophys_Vcmax_25C_R,
                       #Ecophys_Jmax_25C_R,
                       best_Vcmax_25C_B[1:2],
                       
                       gstari_Tleaf.Asat,
                       Kmi.Asat,
                       TleafK.Asat,
                       Gma2[CRef][1],
                       Km2[CRef][1]),nrow=1)


  
  write.table (output, paste(arquivo, ".csv", sep=""), append=TRUE, sep=",", row.names=FALSE, col.names=F)
  #start_time <- Sys.time() 
  #end_time <- Sys.time()  
  #Tempo<- difftime(end_time, start_time, units='mins')*(length(curvas1)-which(curvas1 %in% i))
  # This is just a counter to visualize the progress of the analyzes
  print(paste(basename(sp[i,1])," - "))
  }

dev.off() 

