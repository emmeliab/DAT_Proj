# Clean work area and close graphics
rm(list = ls()) # Remove all objects
graphics.off()  # Remove all graphics
cat("\014")     # Remove script in windows console

# Load relevant libraries
library(chron)
library(grDevices)
library(optimx)
library(plantecophys)
library(DEoptim)
library(weathermetrics)
library(minpack.lm)
library(readxl)

# Change the working directory and set names for output files
#setwd("/Users/maquellegarcia/Documents/GitHub/DAT_Proj/Results")
arquivo <-"A_ci_fit_DAT_Tapajos_20230206"
pdf(file=paste(arquivo, ".pdf", sep=""),height=10,width=20) # Create pdf to store graphs

# Create a matrix to store the column names of the output text file
output.names     <- matrix (c("Tree_id","DAT","person","Mean.Tleaf","Max.Tleaf",
                              "Asat", "Amax","Amin", "Ci.at.Asat","E.at.Asat","Ca","A/Ci","Ci/Ca",
                              "gs.at.Asat","gs.index","Mean.gs","mean.RH","Flow",
                              "mean.Press","Max.VPD","Ciindex",
                              "vcmax_Best_Model","Jmax_Best","Rd_Best",
                              "Best_Vcmax_25C","Best Jmax_25C",
                              "GammaStar","Kmi","TleafK.Asat"),nrow=1)
colnames (output.names) <- output.names
write.table (output.names, paste(arquivo, ".csv", sep=""), append=TRUE, sep=",",
             row.names=FALSE, col.names=FALSE)

# Function to fit Farquhar's biochemical model based on A-Ci curves
A_Ci_Fit1  <- function(x){
    Vcmax <- x[1]    # Variable 1: Vcmax
    Jmax  <- x[2]    # Variable 2: Jmax
    Rd    <- x[3]		 # Variable 3: Rd
    
    # Electron transport rate (J) when light is not saturating
    J    <- (Qabs + Jmax - sqrt((Qabs + Jmax)^2 - 4 * curv * Qabs * Jmax))/(2*curv) 
    aj   <- ifelse(Curve$Pci>=10, J * (Curve$Pci - gstari_Tleaf) / (4*Curve$Pci + 8 * gstari_Tleaf)-Rd,999)
    av   <- ifelse(Curve$Pci<50, Vcmax * (Curve$Pci - gstari_Tleaf) / (Kmi + Curve$Pci)-Rd, 999)# Carboxylation rate
    a     <- pmin(av,aj) # Modeled assimilation rate
    adiff <- sum((Curve$A - a)^2) # Sum of squared differences between modeled and observed assimilation
  #diffsum <- adiff +  abs(28-Cseq[max(head(sort((av-aj)^2, index.return=TRUE)$ix,2))]) + (Rd*(-1))#essa parte adiciona a penalizacao
}
# Constants used in the Farquhar's model, not considering mesophyll conductance
R <- 0.008314 # Gas constant
R2 <- R * 1000
Kc <- 40.49 # Michaelis-Menten constant for CO2 (Pa) (Bernacchi et al 2001,2002) or 404.9 microbar von Caemmerer et al. (1994)
delta_Kc <- 79430 # (J mol-1) from Medlyn et al 2002
Ko <- 27.84 # Michaelis-Menten constant for CO2 (kPa)(Bernacchi et al 2001,2002) or 278.4 mbar von Caemmerer et al. (1994)
delta_Ko <- 36380 # (J mol-1) from Medlyn et al 2002
gstar <- 4.275 # CO2 compensation point (Pa) (Bernacchi et al 2001,2002) or 36.9 microbar von Caemmerer et al. (1994)
delta_gstar <- 37830 # (J mol-1)

# Other Constants
curv <- 0.85 # Curvature for calculating J from Evans (1989)
Ambient.CO2 <- 400
O2 <- 21 # Estimated Oxygen concentration at chloroplast - (kPa)
K0 <- 0.388 # for accounting for diffusion of CO2 through the gasket (Licor 6400 only), maybe this part can be removed, as the new machine is already accounting for it

# Function to convert parameters to standard temperature 25°C according to Bernacchi
Temp_Correction <- function (p1,p2,p3) {
    a <- p1/(exp(26.355 - (65.33 /(R*meanTleafK)))) # Vcmax as in Bernacchi
    b <- p2/(exp( 17.71 - (43.9  /(R*meanTleafK)))) # this is different from Bernacchi
    # c <- p3/((exp(21.46 - 53.1   /(R*meanTleafK)))/(1+exp((0.65*meanTleafK-201.8)/(R*meanTleafK))))
    d <- p3/(exp(18.715 - (46.39 /(R*meanTleafK)))) # Rd as in Bernacchi
    
    return(c(a,b,d))
}

### Peaked function Medlin et al. 2002PCE
### Ea = Activation Energy, DeltaS = entropy factor = 200, Hd = Deactivation Energy)
### updated to Kumarathunge et al 2019, New Phytologist.(222: 768-784) doi: 10.1111/nph.15668
Ea_V     <- 82620.87 #58550 or 82620.87 or ---- these values are slight different from Kamarathunge
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
modeled_lines <- function (p1,p2,p3){
  av   <- p1 * (Cseq - mean(gstari_Tleaf))/(mean(Kmi) + Cseq ) - p3 
  J    <- mean((Qabs + p2 - sqrt((Qabs + p2)^2 - 4 * curv * Qabs * p2))/(2*curv) )
  aj   <- J * (Cseq - mean(gstari_Tleaf))/(4*Cseq + 8*mean(gstari_Tleaf)) - p3 
  #atpu <- p3 * 3 - p4 + (Cseq /10^10) 
  amin <- pmin(av,aj)
  return(cbind(av,aj,amin))
}
modeled_points <- function (p1,p2,p3){
  avb   <- p1 * (Curve$Pci - mean(gstari_Tleaf))/(mean(Kmi) + Curve$Pci ) - p3
  J    <- mean((Qabs + p2 - sqrt((Qabs + p2)^2 - 4 * curv * Qabs * p2))/(2*curv) )
  ajb   <- J * (Curve$Pci - mean(gstari_Tleaf)) / (4*Curve$Pci + 8 * mean(gstari_Tleaf)) - p3
  #atpub <- p3 * 3 - p4 + (Curve$Pci /10^10)
  aminb <- pmin(avb,ajb)
  return(cbind(avb,ajb,aminb))
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
  
 
  # converts Ci from ppm to Pa
  Curve$Pci <- Curve$Ci * Curve$Pa * 0.001
  
  # Convert Tleaf to kelvin
  Tleafk <- Curve$Tleaf + 273.15
  
  # Calculates a standardized measure of stomatal conductance change during measurements
  gsindex <- (max(Curve$gsw) - min(Curve$gsw)) / max(Curve$gsw)
  Ciindex <- (max(Curve$Ci) - min(Curve$Ci)) / max(Curve$Ci)
  
  # Enzyme kinetics at actual leaf temperature
  Kci_Tleaf <- Kc * exp((Tleafk - 298.15) * delta_Kc / (298.15 * 8.314 * Tleafk))
  Koi_Tleaf <- Ko * exp((Tleafk - 298.15) * delta_Ko / (298.15 * 8.314 * Tleafk))
  gstari_Tleaf <- gstar * exp((Tleafk - 298.15) * delta_gstar / (298.15 * 8.314 * Tleafk))
  Kmi <- Kci_Tleaf * (1 + O2 / Koi_Tleaf)
  
  # Light absorbed
  Qabs <- Curve$Qin * 0.85 * (1 - 0.15) / 2
  Curve$Qabs <- Curve$Qin * 0.85 * (1 - 0.15) / 2
  
  # Variables to be exported to final results file
  Curve$Ci_Ca <- Curve$Ci / Curve$CO2_s
  Curve$WUE <- Curve$A / Curve$gsw
  Curve$CUE <- Curve$A / Curve$Ci
  Curve$VUE <- Curve$gsw / Curve$Ci
  
  # Finds the best Ci/Ca ratio (based on an expected ratio of 0.7)
  find.best <- try(nlsLM(VUE ~ a * Ci^-b, control = nls.control(maxiter = 100, warnOnly = T), data = Curve, start = c(a = .2, b = 1.5)))
  
  # If the optimization does not work, select VUE at 350-450 ppm CO2
  if (is(find.best, "try-error")) {
      CRef <- if (isTRUE(length(Curve$VUE[Curve$CO2_r < 450 & Curve$CO2_r > 350]) >= 1)) {
          which(Curve$VUE == max(Curve$VUE[Curve$CO2_r < 450 & Curve$CO2_r > 350]))
      } else {
          NA
      }
  } else {
      # If the optimization works, select the minimum residual value at 350-450 ppm CO2
      CRef <- if (isTRUE(length(Curve$VUE[Curve$CO2_r < 450 & Curve$CO2_r > 350]) >= 1)) {
          which(abs(resid(find.best)) == min(abs(resid(find.best)[Curve$CO2_r < 450 & Curve$CO2_r > 350])))
      } else {
          NA
      }
  }
  
  # Variables to be exported to final results file
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

  # Starting values for the "optim" function, based on relationships between Asat or Amax  and Vcmax, Jmax, TPU, and Rd
  # these are just coefficients from fitted linear regressions derived from previous fittings of TROBIT curves
  MeanA <- mean(Curve$A)
  start.values1 <- c(MeanA * 7.4318+16.927, MeanA * 6.0787+58.343, 0.8)
  

 # parameters (Vcmax, Jmax, TPU and Rd) minimum and maximum values for the "DEoptim" function
 lower.bound  <- c(-15, -15,-15)     
 upper.bound  <- c(500, 500, 15)
 

   
# Fitting the A-ci model  using the "optimx" minimization algorithms
  first    <-  subset(optimx(start.values1, A_Ci_Fit1,control=list(all.methods=TRUE)), convcode < 1)

# organizing results from the curve fitting
  fits1     <- summary(if(sum( first[3]<0) < length( first$p1)) subset( first, p3 > 0 ) else  first, order= c(value))[1,c(1:7)]
  fits1_algorithim  <- row.names(fits1)
  fits1 <- unlist(fits1)
#################### calculating rates from the curve fitting results in order to build figures
  ## CO2 sequence (will be needed for building the figures)
  Cseq <- seq(0,200,0.5)
  
  modeled_1 <- do.call(modeled_points,as.list(fits1[1:3]))
  fits1_erro <- sqrt(sum((Curve$A - modeled_1[,3])^2))/length(Curve$A)
  result_lines_1 <- do.call(modeled_lines,as.list(fits1[1:3]))
 
  result_lines<-data.frame(result_lines_1)
  
###################################
  Trans_fit1 <- max(head(sort((result_lines_1[,1]-result_lines_1[,2])^2, index.return=TRUE)$ix,2))
  best        <- fits1[1:7]
  
###########################
  Tleaf_k <- Curve$Tleaf +273.15
  meanTleafK    <- mean(Tleafk)
 

##### temperature conversions to a reference temperature    

  best_Vcmax_25C        <- Vc_peaked_25C(       best[[1]],Ea_V, Delta_V, Ed)
  best_Jmax_25C         <- J_peaked_25C(       best[[2]], Ea_J, Delta_J, Ed)

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
  title(cex.main = 1.5,main = c (as.character(unique(Curve[,"unique_id"]))))
  lines (Cseq, result_lines$amin,   col = "sienna2", 	lwd = "4")
  lines (Cseq, result_lines$av,   col = "red",   lwd = "4")
 # lines (Cseq, result_lines$atpu, col = "plum", lwd = "4")
  lines (c(27,27),c(-5,45), col = "black", lty = 3)
  lines (c(45,45),c(-5,45), col = "black", lty = 3)
  text (190,f5,paste("Vcmax =",round(best[1],digits=1)),pos=2,col = "red",cex=1.5)
  text (190,f4,paste("Jmax =" ,round(best[2],digits=1)),pos=2,col = "sienna2",cex=1.5)
  text (190,f2,paste("Rd ="   ,round(best[3],digits=2)),pos=2,cex=1.5)
  text ( 50,f3,paste("error =",round(best[4],digits=3)),pos=4,col = "red", cex = 1.5)
  text ( 50,f2,paste("gs ="   ,round(mean(Curve$gsw),digits=2)),pos=4,cex=1.5)
  text ( 50,f1,paste("Tleaf =",round(mean(Curve$Tleaf),digits=1)),pos=4,cex=1.5)
  text ( 50,f4,paste("Ci/Ca =",round(Ci.Ca,digits=2)),pos=4,cex=1.5)
  #points(Cseq[df2$Transition[indice[1]]],result_lines$av[df2$Transition[indice[1]]],col ="green",pch=9, cex = 2)
  
  
  names(Curve)
  
  

  #####################  export analyses to the results table
  output <- matrix (c (unique(Curve[,"unique_id"]),unique(Curve[,"Data_point"]),unique(Curve[,"Scientist"]),mean(Curve$Tleaf),max(Curve$Tleaf),
                       ASat, Amax,min(Curve$A),
                       
                       Ci.Asat,E.Asat,Ca.Asat,(ASat/Ci.Asat),Ci.Ca,gs.Asat,
                       
                       gsindex,mean(Curve$gsw), 
                       
                       mean(Curve$RHcham),mean(Curve$Flow),mean(Curve$Pa),
                       max(Curve$VPDleaf),Ciindex,
                       
                     
                       
                       best[1:3],
                  
                       
                       #Ecophys_Vcmax_25C,
                       #Ecophys_Jmax_25C,
                       #Ecophys_Vcmax_25C_R,
                       #Ecophys_Jmax_25C_R,
                       best_Vcmax_25C,
                       best_Jmax_25C,
                       
                       gstari_Tleaf.Asat,
                       Kmi.Asat,
                       TleafK.Asat),nrow=1)


  
  write.table (output, paste(arquivo, ".csv", sep=""), append=TRUE, sep=",", row.names=FALSE, col.names=F)
  #start_time <- Sys.time() 
  #end_time <- Sys.time()  
  #Tempo<- difftime(end_time, start_time, units='mins')*(length(curvas1)-which(curvas1 %in% i))
  # This is just a counter to visualize the progress of the analyzes
  print(paste(basename(sp[i,1])," - "))
  }

dev.off() 

