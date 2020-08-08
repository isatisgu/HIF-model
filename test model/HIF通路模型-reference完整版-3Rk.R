#modified purpose:
#adding a reverse V3R to simulate the dehydroxylation of HIF-pOH

#initial conditions

#Velocity constant
#Note: ** in the end indicates the Velocity constant has a preVious referential literature
  k1 <- 0.005
  k2 <- 0.0002
  k3 <- 0.045
  km3a <- 250000 #** reference:o2,Km3a&Km5a are in ¦ÌM,others are in nM
  km3b <- 100
  k4 <- 0.1
  km4 <- 150
  k5 <- 0.4
  km5a <- 90000 #**
  km5b <- 20
  k6 <- 0.001
  k7 <- 0.003
  k8 <- 0.01
  k9 <- 0.001
  k10 <- 0.001
  k11 <- 0.0001
  k12 <- 0.0001
  k13 <- 0.0002
  k14 <- 0.0001
  k21f <- 0.0005
  k21r <- 0.01
  k22f <- 0.001
  k22r <- 0.01
  k23 <- 0.0016 #**
  k24 <- 0.002
  k25 <- 0.00001 #based on their >20 hrs measured half-life
  k26 <- 0.00038 #based on their ~30 mins measured half-life
  k27 <- 0.009 #**
  k28 <- 0.0001 #based on their ~2 hrs measured half-life
  k29 <- 0.02
  km7a = km3a
  km7b = km3b
  km8 = km4
  k15 = k3
  km15a = km3a
  km15b = km3b
  k16 = k4
  km16 = km4
  k17 = k5
  km17a = km5a
  km17b = km5b
  k18 = k6
  k19 = k7
  km19a = km15a
  km19b = km15b
  k20 = k8
  km20 = km16
  
  #set k for reverse reaction No.3
  k3R <- 0.01
  k7R <- 0.001
  
#initial concentration 
  HIFa <- 0.05
  HIFa_pOH <- 0
  HIFa_aOH <- 0
  HIFa_aOHpOH <- 0
  HIFan_pOH <- 0
  HIFan <- 0
  HIFd <- 0
  HIFd_HRE <- 0
  HIFan_aOH <- 0
  HIFan_aOHpOH <- 0
  PHD <- 100
  PHDn <- 0
  HIFb <- 170
  HRE <- 50
  mRNA <- 0
  protein <- 0 # same as EGFP in the simplified model
  luciferase <- 0
  
  o2 <- 100000 #** Oxygen leVel at normoxia
  FIH <- 110
  FIHn <- 40
  
  #simulate VHL deficiency
  VHL <- 50*0.1
  VHLn <- 50*0.1
  
  
#set time interVal [unit="s"]
  dt <- 0.5  
  
#set loop to simulate the pathway
  #set timelength
  timesteps <- 2000000
  #set empty Vectors to store the Variables
  s_HIFa <- numeric(length = timesteps)
  s_HIFa[1] <- HIFa
  s_HIFa_pOH <- numeric(length = timesteps)
  s_HIFa_pOH[1] <- HIFa_pOH
  s_HIFa_aOH <- numeric(length = timesteps)
  s_HIFa_aOH[1] <- HIFa_aOH
  s_HIFa_aOHpOH <- numeric(length = timesteps)
  s_HIFa_aOHpOH[1] <- HIFa_aOHpOH
  s_HIFan_pOH <- numeric(length = timesteps)
  s_HIFan_pOH[1] <- HIFan_pOH
  s_HIFan <- numeric(length = timesteps)
  s_HIFan[1] <- HIFan
  s_HIFd <- numeric(length = timesteps)
  s_HIFd[1] <- HIFd
  s_HIFd_HRE <- numeric(length = timesteps)
  s_HIFd_HRE[1] <- HIFd_HRE
  s_HIFan_aOH <- numeric(length = timesteps)
  s_HIFan_aOH[1] <- HIFan_aOH
  s_HIFan_aOHpOH <- numeric(length = timesteps)
  s_HIFan_aOHpOH[1] <- HIFan_aOHpOH
  s_PHD <- numeric(length = timesteps)
  s_PHD[1] <- PHD
  s_PHDn <- numeric(length = timesteps)
  s_PHDn[1] <- PHDn
  s_HIFb <- numeric(length = timesteps)
  s_HIFb[1] <- HIFb
  s_HRE <- numeric(length = timesteps)
  s_HRE[1] <- HRE
  s_mRNA <- numeric(length = timesteps)
  s_mRNA[1] <- mRNA
  s_protein <- numeric(length = timesteps)
  s_protein[1] <- protein
  s_luciferase <- numeric(length = timesteps)
  s_luciferase[1] <- luciferase
  s_VHL <- numeric(length = timesteps)
  s_VHL[1] <- VHL
  # the concentration of VHL is Viewed as a constant in this model
  
  
#set "for" loop
for (i in c(2:timesteps-1)){
 #equation for Velocity
  V1 <- k1
  V2 <- k2*HIFa
  V3 <- k3*PHD*(o2/(km3a+o2))*(HIFa/(km3b+HIFa))
  V4 <- k4*VHL*(HIFa_pOH/(km4+HIFa_pOH))
  V5 <- k5*FIH*(o2/(km5a+o2))*(HIFa/(km5b+HIFa))
  V6 <- k6*HIFa_aOH
  V7 <- k7*PHD*(o2/(km7a+o2))*(HIFa_aOH/(km7b+HIFa_aOH))
  V8 <- k8*VHL*(HIFa_aOHpOH/(km8+HIFa_aOHpOH))
  V9 <- k9*HIFa
  V10 <- k10*HIFan
  V11 <- k11*PHD
  V12 <- k12*PHDn
  V13 <- k13*HIFa_aOH
  V14 <- k14*HIFan_aOH
  V15 <- k15*PHDn*(o2/(km15a+o2))*(HIFan/(km15b+HIFan))
  V16 <- k16*VHLn*(HIFan_pOH/(km16+HIFan_pOH))
  V17 <- k17*FIHn*(o2/(km17a+o2))*(HIFan/(km17b+HIFan))
  V18 <- k18*HIFan_aOH
  V19 <- k19*PHDn*(o2/(km19a+o2))*(HIFan_aOH/(km19b+HIFan_aOH))
  V20 <- k20*VHLn*(HIFan_aOHpOH/(km20+HIFan_aOHpOH))
  V21 <- k21f*HIFan*HIFb-k21r*HIFd
  V22 <- k22f*HIFd*HRE-k22r*HIFd_HRE
  V23 <- k23*HIFd_HRE
  V24 <- k24*HIFd_HRE
  V25 <- k25*PHD
  V26 <- k26*mRNA
  V27 <- k27*mRNA
  V28 <- k28*protein
  V29 <- k29*protein
  
  #set reverse velocity to simulate dehydroxylation
  v3R <- k3R*HIFa_pOH
  v7R <- k7R*HIFa_aOHpOH
  
 #simulate the process
  HIFa <- HIFa + (V1-V2-V9+V10-V3-V5+V6+v3R)*dt
  HIFa_pOH <- HIFa_pOH + (V3-V4-v3R)*dt
  HIFa_aOH <- HIFa_aOH + (V5-V6-V7-V13+V14+v7R)*dt
  HIFa_aOHpOH <- HIFa_aOHpOH + (V7-V8-v7R)*dt
  HIFan_pOH <- HIFan_pOH + (V15-V16)*dt
  HIFan <- HIFan + (V9-V10-V17+V18-V15-V21)*dt
  HIFd <- HIFd + (V21-V22)*dt
  HIFd_HRE <- HIFd_HRE + V22*dt
  HIFan_aOH <- HIFan_aOH + (V17-V18-V19)*dt
  HIFan_aOHpOH <- HIFan_aOHpOH + (V19-V20)*dt
  PHD <- PHD + (V24-V25-V11+V12)*dt
  PHDn <- PHDn + (V11-V12)*dt
  HIFb <- HIFb + (-V21)*dt
  HRE <- HRE + (-V22)*dt
  mRNA <- mRNA + (V23-V26)*dt
  protein <- protein + (V27-V28)*dt
  luciferase <- luciferase + V29*dt
  VHL <- VHL
  s_HIFa[i+1] <- HIFa
  s_HIFa_pOH[i+1] <- HIFa_pOH
  s_HIFa_aOH[i+1] <- HIFa_aOH
  s_HIFa_aOHpOH[i+1] <- HIFa_aOHpOH
  s_HIFan_pOH[i+1] <- HIFan_pOH
  s_HIFan[i+1] <- HIFan
  s_HIFd[i+1] <- HIFd
  s_HIFd_HRE[i+1] <- HIFd_HRE
  s_HIFan_aOH[i+1] <- HIFan_aOH
  s_HIFan_aOHpOH[i+1] <- HIFan_aOHpOH
  s_PHD[i+1] <- PHD
  s_PHDn[i+1] <- PHDn
  s_HIFb[i+1] <- HIFb
  s_HRE[i+1] <- HRE
  s_mRNA[i+1] <- mRNA
  s_protein[i+1] <- protein
  s_luciferase[i+1] <- luciferase
  s_VHL[i+1] <- VHL
}

#plot 
  time <- seq(0, by=dt/3600,length.out = timesteps)
  #plot of "HIF1¦Á & HIF-dimer"
  plot(time,s_HIFa,col ="blue", type = "l",ylim = c(0,0.05),xlab = "time(h)",ylab = "concentration(nM)")
  lines(time,s_HIFd_HRE,col="red")
  legend("topleft",legend = c("HIF¦Á","HIF-HRE"),col = c("blue","red"),lty = 1, cex = 1)
  #plot of "protein & PHD"
  plot(time,s_PHD,col ="red", type = "l", ylim = c(0,50),xlab = "time(h)",ylab = "concentration(nM)")
  lines(time,s_protein,col ="blue",xlab = "time(h)",ylab = "concentration(nM)")
  legend("topleft",legend = c("protein","PHD"),col = c("blue","red"),lty = 1, cex = 1)
  
  plot(time,s_PHD,col ="red", type = "l", xlab = "time(h)",ylab = "concentration(nM)")
  plot(time,s_protein,col ="blue",type = "l",xlab = "time(h)",ylab = "concentration(nM)")
  legend("topleft",legend = c("protein","PHD"),col = c("blue","red"),lty = 1, cex = 1)
   
#additional£ºgraph with ggplot
  #create a dataframe
  HIF_df <- data.frame(time, s_HIFa,s_HIFa_pOH,
                       s_HIFb,s_PHD,s_VHL,s_HIFd,
                       s_HRE,HIFd_HRE,s_mRNA,s_protein)
  
  
  
  #ggplot
  library(ggplot2)
  #plot of "HIF1¦Á & HIF-dimer"
  ggplot(HIF_df)+
    geom_line(aes(x=time,y=s_HIFa,col="s_HIFa"))+
    geom_line(aes(x=time,y=s_HIFd,col="s_HIFd"))+
    scale_colour_discrete(breaks = c("s_HIFa","s_HIFd"), labels = c('HIF1¦Á','HIF-dimer'))+
    labs(x="time(h)",y="concentration(nM)")
  #plot of "protein & PHD"
  ggplot(HIF_df)+
    geom_line(aes(x=time,y=s_protein,col="s_protein"))+
    geom_line(aes(x=time,y=s_PHD,col="s_PHD"))+
    scale_colour_discrete(breaks = c("s_protein","s_PHD"), labels = c('protein','PHD'))+
    labs(x="time(h)",y="concentration(nM)")
  
