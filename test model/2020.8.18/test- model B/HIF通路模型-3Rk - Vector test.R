
# test-A HIFp-OH dehydroxylation simulation:
#adding a reverse V3R to simulate the dehydroxylation of HIF-pOH, 
#the reverse rate is varied from 1/10000; 1/1000; 1/100; 1/10

#condition of control£ºj=0 k=0 a=1 h=1 k_trans=0, uncomment the standard_vector

j <- 0.01 #k3R coefficient(10^?)
k <- 0 #VHL deficiency 0:-2 ~ 50:0.5
a <- 0 #VHL absent = 0
h <- 1 # O2 level
k_trans <- 0.0002 # k for VHL vector transcription (0.002)
p <- 40 # treat time (h)

#initial conditions

#Velocity constant
#Note: ** in the end indicates the Velocity constant has a preVious referential literature
  k1 <- 0.005
  k2 <- 0.0002
  k3 <- 0.045
  km3a <- 250000 #** reference:o2,Km3a&Km3b are in ¦ÌM,others are in nM
  km3b <- 100
  k4 <- 0.1
  km4 <- 150
  k9 <- 0.001
  k10 <- 0.001
  k11 <- 0.0001
  k12 <- 0.0001
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
  k15 = k3
  km15a = km3a
  km15b = km3b
  k16 = k4
  km16 = km4
  
  #set k for reverse reaction No.3
  k3R <- k3*j
  
#initial concentration 
  HIFa <- 0.05
  HIFa_pOH <- 0
  HIFan_pOH <- 0
  HIFan <- 0
  HIFd <- 0
  HIFd_HRE <- 0
  PHD <- 100
  PHDn <- 0
  HIFb <- 170
  HRE <- 50
  mRNA <- 0
  protein <- 0 # same as EGFP in the simplified model
  
  o2 <- 100000*h #** Oxygen leVel at normoxia
  
  #VHL initial concentration
  VHL <- 50*10^k*a
  VHLn <- 50*10^k*a
  
#set time interVal [unit="s"]
  dt <- 5  
  
#set loop to simulate the pathway
  #set timelength
  timesteps <- 100000
  #set empty Vectors to store the Variables
  s_HIFa <- numeric(length = timesteps)
  s_HIFa[1] <- HIFa
  s_HIFa_pOH <- numeric(length = timesteps)
  s_HIFa_pOH[1] <- HIFa_pOH
  s_HIFan_pOH <- numeric(length = timesteps)
  s_HIFan_pOH[1] <- HIFan_pOH
  s_HIFan <- numeric(length = timesteps)
  s_HIFan[1] <- HIFan
  s_HIFd <- numeric(length = timesteps)
  s_HIFd[1] <- HIFd
  s_HIFd_HRE <- numeric(length = timesteps)
  s_HIFd_HRE[1] <- HIFd_HRE
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
  s_VHL <- numeric(length = timesteps)
  s_VHL[1] <- VHL
  s_VHLn <- numeric(length = timesteps)
  s_VHLn[1] <- VHLn
  # the concentration of VHL is Viewed as a constant in this model
  
  
#set "for" loop
for (i in c(2:timesteps-1)){
 #equation for Velocity
  V1 <- k1
  V2 <- k2*HIFa
  V3 <- k3*PHD*(o2/(km3a+o2))*(HIFa/(km3b+HIFa))
  V4 <- k4*VHL*(HIFa_pOH/(km4+HIFa_pOH))
  V9 <- k9*HIFa
  V10 <- k10*HIFan
  V11 <- k11*PHD
  V12 <- k12*PHDn
  V15 <- k15*PHDn*(o2/(km15a+o2))*(HIFan/(km15b+HIFan))
  V16 <- k16*VHLn*(HIFan_pOH/(km16+HIFan_pOH))
  V21 <- k21f*HIFan*HIFb-k21r*HIFd
  V22 <- k22f*HIFd*HRE-k22r*HIFd_HRE
  V23 <- k23*HIFd_HRE
  V24 <- k24*HIFd_HRE
  V25 <- k25*PHD
  V26 <- k26*mRNA
  V27 <- k27*mRNA
  V28 <- k28*protein
  
  #set V3R to simulate dehydroxylation
  v3R <- k3R*HIFa_pOH
  #set v_trans to simulate VHL vector
  v_trans <- k_trans*HIFd_HRE
  
 #simulate the process  (v3R affects HIFa & HIFa_pOH)
  HIFa <- HIFa + (V1-V2-V9+V10-V3+v3R)*dt
  HIFa_pOH <- HIFa_pOH + (V3-V4-v3R)*dt
  HIFan <- HIFan + (V9-V10-V15-V21)*dt
  PHDn <- PHDn + (V11-V12)*dt
  HIFan_pOH <- HIFan_pOH + (V15-V16)*dt
  HIFd <- HIFd + (V21-V22)*dt
  HIFd_HRE <- HIFd_HRE + V22*dt
  PHD <- PHD + (V24-V25-V11+V12)*dt
  HIFb <- HIFb + (-V21)*dt
  HRE <- HRE + (-V22)*dt
  mRNA <- mRNA + (V23-V26)*dt
  protein <- protein + (V27-V28)*dt
  
  #VHL manipulation part
  VHLn <- VHLn
  if(i < p*3600/dt){
    VHL <- VHL
  }else{
    VHL <- VHL + v_trans*dt
  }
  
  s_HIFa[i+1] <- HIFa
  s_HIFa_pOH[i+1] <- HIFa_pOH
  s_HIFan_pOH[i+1] <- HIFan_pOH
  s_HIFan[i+1] <- HIFan
  s_HIFd[i+1] <- HIFd
  s_HIFd_HRE[i+1] <- HIFd_HRE
  s_PHD[i+1] <- PHD
  s_PHDn[i+1] <- PHDn
  s_HIFb[i+1] <- HIFb
  s_HRE[i+1] <- HRE
  s_mRNA[i+1] <- mRNA
  s_protein[i+1] <- protein
  s_VHL[i+1] <- VHL
  s_VHLn[i+1] <- VHLn
}

  #recording of balanced concentration level
    #standard_protein_balance <- s_protein[50000]
    #standard_PHD_balance <- s_PHD[50000]
    #standard_HIFa_balance <- s_HIFa[50000]
    #standard_HIFd_HRE_balance <- s_HIFd_HRE[50000]
  
#plot 
  time <- seq(0, by=dt/3600,length.out = timesteps)
  #plot of "HIF¦Á & HIF-dimer"
  plot(time,s_HIFd_HRE,col="red", type = "l",xlab = "time(h)",ylab = "concentration(nM)")
  abline(v = p, col="grey", lty=3)
  lines(time,s_HIFa,col ="blue")
  abline(h=standard_HIFd_HRE_balance, col="red",lty=2)
  abline(h=standard_HIFa_balance,col="blue",lty=2)
  legend("topright",legend = c("HIF¦Á","HIF-HRE"),col = c("blue","red"),lty = 1, cex = 1)
  
  #plot of "protein & PHD"
  plot(time,s_protein,col ="blue", type = "l", xlab = "time(h)", ylab = "concentration(nM)")
  abline(v = p, col="grey", lty=3)
  lines(time,s_PHD,col ="red",xlab = "time(h)",ylab = "concentration(nM)")
  abline(h = standard_PHD_balance,col = "red",lty = 2)
  abline(h = standard_protein_balance,col = "blue",lty = 2)
  legend("topright",legend = c("protein","PHD"),col = c("blue","red"),lty = 1, cex = 1)
  
  #plot of VHL change
  plot(time,s_VHL,col = "gold",type = "l",lwd = 2, xlab = "time(h)", ylab = "concentration(nM)")
  abline(v = p, col="grey", lty=3)
  legend("topleft",legend = "VHL", col = "gold", lty = 1, cex = 1)
  
#additional£ºgraph with ggplot
  #create a dataframe
  #HIF_df <- data.frame(time, s_HIFa,s_HIFa_pOH,
                      # s_HIFb,s_PHD,s_VHL,s_HIFd,
                      # s_HRE,HIFd_HRE,s_mRNA,s_protein)
  
  
  #ggplot
  #library(ggplot2)
  #plot of "HIF¦Á & HIF-dimer"
  #ggplot(HIF_df)+
   # geom_line(aes(x=time,y=s_HIFa,col="s_HIFa"))+
   # geom_line(aes(x=time,y=s_HIFd,col="s_HIFd"))+
   # scale_colour_discrete(breaks = c("s_HIFa","s_HIFd"), labels = c('HIF¦Á','HIF-dimer'))+
   # labs(x="time(h)",y="concentration(nM)")
  #plot of "protein & PHD"
  #ggplot(HIF_df)+
   # geom_line(aes(x=time,y=s_protein,col="s_protein"))+
   # geom_line(aes(x=time,y=s_PHD,col="s_PHD"))+
   # scale_colour_discrete(breaks = c("s_protein","s_PHD"), labels = c('protein','PHD'))+
   # labs(x="time(h)",y="concentration(nM)")
  
