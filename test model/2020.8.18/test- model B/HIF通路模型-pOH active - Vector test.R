# test-B HIF_p-OH remains active simulation:
#add the reaction of 9p 10p 21p 22p 23p

j <- 0.5 # 0:1 ~ 0%:100% efficiency
k <- 0 #VHL deficiency
a <- 0 #VHL absent = 0
h <- 1 # 1% O2 = 1/21; normoxia = 1
k_trans <- 0.002 # k for VHL vector transcription (0.002)
p <- 40 # treat time (h)

#initial conditions

#velocity constant
#Note: ** in the end indicates the velocity constant has a previous referential literature
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
  
  
#initial concentration 
  HIFa <- 0.05
  HIFa_pOH <- 0
  HIFan_pOH <- 0
  HIFan <- 0
  HIFd <- 0
  HIFd_HRE <- 0
  PHD <- 100
  PHDn <- 0
  HIFb <- 170 #assume to be constant
  HRE <- 50 #assume to be constant
  mRNA <- 0
  protein <- 0 # same as EGFP in the simplified model
  
  o2 <- 100000*h #** Oxygen level at normoxia
  
  #VHL initial concentration
  VHL <- 50*10^k*a
  VHLn <- 50*10^k*a
  
  #SET pOH substance
  HIFd_pOH <- 0
  HIFd_pOH_HRE <- 0
  
#set time interval [unit="s"]
  dt <- 5  
  
#set loop to simulate the pathway
  #set timelength
  timesteps <- 100000
  #set empty vectors to store the variables
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
  # the concentration of VHL is viewed as a constant in this model
  
  #set pOH vectors
  s_HIFd_pOH <- numeric(length = timesteps)
  s_HIFd_pOH[1] <- HIFd_pOH
  s_HIFd_pOH_HRE <- numeric(length = timesteps)
  s_HIFd_pOH_HRE[1] <- HIFd_pOH_HRE
  
  s_HIFa_active <- numeric(length = timesteps)
  s_HIFa_active[1] <- HIFa + HIFa_pOH
  s_HIFd_HRE_active <- numeric(length = timesteps)
  s_HIFd_HRE_active[1] <- HIFd_HRE + HIFd_pOH_HRE
  
#set "for" loop
for (i in c(2:timesteps-1)){
 #equation for velocity
  v1 <- k1
  v2 <- k2*HIFa
  v3 <- k3*PHD*(o2/(km3a+o2))*(HIFa/(km3b+HIFa))
  v4 <- k4*VHL*(HIFa_pOH/(km4+HIFa_pOH))
  v9 <- k9*HIFa
  v10 <- k10*HIFan
  v11 <- k11*PHD
  v12 <- k12*PHDn
  v15 <- k15*PHDn*(o2/(km15a+o2))*(HIFan/(km15b+HIFan))
  v16 <- k16*VHLn*(HIFan_pOH/(km16+HIFan_pOH))
  v21 <- k21f*HIFan*HIFb-k21r*HIFd
  v22 <- k22f*HIFd*HRE-k22r*HIFd_HRE
  v23 <- k23*HIFd_HRE
  v24 <- k24*HIFd_HRE
  v25 <- k25*PHD
  v26 <- k26*mRNA
  v27 <- k27*mRNA
  v28 <- k28*protein
  
  #set v9p v10p v21p to simulate HIFan-pOH remains active
  v9P <- k9*HIFa_pOH
  v10P <- k10*HIFan_pOH
  v21P <- k21f*HIFan_pOH*HIFb-k21r*HIFd_pOH
  v22P <- k22f*HIFd_pOH*HRE-k22r*HIFd_pOH_HRE
  v23P <- k23*HIFd_pOH_HRE*j
  #set v_trans to simulate VHL vector
  v_trans <- k_trans*HIFd_HRE #only HIFd_HRE is counted, causing the causation of delayed boost
  
 #simulate the process  (v3R affects HIFa & HIFa_pOH)
  HIFa <- HIFa + (v1-v2-v9+v10-v3)*dt
  HIFa_pOH <- HIFa_pOH + (v3-v4-v9P+v10P)*dt #add V9P,V10P
  HIFan <- HIFan + (v9-v10-v15-v21)*dt
  PHDn <- PHDn + (v11-v12)*dt
  HIFan_pOH <- HIFan_pOH + (v15-v16+v9P-v10P-v21P)*dt #add V9P,V10P,V21p
  HIFd <- HIFd + (v21-v22)*dt
  HIFd_pOH <- HIFd_pOH +(v21P-v22P)*dt #new variable (1)
  HIFd_HRE <- HIFd_HRE + v22*dt
  HIFd_pOH_HRE <- HIFd_pOH_HRE + v22P*dt #new variable(2)
  PHD <- PHD + (v24-v25-v11+v12)*dt #Assume HIFd_pOH_HRE does not activate PHD transcription
  HIFb <- HIFb + (-v21-v21P)*dt #add v21P [maximum of protein is confined by HIFb and HRE]
  HRE <- HRE + (-v22-v22P)*dt #add v22P  
  mRNA <- mRNA + (v23-v26+v23P)*dt # add v23P
  protein <- protein + (v27-v28)*dt
  
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
  
  #new vectors
  s_HIFd_pOH[i+1] <- HIFd_pOH
  s_HIFd_pOH_HRE[i+1] <- HIFd_pOH_HRE
  
  s_HIFa_active[i+1] <- HIFa + HIFa_pOH
  s_HIFd_HRE_active[i+1] <- HIFd_HRE + HIFd_pOH_HRE
  
}
  
  #recording of balanced concentration level
    #standard_protein_balance <- s_protein[250000]
    #standard_PHD_balance <- s_PHD[250000]
    
    #standard_HIFd_HRE_active_balance <- s_HIFd_HRE_active[250000]
    #standard_HIFd_HRE_balance <- s_HIFd_HRE[250000]
    #standard_HIFa_balance <- s_HIFa[250000]
    #standard_HIFa_active_balance <- s_HIFa_active[250000]
  
#plot 
  time <- seq(0, by=dt/3600,length.out = timesteps)
  #plot of "HIF1¦Á & HIF-dimer"
  plot(time,s_HIFd_HRE_active,col="red", type = "l",xlab = "time(h)",ylab = "concentration(nM)")
  lines(time,s_HIFa_active,col ="green")
  #plot(time,s_HIFa_active,col ="green", type = "l",xlab = "time(h)",ylab = "concentration(nM)")
  #lines(time,s_HIFd_HRE_active,col="red")
  abline(v = p, col="grey", lty=3)
  lines(time,s_HIFa,col ="blue")
  lines(time,s_HIFd_HRE,col = "grey")
  abline(h=standard_HIFd_HRE_active_balance, col="red",lty=2)
  abline(h=standard_HIFa_balance,col="blue",lty=2)
  abline(h=standard_HIFd_HRE_balance, col = "grey", lty=2)
  abline(h=standard_HIFa_active_balance, col="green", lty=2)
  legend("topright",legend = c("HIF¦Á","HIF¦Á + HIF¦Á_pOH","HIF_HRE","HIF_HRE + HIF_pOH_HRE"),col = c("blue","green","grey","red"),lty = 1, cex = 1)
  
  #plot of "protein & PHD"
  plot(time,s_protein,col ="blue", type = "l",xlab = "time(h)",ylab = "concentration(nM)")
  lines(time,s_PHD,col ="red",xlab = "time(h)",ylab = "concentration(nM)")
  abline(v = p, col="grey", lty=3)
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
  #plot of "HIF1¦Á & HIF-dimer"
  #ggplot(HIF_df)+
   # geom_line(aes(x=time,y=s_HIFa,col="s_HIFa"))+
   # geom_line(aes(x=time,y=s_HIFd,col="s_HIFd"))+
   # scale_colour_discrete(breaks = c("s_HIFa","s_HIFd"), labels = c('HIF1¦Á','HIF-dimer'))+
   # labs(x="time(h)",y="concentration(nM)")
  #plot of "protein & PHD"
  #ggplot(HIF_df)+
   # geom_line(aes(x=time,y=s_protein,col="s_protein"))+
   # geom_line(aes(x=time,y=s_PHD,col="s_PHD"))+
   # scale_colour_discrete(breaks = c("s_protein","s_PHD"), labels = c('protein','PHD'))+
   # labs(x="time(h)",y="concentration(nM)")
  
