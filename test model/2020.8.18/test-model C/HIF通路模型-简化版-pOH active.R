# test-B HIF_p-OH remains active simulation:
#add the reaction of 21p 22p 23p

j <- 0 # 0:1 ~ 0%:100% efficiency
k <- 0 #vHL deficiency

#initial conditions

#velocity constant
  k1 <- 0.005
  k2 <- 0.0002
  k3 <- 0.045
  km3a <- 250000 
  km3b <- 100
  # reference:O2,km3a&km5a are in ¦ÌM,others are in nM
  k4 <- 0.1
  km4 <- 150
  k5 <- 0.001*0.0005 #k9*k21f
  k6f <- 0.001 #k22
  k6r <- 0.01
  k7 <- 0.00001 #k25
  # v8&v9 is the velocity of pVHL catalysis and degration, which is not shown in the reference.
  k10 <- 0.009 #k27
  k11 <- 0.0001 #k28
  k12 <- 0.002 #k24 
  k13 <- 0.001*0.01 #k10*k21r
  k14 <- 0.0016 #k23
  k15 <- 0.00038 #k26
  
#initial concentration 
  HIFa <- 0.05
  O2 <- 100000
  HIFa_pOH <- 0
  HIFb <- 170
  PHD <- 100
  HIFd <- 0
  HRE <- 50
  HIFd_HRE <- 0
  mRNA <- 0
  protein <- 0 # same as protein in reference
  
  #SET pOH substance
  HIFd_pOH <- 0
  HIFd_pOH_HRE <- 0
  
  #VHL deficiency
  VHL <- 50*10^k
  
#set time interval [unit="s"]
  dt <- 10
  
#set loop to simulate the pathway
  #set timelength
  timesteps <- 50000
  #set empty vectors to store the variables
  s_HIFa <- numeric(length = timesteps) #equilibrium around 0.4
  s_HIFa[1] <- HIFa
  s_HIFa_pOH <- numeric(length = timesteps)
  s_HIFa_pOH[1] <- HIFa_pOH
  s_HIFb <- numeric(length = timesteps)
  s_HIFb[1] <- HIFb
  s_PHD <- numeric(length = timesteps)
  s_PHD[1] <- PHD
  s_VHL <- numeric(length = timesteps)
  s_VHL[1] <- VHL
  s_HIFd <- numeric(length = timesteps)
  s_HIFd[1] <- HIFd
  s_HRE <- numeric(length = timesteps)
  s_HRE[1] <- HRE
  s_HIFd_HRE <- numeric(length = timesteps)
  s_HIFd_HRE[1] <- HIFd_HRE
  s_mRNA <- numeric(length = timesteps)
  s_mRNA[1] <- mRNA
  s_protein <- numeric(length = timesteps)
  s_protein[1] <- protein 
  
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
  v3 <- k3*PHD*(O2/(km3a+O2))*(HIFa/(km3b+HIFa))
  v4 <- k4*VHL*(HIFa_pOH/(km4+HIFa_pOH))
  v5 <- k5*HIFa*HIFb
  v6 <- k6f*HIFd*HRE-k6r*HIFd_HRE
  v7 <- k7*PHD
  #v8&v9 currently empty
  v10 <- k10*mRNA
  v11 <- k11*protein
  v12 <- k12*HIFd_HRE
  v13 <- k13*HIFd
  v14 <- k14*HIFd_HRE
  v15 <- k15*mRNA
  
  #set v5p v13p v6p v14p to simulate HIFan-pOH remains active 
  v5p <- k5*HIFa_pOH*HIFb
  v13p <- k13*HIFd_pOH
  v6p <- k6f*HIFd_pOH*HRE-k6r*HIFd_pOH_HRE
  v14p <- k14*HIFd_pOH_HRE*j
  
 #simulate the process
  HIFa <- HIFa + (v1-v2-v3-v5+v13)*dt
  HIFa_pOH <- HIFa_pOH + (v3-v4-v5p+v13p)*dt #add v5p v13p
  HIFd_pOH <-HIFd_pOH + (v5p-v13p-v6p)*dt #new variable (1)
  HIFb <- HIFb + (v13-v5+v13p-v5p)*dt #add v5p v13p [maximum of protein is confined by HIFb and HRE]
  HRE <- HRE + (-v6-v6p)*dt #add v6p
  HIFd_pOH_HRE <- HIFd_pOH_HRE + v6p*dt #new variable(2)
  PHD <- PHD + (v12-v7)*dt #Assume HIFd_pOH_HRE does not activate PHD transcription
  VHL <- VHL
  HIFd <- HIFd + (v5-v13-v6)*dt
  HIFd_HRE <- HIFd_HRE + v6*dt
  mRNA <- mRNA + (v14-v15+v14p)*dt #add v14p
  protein <- protein + (v10-v11)
  s_HIFa[i+1] <- HIFa
  s_HIFa_pOH[i+1] <- HIFa_pOH
  s_HIFb[i+1] <- HIFb
  s_PHD[i+1] <- PHD
  s_VHL[i+1] <- VHL
  s_HIFd[i+1] <- HIFd
  s_HRE[i+1] <- HRE
  s_HIFd_HRE[i+1] <- HIFd_HRE
  s_mRNA[i+1] <- mRNA
  s_protein[i+1] <- protein
  
  #new vectors
  s_HIFd_pOH[i+1] <- HIFd_pOH
  s_HIFd_pOH_HRE[i+1] <- HIFd_pOH_HRE
  
  s_HIFa_active[i+1] <- HIFa + HIFa_pOH
  s_HIFd_HRE_active[i+1] <- HIFd_HRE + HIFd_pOH_HRE
}
  
 #recording of balanced concentration level
  #standard_protein_balance <- s_protein[50000]
  #standard_PHD_balance <- s_PHD[50000]
  
  #standard_HIFd_HRE_active_balance <- s_HIFd_HRE_active[50000]
  #standard_HIFd_HRE_balance <- s_HIFd_HRE[50000]
  #standard_HIFa_balance <- s_HIFa[50000]
  #standard_HIFa_active_balance <- s_HIFa_active[50000]
  
#plot 
  time <- seq(0, by=dt/3600,length.out = timesteps)
  #plot of "HIF¦Á & HIF-dimer"
  plot(time,s_HIFa_active,col ="green", type = "l",xlab = "time(h)",ylab = "concentration(nM)")
  lines(time,s_HIFd_HRE_active,col="red")
  lines(time,s_HIFa,col ="blue")
  lines(time,s_HIFd_HRE,col = "grey")
  abline(h=standard_HIFd_HRE_active_balance, col="red",lty=2)
  abline(h=standard_HIFa_balance,col="blue",lty=2)
  abline(h=standard_HIFd_HRE_balance, col = "grey", lty=2)
  abline(h=standard_HIFa_active_balance, col="green", lty=2)
  legend("right",legend = c("HIF¦Á","HIF¦Á + HIF¦Á_pOH","HIF_HRE","HIF_HRE + HIF_pOH_HRE"),col = c("blue","green","grey","red"),lty = 1, cex = 0.8)
  
  #plot of "protein & PHD"
  plot(time,s_protein,col ="blue", type = "l",xlab = "time(h)",ylab = "concentration(nM)")
  lines(time,s_PHD,col ="red",xlab = "time(h)",ylab = "concentration(nM)")
  abline(h = standard_PHD_balance,col = "red",lty = 2)
  abline(h = standard_protein_balance,col = "blue",lty = 2)
  legend("topright",legend = c("protein","PHD"),col = c("blue","red"),lty = 1, cex = 1)
  
  
   
#additional£ºgraph with ggplot
  #create a dataframe
  #HIF_df <- data.frame(time, s_HIFa,s_HIFa_pOH,
                       #s_HIFb,s_PHD,s_VHL,s_HIFd,
                       #s_HRE,HIFd_HRE,s_mRNA,s_protein)
  
  
  #ggplot
  #library(ggplot2)
  #plot of "HIF¦Á & HIF-dimer"
  #ggplot(HIF_df)+
    #geom_line(aes(x=time,y=s_HIFa,col="s_HIFa"))+
    #geom_line(aes(x=time,y=s_HIFd,col="s_HIFd"))+
    #scale_colour_discrete(breaks = c("s_HIFa","s_HIFd"), labels = c('HIF¦Á','HIF-dimer'))+
    #labs(x="time(h)",y="concentration(nM)")
  #plot of "protein & PHD"
  #ggplot(HIF_df)+
    #geom_line(aes(x=time,y=s_protein,col="s_protein"))+
    #geom_line(aes(x=time,y=s_PHD,col="s_PHD"))+
    #scale_colour_discrete(breaks = c("s_protein","s_PHD"), labels = c('protein','PHD'))+
    #labs(x="time(h)",y="concentration(nM)")
  
