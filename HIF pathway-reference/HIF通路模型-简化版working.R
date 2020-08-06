#initial conditions

#velocity constant
  k1 <- 0.005
  k2 <- 0.0002
  k3 <- 0.045
  km3a <- 250000 
  km3b <- 100
  # reference:O2,Km3a&Km5a are in ¦ÌM,others are in nM
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
  HIF1a <- 0.05
  O2 <- 100000
  HIF1a_pOH <- 0
  HIFb <- 170
  PHD <- 100
  VHL <- 50
  HIFd <- 0
  HRE <- 50
  HIFd_HRE <- 0
  mRNA <- 0
  EGFP <- 0 # same as protein in reference
  
#set time interval [unit="s"]
  dt <- 10
  
#set loop to simulate the pathway
  #set timelength
  timesteps <- 50000
  #set empty vectors to store the variables
  s_HIF1a <- numeric(length = timesteps) #equilibrium around 0.4
  s_HIF1a[1] <- HIF1a
  s_HIF1a_pOH <- numeric(length = timesteps)
  s_HIF1a_pOH[1] <- HIF1a_pOH
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
  s_EGFP <- numeric(length = timesteps)
  s_EGFP[1] <- EGFP 
  
#set "for" loop
for (i in c(2:timesteps-1)){
 #equation for velocity
  v1 <- k1
  v2 <- k2*HIF1a
  v3 <- k3*PHD*(O2/(km3a+O2))*(HIF1a/(km3b+HIF1a))
  v4 <- k4*VHL*(HIF1a_pOH/(km4+HIF1a_pOH))
  v5 <- k5*HIF1a*HIFb
  v6 <- k6f*HIFd*HRE-k6r*HIFd_HRE
  v7 <- k7*PHD
  #v8&v9 currently empty
  v10 <- k10*mRNA
  v11 <- k11*EGFP
  v12 <- k12*HIFd_HRE
  v13 <- k13*HIFd
  v14 <- k14*HIFd_HRE
  v15 <- k15*mRNA
 #simulate the process
  HIF1a <- HIF1a + (v1-v2-v3-v5+v13)*dt
  HIF1a_pOH <- HIF1a_pOH + (v3-v4)*dt
  HIFb <- HIFb + (v13-v5)*dt
  PHD <- PHD + (v12-v7)*dt
  VHL <- VHL
  HIFd <- HIFd + (v5-v13-v6)*dt
  HRE <- HRE + (-v6)*dt
  HIFd_HRE <- HIFd_HRE + v6*dt
  mRNA <- mRNA + (v14-v15)*dt
  EGFP <- EGFP + (v10-v11)*dt
  s_HIF1a[i+1] <- HIF1a
  s_HIF1a_pOH[i+1] <- HIF1a_pOH
  s_HIFb[i+1] <- HIFb
  s_PHD[i+1] <- PHD
  s_VHL[i+1] <- VHL
  s_HIFd[i+1] <- HIFd
  s_HRE[i+1] <- HRE
  s_HIFd_HRE[i+1] <- HIFd_HRE
  s_mRNA[i+1] <- mRNA
  s_EGFP[i+1] <- EGFP
}

#plot 
  time <- seq(0, by=dt/3600,length.out = timesteps)
  #plot of "HIF1¦Á & HIF-dimer"
  plot(time,s_HIF1a,col ="blue", type = "l",xlab = "time(h)",ylab = "concentration(nM)")
  lines(time,s_HIFd_HRE,col="red")
  legend("topleft",legend = c("HIF1¦Á","HIF-HRE"),col = c("blue","red"),lty = 1, cex = 1)
  #plot of "EGFP & PHD"
  plot(time,s_EGFP,col ="blue", type = "l",xlab = "time(h)",ylab = "concentration(nM)")
  lines(time,s_PHD,col ="red",xlab = "time(h)",ylab = "concentration(nM)")
  legend("topleft",legend = c("EGFP","PHD"),col = c("blue","red"),lty = 1, cex = 1)
 
  
   
#additional£ºgraph with ggplot
  #create a dataframe
  HIF_df <- data.frame(time, s_HIF1a,s_HIF1a_pOH,
                       s_HIFb,s_PHD,s_VHL,s_HIFd,
                       s_HRE,HIFd_HRE,s_mRNA,s_EGFP)
  
  
  #ggplot
  library(ggplot2)
  #plot of "HIF1¦Á & HIF-dimer"
  ggplot(HIF_df)+
    geom_line(aes(x=time,y=s_HIF1a,col="s_HIF1a"))+
    geom_line(aes(x=time,y=s_HIFd,col="s_HIFd"))+
    scale_colour_discrete(breaks = c("s_HIF1a","s_HIFd"), labels = c('HIF1¦Á','HIF-dimer'))+
    labs(x="time(h)",y="concentration(nM)")
  #plot of "EGFP & PHD"
  ggplot(HIF_df)+
    geom_line(aes(x=time,y=s_EGFP,col="s_EGFP"))+
    geom_line(aes(x=time,y=s_PHD,col="s_PHD"))+
    scale_colour_discrete(breaks = c("s_EGFP","s_PHD"), labels = c('EGFP','PHD'))+
    labs(x="time(h)",y="concentration(nM)")
  
