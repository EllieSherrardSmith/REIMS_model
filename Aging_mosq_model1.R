mosq_age_model = odin::odin({
  
  ### Compartments ###
  
  ## susceptible mosquitoes
  deriv(S1) <- emergence - rho1 * alpha * S1 - (1 - rho1)*alpha*S1 - mu*S1 #BORN susceptible Age up to 3 days old
  deriv(S2) <- rho1 * alpha * S1 - rho2 * alpha * S2 - (1 - rho2)*alpha*S2 - mu*S2  #Susceptible aged 4 - 6 days
  deriv(S3) <- rho2 * alpha * S2 - rho3 * alpha * S3 - (1 - rho3)*alpha*S3 - mu*S3  #Susceptible aged 7 - 9 days
  deriv(S4) <- rho3 * alpha * S3 - rho4 * alpha * S4 - (1 - rho4)*alpha*S4 - mu*S4  #Susceptible aged 10 - 12 days
  deriv(S5) <- rho4 * alpha * S4 - rho5 * alpha * S5 - (1 - rho5)*alpha*S5 - mu*S5  #Susceptible aged 13 - 15 days
  deriv(S6) <- rho5 * alpha * S5 - rho6 * alpha * S6 - (1 - rho6)*alpha*S6 - mu*S6  #Susceptible aged 16 - 19 days
  deriv(S7) <- rho6 * alpha * S6 - rho7 * alpha * S7 - (1 - rho7)*alpha*S7 - mu*S7  #Susceptible aged 16 - 19 days
  
  SUSCEPTIBLES = S1 + S2 + S3 + S4 + S5 + S6 + S7
  
  ### exposed mosquitoes
  
  deriv(E1) <- (1 - rho1)*alpha*S1 - nu*E1 - mu*E1  #exposed mosquitoes (aged 4-6 days)
  deriv(E2) <- (1 - rho2)*alpha*S2 - nu*E2 - mu*E2  #exposed mosquitoes (aged 4-6 days)
  deriv(E3) <- (1 - rho3)*alpha*S3 - nu*E3 - mu*E3  #exposed mosquitoes (aged 4-6 days)
  deriv(E4) <- (1 - rho4)*alpha*S4 - nu*E4 - mu*E4  #exposed mosquitoes (aged 4-6 days)
  deriv(E5) <- (1 - rho5)*alpha*S5 - nu*E5 - mu*E5  #exposed mosquitoes (aged 4-6 days)
  deriv(E6) <- (1 - rho6)*alpha*S6 - nu*E6 - mu*E6  #exposed mosquitoes (aged 4-6 days)
  
  EXPOSED = E1 + E2 + E3 + E4 + E5 + E6
  
  ### infectious mosquito

  deriv(I1) <- nu*E4 - (mu + delta)*I1 #exposed mosquitoes (aged 4-6 days)
  deriv(I2) <- nu*E5 - (mu + delta)*I2 #exposed mosquitoes (aged 4-6 days)
  deriv(I3) <- nu*E6 - (mu + delta)*I3 #exposed mosquitoes (aged 4-6 days)
 
  INFECTIOUS = I1 + I2 + I3
  NMOSQ = SUSCEPTIBLES + EXPOSED + INFECTIOUS
  
  
  output(SUSCEPTIBLES_MOSQ) <- SUSCEPTIBLES
  output(EXPOSED_MOSQ) <- EXPOSED
  output(INFECTIOUS_MOSQ) <- INFECTIOUS
  
  output(TOTAL_MOSQ) <- NMOSQ
  
  
  
  # Number of mosquitoes born (depends on PL, number of larvae), or is constant outside of seasonality
  emergence <- 0.5*PL/dPL
  
    ## The force of infection can be defined for each rho[x]
  FOIv <- user() # ultimately this is from the transmission model and is variable, for now, fixed
  nu <- user() # Extrinsic incubation period. (nu)
  delta <- user() # increased death rate due to risky behaviour from infectious mosquitoes
  alpha <- user() # a blood meal each 3 days
  
  ## Research whether mosquitoes bite at different rates at different ages (when uninfected)
  ## There is evidence of infected mosquitoes biting more often...
  rho1 <- 1 - FOIv
  rho2 <- 1 - FOIv
  rho3 <- 1 - FOIv
  rho4 <- 1 - FOIv
  rho5 <- 1 - FOIv
  rho6 <- 1 - FOIv
  rho7 <- 1 - FOIv
  

  # initial conditions of the variables
  
  initial(S1) <- 1850
  initial(S2) <- 0
  initial(S3) <- 0   # num susceptibles
  initial(S4) <- 0
  initial(S5) <- 0
  initial(S6) <- 0
  initial(S7) <- 0
  
  initial(E1) <- 100
  initial(E2) <- 0
  initial(E3) <- 0   # num exposed
  initial(E4) <- 0
  initial(E5) <- 0
  initial(E6) <- 0

  initial(I1) <- 50
  initial(I2) <- 0
  initial(I3) <- 0   

  # parameter values
  
  ##------------------------------------------------------------------------------
  ###################
  ## LARVAL STATES ##
  ###################
  ##------------------------------------------------------------------------------
  
  # Model by White et al.
  # (https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153)
  
  # EL - early larval instar stage
  # LL - late larval instar stage
  # PL - pupal stage
  
  # mean carrying capacity from initial mosquito density:
  dLL <- user() # development time of larvae
  dPL <- user() #development time of pupae
  dEL <- user() #development time of early stage
  muLL <- user() #daily density dep. mortality rate of larvae
  muPL <- user() #daily den. dep. mortality rate of pupae
  muEL <- user() #daily den. dep. mortality rate of early stage
  gammaL <- user() # eff. of den. dep. on late stage relative to early stage
  
  # fitted entomological parameters:
  mv0 <- user() # initial mosquito density
  mu0 <- user() # baseline mosquito death rate
  tau1 <- user() # duration of host-seeking behaviour
  tau2 <- user() # duration of resting behaviour
  betaL <- user() # maximum number of eggs per oviposition per mosq
  
  # Entomological variables:
  p2 <- exp(-mu0 * tau2)  # probability of surviving one resting cycle
  p10 <- exp(-mu0 * tau1)  # probability of surviving one feeding cycle
  p1 <- p10 ## would update with interventions
  
  fv0 <- 1 / (tau1 + tau2)
  eov <- betaL/mu
  beta_larval <- ((INFECTIOUS + E3 + E4 + E5 + E6)/NMOSQ) * (eov*mu) # Number of eggs laid per day adjusted by blood fed and lagged i.e. those that can lay
  b_lambda <- (gammaL*muLL/muEL-dEL/dLL+(gammaL-1)*muLL*dEL)
  lambda <- -0.5*b_lambda + sqrt(0.25*b_lambda^2 + gammaL*beta_larval*muLL*dEL/(2*muEL*mu0*dLL*(1+dPL*muPL)))
  K0 <- 2*mv0*dLL*mu0*(1+dPL*muPL)*gammaL*(lambda+1)/(lambda/(muLL*dEL)-1/(muLL*dLL)-1)
  
  
  ##########################
  ## SEASONALITY FUNCTION ##
  ##########################
  ##------------------------------------------------------------------------------
  
  # Seasonality is added into the model using a Fourier series that was fit to rainfall at every admin 1 level
  pi <- user() # weird quirk, need to pass pi
  
  # The parameters for the fourier series
  ssa0 <- user()
  ssa1 <- user()
  ssa2 <- user()
  ssa3 <- user()
  ssb1 <- user()
  ssb2 <- user()
  ssb3 <- user()
  theta_c <- user()
  # Recreation of the rainfall function
  theta2 <- if(ssa0 == 0 && ssa1  == 0 && ssa2  == 0 && ssb1  == 0 && ssb2  == 0 && ssb3  == 0 && theta_c  == 0)
    1 else max((ssa0+ssa1*cos(2*pi*t/365)+ssa2*cos(2*2*pi*t/365)+ssa3*cos(3*2*pi*t/365)+ssb1*sin(2*pi*t/365)+ssb2*sin(2*2*pi*t/365)+ ssb3*sin(3*2*pi*t/365) ) /theta_c,0.001)
  
  
  # Seasonal carrying capacity KL = base carrying capacity K0 * effect for time of year theta:
  KL <- K0*theta2
  # fv <- 1/( tau1/(1-zbar) + tau2 ) # mosquito feeding rate (zbar from intervention param.)
  mu <- -fv0*log(p1*p2) # mosquito death rate
  
  # finding equilibrium and initial values for EL, LL & PL
  init_PL <- user()
  initial(PL) <- init_PL
  init_LL <- user()
  initial(LL) <- init_LL
  init_EL <- user()
  initial(EL) <- init_EL
  
  # (beta_larval (egg rate) * total mosquito) - den. dep. egg mortality - egg hatching
  deriv(EL) <- beta_larval*NMOSQ-muEL*(1+(EL+LL)/KL)*EL - EL/dEL
  # egg hatching - den. dep. mortality - maturing larvae
  deriv(LL) <- EL/dEL - muLL*(1+gammaL*(EL + LL)/KL)*LL - LL/dLL
  # pupae - mortality - fully developed pupae
  deriv(PL) <- LL/dLL - muPL*PL - PL/dPL
  
})




mod <- mosq_age_model$new(    ## The force of infection can be defined for each rho[x]
  FOIv = 1/5, # ultimately this is from the transmission model and is variable, for now, fixed
  nu = 10, # Extrinsic incubation period. (nu) DAYS
  delta = 1/8, # increased death rate due to risky behaviour from infectious mosquitoes
  alpha = 1/3, # a blood meal each 3 days
  
  muEL = 0.0338,#daily den. dep. mortality rate of early stage
  muLL = 0.0348,#daily density dep. mortality rate of larvae
  muPL = 0.249,#daily den. dep. mortality rate of pupae
  dEL = 6.64,#development time of early stage
  dLL = 3.72,# development time of larvae
  dPL = 0.643,#development time of pupae
  gammaL = 13.25,# eff. of den. dep. on late stage relative to early stage
  # km = 11,
  # cm = 0.05,
  betaL = 21.2,# maximum number of eggs per oviposition per mosq
  
  # fitted entomological parameters:
  mv0 = 2000, # initial mosquito density
  mu0 = 1/12, # baseline mosquito death rate
  tau1 = 0.69, # duration of host-seeking behaviour
  tau2 = 2.31, # duration of resting behaviour
  
  pi = 3.141593, # weird quirk, need to pass pi
  
  # The parameters for the fourier series
  ssa0 = 0.2852297,					
  ssa1 = -0.2952712,
  ssa2 = -0.03408224,	
  ssa3 = 0.07596435,
  ssb1 = -0.1126063,
  ssb2 = 0.07789561,
  ssb3 = -0.007051094,
  theta_c = 0.2852297,
  
  init_PL = 600,
  init_LL = 3000,
  init_EL = 1200
)#,


## Run the model for a series of times from 0 to 10:
t <- seq(0, 3650, length.out = 365*10)
y <- mod$run(t)
dim(y) ## check the dimensions; column 2 == S; column 3 == I; column 4 == R etc

# head(y) ## so you can see the outputs


par(mfrow=c(1,3))
prev_infectious = y[,23]/y[,24] ## infectious
prev_suscept = y[,21]/y[,24] ## susceptibles
prev_exposed = y[,22]/y[,24] ## exposed
plot(prev_infectious ~ t,type="l",col="darkred",xlim=c(1800,3000),yaxt="n",
     ylab="Infectious mosquitoes (%)",ylim=c(0,1),
     xlab = "Time in days")
axis(2,las=2,at=c(0.05,0.10,0.15),labels=c(5,10,15))

plot(prev_suscept ~ t,type="l",col="blue",xlim=c(1800,3000),yaxt="n",
     ylab="Susceptible mosquitoes (%)",ylim=c(0,1),
     xlab = "Time in days")
axis(2,las=2,at=c(0.85,0.90,0.95),labels=c(85,90,95))

age_0_3 = y[,2]
age_4_6 = y[,3]+y[,9]

prop_age_under_6 = (age_0_3 + age_4_6)/y[,24]

plot(prop_age_under_6 ~ t,type="l",col="darkblue",xlim=c(1800,3000),yaxt="n",
     ylab="Proportion mosquitoes (%)",ylim=c(0,1),
     xlab = "Time in days")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))



age_7_9 = y[,4]+y[,10]
age_10_12 = y[,5]+y[,11]
age_13_15 = y[,6]+y[,12]+y[,15]
age_16_18 = y[,7]+y[,13]+y[,16]
age_19_21 = y[,8]+y[,14]+y[,17]

prop_age_over_13 = (age_13_15 + age_16_18 + age_19_21)/y[,24]
lines(prop_age_over_13 ~ t,col="red")

legend("topleft",legend = c("up to 6 days old",
                            "over 12 days old"),
       col=c("darkblue","red"),lty=1,bty="n")


par(mfrow=c(2,3))
td = 2000
barplot(c(age_0_3[td],age_4_6[td],age_7_9[td],age_10_12[td],age_13_15[td],age_16_18[td],age_19_21[td]),
        xlab="Age group (days)",col = adegenet::transp("darkblue",1))
axis(1,at=seq(0.75,8,length=7),labels=c("0-3","4-6","7-9","10-12","13-15","16-19","20-21"))

td = 2060
barplot(c(age_0_3[td],age_4_6[td],age_7_9[td],age_10_12[td],age_13_15[td],age_16_18[td],age_19_21[td]),
        xlab="Age group (days)",col = adegenet::transp("darkblue",0.8))
axis(1,at=seq(0.75,8,length=7),labels=c("0-3","4-6","7-9","10-12","13-15","16-19","20-21"))

td = 2121
barplot(c(age_0_3[td],age_4_6[td],age_7_9[td],age_10_12[td],age_13_15[td],age_16_18[td],age_19_21[td]),
        xlab="Age group (days)",col = adegenet::transp("darkblue",0.6))
axis(1,at=seq(0.75,8,length=7),labels=c("0-3","4-6","7-9","10-12","13-15","16-19","20-21"))

td = 2160
barplot(c(age_0_3[td],age_4_6[td],age_7_9[td],age_10_12[td],age_13_15[td],age_16_18[td],age_19_21[td]),
        xlab="Age group (days)",col = adegenet::transp("darkblue",0.4))
axis(1,at=seq(0.75,8,length=7),labels=c("0-3","4-6","7-9","10-12","13-15","16-19","20-21"))


td = 2243
barplot(c(age_0_3[td],age_4_6[td],age_7_9[td],age_10_12[td],age_13_15[td],age_16_18[td],age_19_21[td]),
        xlab="Age group (days)",col = adegenet::transp("darkblue",0.3))
axis(1,at=seq(0.75,8,length=7),labels=c("0-3","4-6","7-9","10-12","13-15","16-19","20-21"))

td = 2304
barplot(c(age_0_3[td],age_4_6[td],age_7_9[td],age_10_12[td],age_13_15[td],age_16_18[td],age_19_21[td]),
        xlab="Age group (days)",col = adegenet::transp("darkblue",0.2))
axis(1,at=seq(0.75,8,length=7),labels=c("0-3","4-6","7-9","10-12","13-15","16-19","20-21"))




par(mfrow=c(1,1))
plot(prop_age_under_6 ~ t,type="l",col="darkblue",xlim=c(1800,2400),yaxt="n",
     ylab="Proportion mosquitoes (%)",ylim=c(0,1),
     xlab = "Time in days")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
lines(prop_age_over_13 ~ t,col="red")

legend("topleft",legend = c("up to 6 days old",
                            "over 12 days old"),
       col=c("darkblue","red"),lty=1,bty="n")

tds = c(2000,2060,2121,2160,2243)
cols_tds = adegenet::transp("darkblue",c(1,0.8,0.6,0.4,0.3))
for(i in 1:length(tds)){
  abline(v=tds[i],lty=2,col=cols_tds[i])
}










#########################################
##
## Explore 1: increased mortality

mod_delta_1 <- mosq_age_model$new(    ## The force of infection can be defined for each rho[x]
  FOIv = 1/5, # ultimately this is from the transmission model and is variable, for now, fixed
  nu = 10, # Extrinsic incubation period. (nu) DAYS
  delta = 1/8, # increased death rate due to risky behaviour from infectious mosquitoes
  alpha = 1/3, # a blood meal each 3 days
  
  muEL = 0.0338,#daily den. dep. mortality rate of early stage
  muLL = 0.0348,#daily density dep. mortality rate of larvae
  muPL = 0.249,#daily den. dep. mortality rate of pupae
  dEL = 6.64,#development time of early stage
  dLL = 3.72,# development time of larvae
  dPL = 0.643,#development time of pupae
  gammaL = 13.25,# eff. of den. dep. on late stage relative to early stage
  # km = 11,
  # cm = 0.05,
  betaL = 21.2,# maximum number of eggs per oviposition per mosq
  
  # fitted entomological parameters:
  mv0 = 2000, # initial mosquito density
  mu0 = 1/12, # baseline mosquito death rate
  tau1 = 0.69, # duration of host-seeking behaviour
  tau2 = 2.31, # duration of resting behaviour
  
  pi = 3.141593, # weird quirk, need to pass pi
  
  # The parameters for the fourier series
  ssa0 = 0.2852297,					
  ssa1 = -0.2952712,
  ssa2 = -0.03408224,	
  ssa3 = 0.07596435,
  ssb1 = -0.1126063,
  ssb2 = 0.07789561,
  ssb3 = -0.007051094,
  theta_c = 0.2852297,
  
  init_PL = 600,
  init_LL = 3000,
  init_EL = 1200
)#



y_delta0 <- y                  ## 1/12
y_delta1 <- mod_delta_1$run(t) ## 1/10
y_delta2 <- mod_delta_1$run(t) ## 1/8

par(mfrow=c(1,1))
prev_infectious_0 = y_delta0[,23]/y_delta0[,24] ## infectious
prev_infectious_1 = y_delta1[,23]/y_delta1[,24] ## infectious
prev_infectious_2 = y_delta2[,23]/y_delta2[,24] ## infectious

plot(prev_infectious_0 ~ t,ylim = c(0,0.3),type="l",col="darkred",xlim=c(1800,3000))
lines(prev_infectious_1 ~ t, lty=2,col="darkred")
lines(prev_infectious_2 ~ t, lty=3,col="darkred")

## y_delta0
age_0_3 = y[,2]
age_4_6 = y[,3]+y[,9]

prop_age_under_6 = (age_0_3 + age_4_6)/y[,24]

age_7_9 = y[,4]+y[,10]
age_10_12 = y[,5]+y[,11]
age_13_15 = y[,6]+y[,12]+y[,15]
age_16_18 = y[,7]+y[,13]+y[,16]
age_19_21 = y[,8]+y[,14]+y[,17]

prop_age_over_13 = (age_13_15 + age_16_18 + age_19_21)/y_delta0[,24]

## y_delta1
age_0_3_d1 = y_delta1[,2]
age_4_6_d1 = y_delta1[,3]+y_delta1[,9]

prop_age_under_6_d1 = (age_0_3_d1 + age_4_6_d1)/y_delta1[,24]

age_7_9_d1 = y_delta1[,4]+y_delta1[,10]
age_10_12_d1 = y_delta1[,5]+y_delta1[,11]
age_13_15_d1 = y_delta1[,6]+y_delta1[,12]+y_delta1[,15]
age_16_18_d1 = y_delta1[,7]+y_delta1[,13]+y_delta1[,16]
age_19_21_d1 = y_delta1[,8]+y_delta1[,14]+y_delta1[,17]

prop_age_over_13_d1 = (age_13_15_d1 + age_16_18_d1 + age_19_21_d1)/y_delta1[,24]


## y_delta2
age_0_3_d2 = y_delta2[,2]
age_4_6_d2 = y_delta2[,3]+y_delta2[,9]

prop_age_under_6_d2 = (age_0_3_d2 + age_4_6_d2)/y_delta2[,24]

age_7_9_d2 = y_delta2[,4]+y_delta2[,10]
age_10_12_d2 = y_delta2[,5]+y_delta2[,11]
age_13_15_d2 = y_delta2[,6]+y_delta2[,12]+y_delta2[,15]
age_16_18_d2 = y_delta2[,7]+y_delta2[,13]+y_delta2[,16]
age_19_21_d2 = y_delta2[,8]+y_delta2[,14]+y_delta2[,17]

prop_age_over_13_d2 = (age_13_15_d2 + age_16_18_d2 + age_19_21_d2)/y_delta2[,24]

td = 2121
barplot(c(age_0_3[td],age_0_3_d1[td],age_0_3_d2[td],NA,
          age_4_6[td],age_4_6_d1[td],age_4_6_d2[td],NA,
          age_7_9[td],age_7_9_d1[td],age_7_9_d2[td],NA,
          age_10_12[td],age_10_12_d1[td],age_10_12_d2[td],NA,
          age_13_15[td],age_13_15_d1[td],age_13_15_d2[td],NA,
          age_16_18[td],age_16_18_d1[td],age_16_18_d2[td],NA,
          age_19_21[td],age_19_21_d1[td],age_19_21_d2[td]),
        xlab="Age group (days)",col = c(adegenet::transp("darkblue",c(1,0.65,0.25)),NA),
        ylim=c(0,400))
axis(1,at=seq(2,31,length=7),labels=c("0-3","4-6","7-9","10-12","13-15","16-19","20-21"))

legend("topright",
       legend = c("background mortality 1/12",
                  "background mortality 1/10",
                  "background mortality 1/8"),
       col = c(adegenet::transp("darkblue",c(1,0.65,0.25)),NA),
       pch=15,bty="n")

td = 2160
barplot(c(age_0_3[td],age_0_3_d1[td],age_0_3_d2[td],NA,
          age_4_6[td],age_4_6_d1[td],age_4_6_d2[td],NA,
          age_7_9[td],age_7_9_d1[td],age_7_9_d2[td],NA,
          age_10_12[td],age_10_12_d1[td],age_10_12_d2[td],NA,
          age_13_15[td],age_13_15_d1[td],age_13_15_d2[td],NA,
          age_16_18[td],age_16_18_d1[td],age_16_18_d2[td],NA,
          age_19_21[td],age_19_21_d1[td],age_19_21_d2[td]),
        xlab="Age group (days)",col = c(adegenet::transp("darkblue",c(1,0.65,0.25)),NA),
        ylim=c(0,50))
axis(1,at=seq(2,31,length=7),labels=c("0-3","4-6","7-9","10-12","13-15","16-19","20-21"))

# plot(y[1800:3000,24] ~ t[1800:3000])
# lines(1000*prop_age_under_6 ~ t,col="blue")
# lines(1000*prop_age_over_13 ~ t,col="red")




#########################################
##
## Explore 2: increased time between bites

mod_delta_1 <- mosq_age_model$new(    ## The force of infection can be defined for each rho[x]
  FOIv = 1/5, # ultimately this is from the transmission model and is variable, for now, fixed
  nu = 10, # Extrinsic incubation period. (nu) DAYS
  delta = 1/8, # increased death rate due to risky behaviour from infectious mosquitoes
  alpha = 1/3, # a blood meal each 3 days, or 4 or 5
  
  muEL = 0.0338,#daily den. dep. mortality rate of early stage
  muLL = 0.0348,#daily density dep. mortality rate of larvae
  muPL = 0.249,#daily den. dep. mortality rate of pupae
  dEL = 6.64,#development time of early stage
  dLL = 3.72,# development time of larvae
  dPL = 0.643,#development time of pupae
  gammaL = 13.25,# eff. of den. dep. on late stage relative to early stage
  # km = 11,
  # cm = 0.05,
  betaL = 21.2,# maximum number of eggs per oviposition per mosq
  
  # fitted entomological parameters:
  mv0 = 2000, # initial mosquito density
  mu0 = 1/12, # baseline mosquito death rate
  tau1 = 0.69, # duration of host-seeking behaviour ## 0.69, 0.8, 1
  tau2 = 2.31, # duration of resting behaviour      ## 2.31, 3.2, 4 
  
  pi = 3.141593, # weird quirk, need to pass pi
  
  # The parameters for the fourier series
  ssa0 = 0.2852297,					
  ssa1 = -0.2952712,
  ssa2 = -0.03408224,	
  ssa3 = 0.07596435,
  ssb1 = -0.1126063,
  ssb2 = 0.07789561,
  ssb3 = -0.007051094,
  theta_c = 0.2852297,
  
  init_PL = 600,
  init_LL = 3000,
  init_EL = 1200
)#



y_delta0 <- y                  ## 1/12
y_delta1 <- mod_delta_1$run(t) ## 0.8, 3.2 and 4
y_delta2 <- mod_delta_1$run(t) ## 1, 4 and 5

par(mfrow=c(1,1))
prev_infectious_0 = y_delta0[,23]/y_delta0[,24] ## infectious
prev_infectious_1 = y_delta1[,23]/y_delta1[,24] ## infectious
prev_infectious_2 = y_delta2[,23]/y_delta2[,24] ## infectious

plot(prev_infectious_0 ~ t,ylim = c(0,0.3),type="l",col="darkred",xlim=c(1800,3000))
lines(prev_infectious_1 ~ t, lty=2,col="darkred")
lines(prev_infectious_2 ~ t, lty=3,col="darkred")

## y_delta0
age_0_3 = y[,2]
age_4_6 = y[,3]+y[,9]

prop_age_under_6 = (age_0_3 + age_4_6)/y[,24]

age_7_9 = y[,4]+y[,10]
age_10_12 = y[,5]+y[,11]
age_13_15 = y[,6]+y[,12]+y[,15]
age_16_18 = y[,7]+y[,13]+y[,16]
age_19_21 = y[,8]+y[,14]+y[,17]

prop_age_over_13 = (age_13_15 + age_16_18 + age_19_21)/y_delta0[,24]

## y_delta1
age_0_3_d1 = y_delta1[,2]
age_4_6_d1 = y_delta1[,3]+y_delta1[,9]

prop_age_under_6_d1 = (age_0_3_d1 + age_4_6_d1)/y_delta1[,24]

age_7_9_d1 = y_delta1[,4]+y_delta1[,10]
age_10_12_d1 = y_delta1[,5]+y_delta1[,11]
age_13_15_d1 = y_delta1[,6]+y_delta1[,12]+y_delta1[,15]
age_16_18_d1 = y_delta1[,7]+y_delta1[,13]+y_delta1[,16]
age_19_21_d1 = y_delta1[,8]+y_delta1[,14]+y_delta1[,17]

prop_age_over_13_d1 = (age_13_15_d1 + age_16_18_d1 + age_19_21_d1)/y_delta1[,24]


## y_delta2
age_0_3_d2 = y_delta2[,2]
age_4_6_d2 = y_delta2[,3]+y_delta2[,9]

prop_age_under_6_d2 = (age_0_3_d2 + age_4_6_d2)/y_delta2[,24]

age_7_9_d2 = y_delta2[,4]+y_delta2[,10]
age_10_12_d2 = y_delta2[,5]+y_delta2[,11]
age_13_15_d2 = y_delta2[,6]+y_delta2[,12]+y_delta2[,15]
age_16_18_d2 = y_delta2[,7]+y_delta2[,13]+y_delta2[,16]
age_19_21_d2 = y_delta2[,8]+y_delta2[,14]+y_delta2[,17]

prop_age_over_13_d2 = (age_13_15_d2 + age_16_18_d2 + age_19_21_d2)/y_delta2[,24]

td = 2121
barplot(c(age_0_3[td],age_0_3_d1[td],age_0_3_d2[td],NA,
          age_4_6[td],age_4_6_d1[td],age_4_6_d2[td],NA,
          age_7_9[td],age_7_9_d1[td],age_7_9_d2[td],NA,
          age_10_12[td],age_10_12_d1[td],age_10_12_d2[td],NA,
          age_13_15[td],age_13_15_d1[td],age_13_15_d2[td],NA,
          age_16_18[td],age_16_18_d1[td],age_16_18_d2[td],NA,
          age_19_21[td],age_19_21_d1[td],age_19_21_d2[td]),
        xlab="Age group (days)",col = c(adegenet::transp("darkblue",c(1,0.65,0.25)),NA),
        ylim=c(0,400))
axis(1,at=seq(2,31,length=7),labels=c("0-3","4-6","7-9","10-12","13-15","16-19","20-21"))

legend("topright",
       legend = c("biting rate 1 every 3 days",
                  "biting rate 1 every 4 days",
                  "biting rate 1 every 5 days"),
       col = c(adegenet::transp("darkblue",c(1,0.65,0.25)),NA),
       pch=15,bty="n")

td = 2160
barplot(c(age_0_3[td],age_0_3_d1[td],age_0_3_d2[td],NA,
          age_4_6[td],age_4_6_d1[td],age_4_6_d2[td],NA,
          age_7_9[td],age_7_9_d1[td],age_7_9_d2[td],NA,
          age_10_12[td],age_10_12_d1[td],age_10_12_d2[td],NA,
          age_13_15[td],age_13_15_d1[td],age_13_15_d2[td],NA,
          age_16_18[td],age_16_18_d1[td],age_16_18_d2[td],NA,
          age_19_21[td],age_19_21_d1[td],age_19_21_d2[td]),
        xlab="Age group (days)",col = c(adegenet::transp("darkblue",c(1,0.65,0.25)),NA),
        ylim=c(0,50))
axis(1,at=seq(2,31,length=7),labels=c("0-3","4-6","7-9","10-12","13-15","16-19","20-21"))

# plot(y[1800:3000,24] ~ t[1800:3000])
# lines(1000*prop_age_under_6 ~ t,col="blue")
# lines(1000*prop_age_over_13 ~ t,col="red")






#########################################
##
## Explore 2: increased time between bites

mod_delta_1 <- mosq_age_model$new(    ## The force of infection can be defined for each rho[x]
  FOIv = 1/5, # ultimately this is from the transmission model and is variable, for now, fixed 1/3, 1/5, 1/8
  nu = 10, # Extrinsic incubation period. (nu) DAYS
  delta = 1/8, # increased death rate due to risky behaviour from infectious mosquitoes
  alpha = 1/3, # a blood meal each 3 days, or 4 or 5
  
  muEL = 0.0338,#daily den. dep. mortality rate of early stage
  muLL = 0.0348,#daily density dep. mortality rate of larvae
  muPL = 0.249,#daily den. dep. mortality rate of pupae
  dEL = 6.64,#development time of early stage
  dLL = 3.72,# development time of larvae
  dPL = 0.643,#development time of pupae
  gammaL = 13.25,# eff. of den. dep. on late stage relative to early stage
  # km = 11,
  # cm = 0.05,
  betaL = 21.2,# maximum number of eggs per oviposition per mosq
  
  # fitted entomological parameters:
  mv0 = 2000, # initial mosquito density
  mu0 = 1/12, # baseline mosquito death rate
  tau1 = 0.69, # duration of host-seeking behaviour ## 0.69, 0.8, 1
  tau2 = 2.31, # duration of resting behaviour      ## 2.31, 3.2, 4 
  
  pi = 3.141593, # weird quirk, need to pass pi
  
  # The parameters for the fourier series
  ssa0 = 0.2852297,					
  ssa1 = -0.2952712,
  ssa2 = -0.03408224,	
  ssa3 = 0.07596435,
  ssb1 = -0.1126063,
  ssb2 = 0.07789561,
  ssb3 = -0.007051094,
  theta_c = 0.2852297,
  
  init_PL = 600,
  init_LL = 3000,
  init_EL = 1200
)#



y_delta0 <- y                  ## 1/5
y_delta1 <- mod_delta_1$run(t) ## 1/3
y_delta2 <- mod_delta_1$run(t) ## 1/8

par(mfrow=c(1,1))
prev_infectious_0 = y_delta0[,23]/y_delta0[,24] ## infectious
prev_infectious_1 = y_delta1[,23]/y_delta1[,24] ## infectious
prev_infectious_2 = y_delta2[,23]/y_delta2[,24] ## infectious

plot(prev_infectious_0 ~ t,ylim = c(0,0.3),type="l",col="darkred",xlim=c(1800,3000))
lines(prev_infectious_1 ~ t, lty=2,col="darkred")
lines(prev_infectious_2 ~ t, lty=3,col="darkred")

## y_delta0
age_0_3 = y[,2]
age_4_6 = y[,3]+y[,9]

prop_age_under_6 = (age_0_3 + age_4_6)/y[,24]

age_7_9 = y[,4]+y[,10]
age_10_12 = y[,5]+y[,11]
age_13_15 = y[,6]+y[,12]+y[,15]
age_16_18 = y[,7]+y[,13]+y[,16]
age_19_21 = y[,8]+y[,14]+y[,17]

prop_age_over_13 = (age_13_15 + age_16_18 + age_19_21)/y_delta0[,24]

## y_delta1
age_0_3_d1 = y_delta1[,2]
age_4_6_d1 = y_delta1[,3]+y_delta1[,9]

prop_age_under_6_d1 = (age_0_3_d1 + age_4_6_d1)/y_delta1[,24]

age_7_9_d1 = y_delta1[,4]+y_delta1[,10]
age_10_12_d1 = y_delta1[,5]+y_delta1[,11]
age_13_15_d1 = y_delta1[,6]+y_delta1[,12]+y_delta1[,15]
age_16_18_d1 = y_delta1[,7]+y_delta1[,13]+y_delta1[,16]
age_19_21_d1 = y_delta1[,8]+y_delta1[,14]+y_delta1[,17]

prop_age_over_13_d1 = (age_13_15_d1 + age_16_18_d1 + age_19_21_d1)/y_delta1[,24]


## y_delta2
age_0_3_d2 = y_delta2[,2]
age_4_6_d2 = y_delta2[,3]+y_delta2[,9]

prop_age_under_6_d2 = (age_0_3_d2 + age_4_6_d2)/y_delta2[,24]

age_7_9_d2 = y_delta2[,4]+y_delta2[,10]
age_10_12_d2 = y_delta2[,5]+y_delta2[,11]
age_13_15_d2 = y_delta2[,6]+y_delta2[,12]+y_delta2[,15]
age_16_18_d2 = y_delta2[,7]+y_delta2[,13]+y_delta2[,16]
age_19_21_d2 = y_delta2[,8]+y_delta2[,14]+y_delta2[,17]

prop_age_over_13_d2 = (age_13_15_d2 + age_16_18_d2 + age_19_21_d2)/y_delta2[,24]

td = 2121
barplot(c(age_0_3[td],age_0_3_d1[td],age_0_3_d2[td],NA,
          age_4_6[td],age_4_6_d1[td],age_4_6_d2[td],NA,
          age_7_9[td],age_7_9_d1[td],age_7_9_d2[td],NA,
          age_10_12[td],age_10_12_d1[td],age_10_12_d2[td],NA,
          age_13_15[td],age_13_15_d1[td],age_13_15_d2[td],NA,
          age_16_18[td],age_16_18_d1[td],age_16_18_d2[td],NA,
          age_19_21[td],age_19_21_d1[td],age_19_21_d2[td]),
        xlab="Age group (days)",col = c(adegenet::transp("darkblue",c(1,0.65,0.25)),NA),
        ylim=c(0,400))
axis(1,at=seq(2,31,length=7),labels=c("0-3","4-6","7-9","10-12","13-15","16-19","20-21"))

legend("topright",
       legend = c("FOIv - 1/5",
                  "FOIv - 1/3",
                  "FOIv - 1/8"),
       col = c(adegenet::transp("darkblue",c(1,0.65,0.25)),NA),
       pch=15,bty="n")

td = 2160
barplot(c(age_0_3[td],age_0_3_d1[td],age_0_3_d2[td],NA,
          age_4_6[td],age_4_6_d1[td],age_4_6_d2[td],NA,
          age_7_9[td],age_7_9_d1[td],age_7_9_d2[td],NA,
          age_10_12[td],age_10_12_d1[td],age_10_12_d2[td],NA,
          age_13_15[td],age_13_15_d1[td],age_13_15_d2[td],NA,
          age_16_18[td],age_16_18_d1[td],age_16_18_d2[td],NA,
          age_19_21[td],age_19_21_d1[td],age_19_21_d2[td]),
        xlab="Age group (days)",col = c(adegenet::transp("darkblue",c(1,0.65,0.25)),NA),
        ylim=c(0,50))
axis(1,at=seq(2,31,length=7),labels=c("0-3","4-6","7-9","10-12","13-15","16-19","20-21"))

# plot(y[1800:3000,24] ~ t[1800:3000])
# lines(1000*prop_age_under_6 ~ t,col="blue")
# lines(1000*prop_age_over_13 ~ t,col="red")
