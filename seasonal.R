ssa0 <- 0.2852297					
ssa1 <- -0.2952712
ssa2 <- -0.03408224	
ssa3 = 0.07596435
ssb1 = -0.1126063
ssb2 = 0.07789561
ssb3 = -0.007051094


TIME = 1:365
data = (ssa0+ssa1*cos(2*pi*TIME/365)+ssa2*cos(2*2*pi*TIME/365)+
          ssa3*cos(3*2*pi*TIME/365)+ssb1*sin(2*pi*TIME/365)+
          ssb2*sin(2*2*pi*TIME/365)+ ssb3*sin(3*2*pi*TIME/365) )

data[which(data < 0)] = 0.001

theta_c =  mean(data)
