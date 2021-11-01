#############################
##
##  What sample size do we need to reflect age structure?
##  How quickly do we need to collect mosquitoes for a snap shot? 
##
##


## Using model in "Aging_mosq_model1.R"

age_0_3 = y[,2]
age_4_6 = y[,3]+y[,9]
age_7_9 = y[,4]+y[,10]
age_10_12 = y[,5]+y[,11]
age_13_15 = y[,6]+y[,12]+y[,15]
age_16_18 = y[,7]+y[,13]+y[,16]
age_19_21 = y[,8]+y[,14]+y[,17]

par(mfrow=c(1,1))
## 1 day
td = 2000
tw = 2000:2006
tm = 2000:2030

totd = age_0_3[td] + age_4_6[td] + age_7_9[td] + age_10_12[td] + 
  age_13_15[td] + age_16_18[td] + age_19_21[td]

totw = sum(age_0_3[tw],age_4_6[tw],age_7_9[tw],age_10_12[tw],
  age_13_15[tw],age_16_18[tw],age_19_21[tw])

totm = sum(age_0_3[tm],age_4_6[tm],age_7_9[tm],age_10_12[tm],
           age_13_15[tm],age_16_18[tm],age_19_21[tm])

barplot(c(age_0_3[td]/totd,   sum(age_0_3[tw])/totw,   sum(age_0_3[tm])/totm,NA,
          age_4_6[td]/totd,   sum(age_4_6[tw])/totw,   sum(age_4_6[tm])/totm,NA,
          age_7_9[td]/totd,   sum(age_7_9[tw])/totw,   sum(age_7_9[tm])/totm,NA,
          age_10_12[td]/totd, sum(age_10_12[tw])/totw, sum(age_10_12[tm])/totm,NA,
          age_13_15[td]/totd, sum(age_13_15[tw])/totw, sum(age_13_15[tm])/totm,NA,
          age_16_18[td]/totd, sum(age_16_18[tw])/totw, sum(age_16_18[tm])/totm,NA,
          age_19_21[td]/totd, sum(age_19_21[tw])/totw, sum(age_19_21[tm])/totm),
        xlab="Age group (days)",col = c(adegenet::transp("darkblue",c(1,0.6,0.2)),NA))
axis(1,at=seq(1.5,31,length=7),labels=c("0-3","4-6","7-9","10-12","13-15","16-19","20-21"))




td = 2160
tw = 2160:2166
tm = 2160:2190

totd = age_0_3[td] + age_4_6[td] + age_7_9[td] + age_10_12[td] + 
  age_13_15[td] + age_16_18[td] + age_19_21[td]

totw = sum(age_0_3[tw],age_4_6[tw],age_7_9[tw],age_10_12[tw],
           age_13_15[tw],age_16_18[tw],age_19_21[tw])

totm = sum(age_0_3[tm],age_4_6[tm],age_7_9[tm],age_10_12[tm],
           age_13_15[tm],age_16_18[tm],age_19_21[tm])

barplot(c(age_0_3[td]/totd,   sum(age_0_3[tw])/totw,   sum(age_0_3[tm])/totm,NA,
          age_4_6[td]/totd,   sum(age_4_6[tw])/totw,   sum(age_4_6[tm])/totm,NA,
          age_7_9[td]/totd,   sum(age_7_9[tw])/totw,   sum(age_7_9[tm])/totm,NA,
          age_10_12[td]/totd, sum(age_10_12[tw])/totw, sum(age_10_12[tm])/totm,NA,
          age_13_15[td]/totd, sum(age_13_15[tw])/totw, sum(age_13_15[tm])/totm,NA,
          age_16_18[td]/totd, sum(age_16_18[tw])/totw, sum(age_16_18[tm])/totm,NA,
          age_19_21[td]/totd, sum(age_19_21[tw])/totw, sum(age_19_21[tm])/totm),
        xlab="Age group (days)",col = c(adegenet::transp("darkblue",c(1,0.6,0.2)),NA))
axis(1,at=seq(1.5,31,length=7),labels=c("0-3","4-6","7-9","10-12","13-15","16-19","20-21"))



############################################################
##
## Given the above, how big relative to all population
## does our sample size need to be? 

## Let's start with a weekly sampling strategy, in 
## the season when mosquitoes start to increase in number

## days 2300-2306

samp_all = sum(round(y[2300:2306,24],0))
age_0_3_samp_all = sum(round(y[2300:2306,2],0))
age_4_6_samp_all = sum(round(y[2300:2306,3]+y[2300:2306,9],0))
age_7_9_samp_all = sum(round(y[2300:2306,4]+y[2300:2306,10],0))
age_10_12_samp_all = sum(round(y[2300:2306,5]+y[2300:2306,11],0))
age_13_15_samp_all = sum(round(y[2300:2306,6]+y[2300:2306,12]+y[2300:2306,15],0))
age_16_18_samp_all = sum(round(y[2300:2306,7]+y[2300:2306,13]+y[2300:2306,16],0))
age_19_21_samp_all = sum(round(y[2300:2306,8]+y[2300:2306,14]+y[2300:2306,17],0))


## Simulated population
sim_pop = c(sample(0:3,replace=TRUE,size=age_0_3_samp_all),
            sample(4:6,replace=TRUE,size=age_4_6_samp_all),
            sample(7:9,replace=TRUE,size=age_7_9_samp_all),
            sample(10:12,replace=TRUE,size=age_10_12_samp_all),
            sample(13:15,replace=TRUE,size=age_13_15_samp_all),
            sample(16:18,replace=TRUE,size=age_16_18_samp_all),
            sample(19:21,replace=TRUE,size=age_19_21_samp_all))
length(sim_pop)

hist(sim_pop)



store = array(dim=c(4,4,7))
sampA = 50  ## 2%
sampB = 100 ## 4%
sampC = 250 ## 10%
sampD = 500 ## 20%
for(i in 1:10){
  ## Assuming entirely random distribution of the population
  test1 = sample(sim_pop, replace=FALSE,size = sampA) ## approx 2%
  test2 = sample(sim_pop, replace=FALSE,size = sampB) ## approx 4%
  test3 = sample(sim_pop, replace=FALSE,size = sampC) ## approx 10%
  test4 = sample(sim_pop, replace=FALSE,size = sampD) ## approx 40%
  
  age_0_3_test1 = length(test1[test1 < 4])/sampA
  age_4_6_test1 = length(test1[test1 > 3 & test1 < 7])/sampA
  age_7_9_test1 = length(test1[test1 > 6 & test1 < 10])/sampA
  age_10_12_test1 = length(test1[test1 > 9 & test1 < 13])/sampA
  age_13_15_test1 = length(test1[test1 > 12 & test1 < 16])/sampA
  age_16_18_test1 = length(test1[test1 > 15 & test1 < 19])/sampA
  age_19_21_test1 = length(test1[test1 > 18])/sampA
  
  
  age_0_3_test2 = length(test2[test2 < 4])/sampB
  age_4_6_test2 = length(test2[test2 > 3 & test2 < 7])/sampB
  age_7_9_test2 = length(test2[test2 > 6 & test2 < 10])/sampB
  age_10_12_test2 = length(test2[test2 > 9 & test2 < 13])/sampB
  age_13_15_test2 = length(test2[test2 > 12 & test2 < 16])/sampB
  age_16_18_test2 = length(test2[test2 > 15 & test2 < 19])/sampB
  age_19_21_test2 = length(test2[test2 > 18])/sampB
  
  
  age_0_3_test3 = length(test3[test3 < 4])/sampC
  age_4_6_test3 = length(test3[test3 > 3 & test3 < 7])/sampC
  age_7_9_test3 = length(test3[test3 > 6 & test3 < 10])/sampC
  age_10_12_test3 = length(test3[test3 > 9 & test3 < 13])/sampC
  age_13_15_test3 = length(test3[test3 > 12 & test3 < 16])/sampC
  age_16_18_test3 = length(test3[test3 > 15 & test3 < 19])/sampC
  age_19_21_test3 = length(test3[test3 > 18])/sampC
  
  
  age_0_3_test4 = length(test4[test4 < 4])/sampD
  age_4_6_test4 = length(test4[test4 > 3 & test4 < 7])/sampD
  age_7_9_test4 = length(test4[test4 > 6 & test4 < 10])/sampD
  age_10_12_test4 = length(test4[test4 > 9 & test4 < 13])/sampD
  age_13_15_test4 = length(test4[test4 > 12 & test4 < 16])/sampD
  age_16_18_test4 = length(test4[test4 > 15 & test4 < 19])/sampD
  age_19_21_test4 = length(test4[test4 > 18])/sampD
  
  store[,i,1] = c(age_0_3_test1,age_0_3_test2,age_0_3_test3,age_0_3_test4)
  store[,i,2] = c(age_4_6_test1,age_4_6_test2,age_4_6_test3,age_4_6_test4)
  store[,i,3] = c(age_7_9_test1,age_7_9_test2,age_7_9_test3,age_7_9_test4)
  store[,i,4] = c(age_10_12_test1,age_10_12_test2,age_10_12_test3,age_10_12_test4)
  store[,i,5] = c(age_13_15_test1,age_13_15_test2,age_13_15_test3,age_13_15_test4)
  store[,i,6] = c(age_16_18_test1,age_16_18_test2,age_16_18_test3,age_16_18_test4)
  store[,i,7] = c(age_19_21_test1,age_19_21_test2,age_19_21_test3,age_19_21_test4)
  
}





TRUE_proportionage_0_3 = age_0_3_samp_all/samp_all
TRUE_proportionage_4_6 = age_4_6_samp_all/samp_all
TRUE_proportionage_7_9 = age_7_9_samp_all/samp_all
TRUE_proportionage_10_12 = age_10_12_samp_all/samp_all
TRUE_proportionage_13_15 = age_13_15_samp_all/samp_all
TRUE_proportionage_16_18 = age_16_18_samp_all/samp_all
TRUE_proportionage_19_21 = age_19_21_samp_all/samp_all


age_0_3_test1 = mean(store[1,,1])
age_0_3_test2 = mean(store[2,,1])
age_0_3_test3 = mean(store[3,,1])
age_0_3_test4 = mean(store[4,,1])

age_4_6_test1 = mean(store[1,,2])
age_4_6_test2 = mean(store[2,,2])
age_4_6_test3 = mean(store[3,,2])
age_4_6_test4 = mean(store[4,,2])

age_7_9_test1 = mean(store[1,,3])
age_7_9_test2 = mean(store[2,,3])
age_7_9_test3 = mean(store[3,,3])
age_7_9_test4 = mean(store[4,,3])

age_10_12_test1 = mean(store[1,,4])
age_10_12_test2 = mean(store[2,,4])
age_10_12_test3 = mean(store[3,,4])
age_10_12_test4 = mean(store[4,,4])

age_13_15_test1 = mean(store[1,,5])
age_13_15_test2 = mean(store[2,,5])
age_13_15_test3 = mean(store[3,,5])
age_13_15_test4 = mean(store[4,,5])

age_16_18_test1 = mean(store[1,,6])
age_16_18_test2 = mean(store[2,,6])
age_16_18_test3 = mean(store[3,,6])
age_16_18_test4 = mean(store[4,,6])

age_19_21_test1 = mean(store[1,,7])
age_19_21_test2 = mean(store[2,,7])
age_19_21_test3 = mean(store[3,,7])
age_19_21_test4 = mean(store[4,,7])

barplot(c(TRUE_proportionage_0_3,age_0_3_test1,age_0_3_test2,age_0_3_test3,age_0_3_test4,NA,
          TRUE_proportionage_4_6,age_4_6_test1,age_4_6_test2,age_4_6_test3,age_4_6_test4,NA,
          TRUE_proportionage_7_9,age_7_9_test1,age_7_9_test2,age_7_9_test3,age_7_9_test4,NA,
          TRUE_proportionage_10_12,age_10_12_test1,age_10_12_test2,age_10_12_test3,age_10_12_test4,NA,
          TRUE_proportionage_13_15,age_13_15_test1,age_13_15_test2,age_13_15_test3,age_13_15_test4,NA,
          TRUE_proportionage_16_18,age_16_18_test1,age_16_18_test2,age_16_18_test3,age_16_18_test4,NA,
          TRUE_proportionage_19_21,age_19_21_test1,age_19_21_test2,age_19_21_test3,age_19_21_test4),
        col = c("orange",adegenet::transp("darkblue",c(0.2,0.5,0.7,1)),NA),
        ylim=c(0,0.5),xaxt="n",xlab = "Age category (days)")

axis(1,at=seq(3,46,length=7),labels=c("0-3","4-6","7-9","10-12","13-15","16-19","20-21"))



age_0_3_test1u = as.numeric(quantile(store[1,,1],0.975))
age_0_3_test2u = as.numeric(quantile(store[2,,1],0.975))
age_0_3_test3u = as.numeric(quantile(store[3,,1],0.975))
age_0_3_test4u = as.numeric(quantile(store[4,,1],0.975))

age_4_6_test1u = as.numeric(quantile(store[1,,2],0.975))
age_4_6_test2u = as.numeric(quantile(store[2,,2],0.975))
age_4_6_test3u = as.numeric(quantile(store[3,,2],0.975))
age_4_6_test4u = as.numeric(quantile(store[4,,2],0.975))

age_7_9_test1u = as.numeric(quantile(store[1,,3],0.975))
age_7_9_test2u = as.numeric(quantile(store[2,,3],0.975))
age_7_9_test3u = as.numeric(quantile(store[3,,3],0.975))
age_7_9_test4u = as.numeric(quantile(store[4,,3],0.975))

age_10_12_test1u = as.numeric(quantile(store[1,,4],0.975))
age_10_12_test2u = as.numeric(quantile(store[2,,4],0.975))
age_10_12_test3u = as.numeric(quantile(store[3,,4],0.975))
age_10_12_test4u = as.numeric(quantile(store[4,,4],0.975))

age_13_15_test1u = as.numeric(quantile(store[1,,5],0.975))
age_13_15_test2u = as.numeric(quantile(store[2,,5],0.975))
age_13_15_test3u = as.numeric(quantile(store[3,,5],0.975))
age_13_15_test4u = as.numeric(quantile(store[4,,5],0.975))

age_16_18_test1u = as.numeric(quantile(store[1,,6],0.975))
age_16_18_test2u = as.numeric(quantile(store[2,,6],0.975))
age_16_18_test3u = as.numeric(quantile(store[3,,6],0.975))
age_16_18_test4u = as.numeric(quantile(store[4,,6],0.975))

age_19_21_test1u = as.numeric(quantile(store[1,,7],0.975))
age_19_21_test2u = as.numeric(quantile(store[2,,7],0.975))
age_19_21_test3u = as.numeric(quantile(store[3,,7],0.975))
age_19_21_test4u = as.numeric(quantile(store[4,,7],0.975))



age_0_3_test1l = as.numeric(quantile(store[1,,1],0.025))
age_0_3_test2l = as.numeric(quantile(store[2,,1],0.025))
age_0_3_test3l = as.numeric(quantile(store[3,,1],0.025))
age_0_3_test4l = as.numeric(quantile(store[4,,1],0.025))

age_4_6_test1l = as.numeric(quantile(store[1,,2],0.025))
age_4_6_test2l = as.numeric(quantile(store[2,,2],0.025))
age_4_6_test3l = as.numeric(quantile(store[3,,2],0.025))
age_4_6_test4l = as.numeric(quantile(store[4,,2],0.025))

age_7_9_test1l = as.numeric(quantile(store[1,,3],0.025))
age_7_9_test2l = as.numeric(quantile(store[2,,3],0.025))
age_7_9_test3l = as.numeric(quantile(store[3,,3],0.025))
age_7_9_test4l = as.numeric(quantile(store[4,,3],0.025))

age_10_12_test1l = as.numeric(quantile(store[1,,4],0.025))
age_10_12_test2l = as.numeric(quantile(store[2,,4],0.025))
age_10_12_test3l = as.numeric(quantile(store[3,,4],0.025))
age_10_12_test4l = as.numeric(quantile(store[4,,4],0.025))

age_13_15_test1l = as.numeric(quantile(store[1,,5],0.025))
age_13_15_test2l = as.numeric(quantile(store[2,,5],0.025))
age_13_15_test3l = as.numeric(quantile(store[3,,5],0.025))
age_13_15_test4l = as.numeric(quantile(store[4,,5],0.025))

age_16_18_test1l = as.numeric(quantile(store[1,,6],0.025))
age_16_18_test2l = as.numeric(quantile(store[2,,6],0.025))
age_16_18_test3l = as.numeric(quantile(store[3,,6],0.025))
age_16_18_test4l = as.numeric(quantile(store[4,,6],0.025))

age_19_21_test1l = as.numeric(quantile(store[1,,7],0.025))
age_19_21_test2l = as.numeric(quantile(store[2,,7],0.025))
age_19_21_test3l = as.numeric(quantile(store[3,,7],0.025))
age_19_21_test4l = as.numeric(quantile(store[4,,7],0.025))

# axis(1,at=seq(0.8,48.9,length=41))

x0_seg = seq(0.8,48.9,length=41) 
x1_seg = seq(0.8,48.9,length=41) 

y0_seg =c(NA,age_0_3_test1l,age_0_3_test2l,age_0_3_test3l,age_0_3_test4l,
          NA,NA,age_4_6_test1l,age_4_6_test2l,age_4_6_test3l,age_4_6_test4l,
          NA,NA,age_7_9_test1l,age_7_9_test2l,age_7_9_test3l,age_7_9_test4l,
          NA,NA,age_10_12_test1l,age_10_12_test2l,age_10_12_test3l,age_10_12_test4l,
          NA,NA,age_13_15_test1l,age_13_15_test2l,age_13_15_test3l,age_13_15_test4l,
          NA,NA,age_16_18_test1l,age_16_18_test2l,age_16_18_test3l,age_16_18_test4l,
          NA,NA,age_19_21_test1l,age_19_21_test2l,age_19_21_test3l,age_19_21_test4l)

y1_seg =c(NA,age_0_3_test1u,age_0_3_test2u,age_0_3_test3u,age_0_3_test4u,
          NA,NA,age_4_6_test1u,age_4_6_test2u,age_4_6_test3u,age_4_6_test4u,
          NA,NA,age_7_9_test1u,age_7_9_test2u,age_7_9_test3u,age_7_9_test4u,
          NA,NA,age_10_12_test1u,age_10_12_test2u,age_10_12_test3u,age_10_12_test4u,
          NA,NA,age_13_15_test1u,age_13_15_test2u,age_13_15_test3u,age_13_15_test4u,
          NA,NA,age_16_18_test1u,age_16_18_test2u,age_16_18_test3u,age_16_18_test4u,
          NA,NA,age_19_21_test1u,age_19_21_test2u,age_19_21_test3u,age_19_21_test4u)

for(i in 1:41){
  segments(x0 = x0_seg[i], x1 = x1_seg[i],
           y0 = y0_seg[i], y1 = y1_seg[i])
}
