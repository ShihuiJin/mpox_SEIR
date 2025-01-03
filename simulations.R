#main analysis script
library(dplyr)
library(stringr)
library(parallel)
library(deSolve)
n.cores=detectCores()

#SEIR models
source('models.R')
#simualtion functions
source('sim_funcs.R')

#city-specific information
city_info=read.csv('city_info.csv')
age_dist=read.csv('national_age_dist.csv') #age distribution
immunity=read.csv("immunity.csv")


Time=360*5
paras=list(
  N=5e6,
  time=Time,
  
  #population involved
  p_susceptible=0.9, #people not immuned
  p_susceptible1=0.027, #people aged between 20-29 and not immuned
  p_susceptible_a=0.9, #people not immuned among high risk
  p_susceptible_a0=0.03, #people high risk
  p_active=.025,  #highly sexually active proportion among susceptible
  
  #transmissibility parameter
  incub=7, #incubation period
  recovery_I=21,
  recovery_I0=21,
  recovery_D=15, #death pathway: icu->death
  
  #case-hospitalization/icu/death rate
  pD=.03,
  
  #R0
  R0_c=1.021, #community transmission
  R0_s=1.617, #sextual transmission
  
  #diagnosis rate
  dr=0.25,
  
  #surveillance efforts
  ct_max=0, #number of contact tracing
  ct_connect=1, #average number of contacts per infected
  notif=6, #notification period (i.e., from diagnosis to isolation)
  c_iso=1,  #reduced capability of infecting others
  t_quar=3, #contact tracing duration
  
  #importation on each day
  imp=rep(0, Time),
  
  #initial exposure
  high_init=1,
  low_init=0
)



#analysis code------------
sens=F; risky=NA
#sens=T; risky=c(0.02,0.05,0.1,'low','high')[1]
#'low' stands for the lower-risk scenario where undiagnosed infections would recover in 14 days
#'high' stands for the high-risk scenario where intial importations are all high-risk individuals

#one initial exposure scenario-----------
n.imp=1
Time=360*5
out_imp1=list()
for(k in 0:n.imp){
  out_full=list()
  paras$c_iso=1
  paras$ct_max=0
  paras$high_init=k;paras$low_init=n.imp-k
  paras$imp=rep(0,Time)
  paras$recovery_I0=ifelse(risky%in%'low',14,21)
  for(c in city_info$ISOALPHA){
    if(sens&is.numeric(risky)){
      res=mpox_sim(c, Time=Time, sens=T, risky=risky)
    }else{
      res=mpox_sim(c, Time=Time)
    }
    out_full[[c]]=res
  }
  out_imp1[[k+1]]=out_full
}


#varying initial exposure sizes------------
out_imp=list()
for(n.imp in 1:3*2-1)
{
  {
    print(paste0(n.imp,': ',Sys.time()))
    out_full=list()
    paras$c_iso=1
    paras$ct_max=0
    paras$imp=c(n.imp,rep(0,Time-1))
    paras$high_init=paras$low_init=0
    for(c in city_info$ISOALPHA){
      if(sens&is.numeric(risky)){
        res=mpox_sim(c, Time=Time, sens=T, risky=risky)
      }else if(sens&(risky=='high')){
        paras$imp=rep(0,Time)
        paras$high_init=n.imp; paras$low_init=0
        res=mpox_sim(c, Time=Time)
      }else{
        res=mpox_sim(c, Time=Time)
      }
      out_full[[c]]=res
    }
  }
  out_imp[[(n.imp+1)/2]]=out_full
}

out_imp_summ0=out_sum(out_imp, 12, T)


#intervention strategies------------
n.imp=3; paras$imp=c(n.imp,rep(0,Time-1))
paras$high_init=paras$low_init=0
iso=1-rep(0:4/4,2)
quar=c(rep(0,5),rep(1,5))
out_control=list()
for(i in seq_along(iso)){
  print(paste0(i,': ',Sys.time()))
    paras$recovery_I0=ifelse(risky%in%'low',14,21)
    paras$c_iso=iso[i]
    paras$ct_max=ifelse(quar[i]==1,30,0) #[edit]
    out_full=list()
    for(c in city_info$ISOALPHA){
      if(sens&is.numeric(risky)){
        res=mpox_sim(c, Time=Time, sens=T, risky=risky)
      }else if(sens&(risky=='high')){
        paras$imp=rep(0,Time)
        paras$high_init=n.imp; paras$low_init=0
        res=mpox_sim(c, Time=Time)
      }else{
        res=mpox_sim(c, Time=Time)
      }
      out_full[[c]]=res
    }
    out_control[[i]]=out_full
}

out_control_summ0=out_sum(out_control, 12, F)
out_control_comp=out_sum_compare(out_control_summ0[[2]])




