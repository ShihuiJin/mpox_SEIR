#SEIR models 
#leaky vaccine model
quarantine_sim=function(paras)
{
  #SEIR model----------------
  seir <- function(t, state, parameters){
    with(as.list(c(state, parameters)),{
      dS1 <- -(beta_c*(I1_u+I2_u+I1_d+I2_d+c_iso*(QI1+QI2+I1_iso+I2_iso+I1_q+I2_q+I1_s+I2_s))+
                 beta_s*(I1_u+I1_d+c_iso*(I1_iso+I1_q+I1_s+QI1)))*S1
      dS2 <- -beta_c*(I1_u+I2_u+I1_d+I2_d+c_iso*(QI1+QI2+I1_iso+I2_iso+I1_q+I2_q+I1_s+I2_s))*S2
      
      dE1 <- (beta_c*(I1_u+I2_u+I1_d+I2_d+c_iso*(QI1+QI2+I1_iso+I2_iso+I1_q+I2_q+I1_s+I2_s))+
                beta_s*(I1_u+I1_d+c_iso*(I1_iso+I1_q+I1_s+QI1)))*S1-alpha*E1 #
      dE2 <- beta_c*(I1_u+I2_u+I1_d+I2_d+c_iso*(QI1+QI2+I1_iso+I2_iso+I1_q+I2_q+I1_s+I2_s))*S2-alpha*E2 #
      
      #isolation due to contact tracing (assuming they are exposed)
      #phi_Q: inverse tracing duration
      dQ1 <- -phi_Q*Q1
      dQ2 <- -phi_Q*Q2
      #remaining incubation
      dQ11 <- phi_Q*Q1-alpha_q*Q11
      dQ12 <- phi_Q*Q2-alpha_q*Q12
      
      #infected after quarantine
      #alpha_q: alpha_q^(-1)=alpha^(-1)-phi_Q^(-1)
      dQI1 <- alpha_q*Q11*dr-phi_I*QI1
      dQI2 <- alpha_q*Q12*dr-phi_I*QI2
      dI1_q <- alpha_q*Q11*(1-dr)-gamma_I0*I1_q
      dI2_q <- alpha_q*Q12*(1-dr)-gamma_I0*I2_q
      
      #distinguishing infections (undiagnosed vs diagnosed)
      #to be diagnosed
      dI1_d <- alpha*E1*dr-phi_I*I1_d
      dI2_d <- alpha*E2*dr-phi_I*I2_d
      #undiagnosed: directly to recovery
      dI1_u <- alpha*E1*(1-dr)-gamma_I0*I1_u
      dI2_u <- alpha*E2*(1-dr)-gamma_I0*I2_u
      
      #from diagnosis to recovery/death
      #add two compartments to balance the hospitalization rate for quarantined ppl
      dI1_iso <- phi_I*(I1_d+QI1)*(1-pD)-gamma_II*I1_iso 
      dI2_iso <- phi_I*(I2_d+QI2)*(1-pD)-gamma_II*I2_iso
      dI1_s <- phi_I*(I1_d+QI1)*pD-gamma_D*I1_s 
      dI2_s <- phi_I*(I2_d+QI2)*pD-gamma_D*I2_s
      
      #recovered for ordinary infections
      dRI_q <- gamma_II*(I1_iso+I2_iso) #diagnosed
      dRI_u <- gamma_I0*(I1_u+I2_u) #undiagnosed
      dRI_qu <- gamma_I0*(I1_q+I2_q) #mild quarantined
      
      #death pathway
      dRD <- gamma_D*(I1_s+I2_s)
      
      list(c(dS1, dS2, dE1, dE2,
             dQ1, dQ2, dQ11, dQ12, dQI1, dQI2,
             dI1_d, dI2_d, dI1_u, dI2_u,
             dI1_iso, dI2_iso, dI1_q, dI2_q, dI1_s, dI2_s,
             dRI_q, dRI_u, dRI_qu,dRD))
    })
  }
  
  #parameters---------------
  {
    N=paras$N
    time=paras$time
    
    #population involved
    p_susceptible=paras$p_susceptible #people not immuned
    #p_active=paras$p_active #highly sexually active proportion among suscepted population
    p_susceptible_a=paras$p_susceptible_a # % not immuned among highly sexually active population
    p_susceptible_a0=paras$p_susceptible_a0 # % highly sexually active population
    #p_susceptible1=paras$p_susceptible1 # highly sexually active susceptible population
    #p_susceptible2=p_susceptible-p_susceptible1
    
    #transmissibility parameter
    incub=paras$incub #incubation period
    recovery_I=paras$recovery_I
    recovery_I0=paras$recovery_I0 #recovery period for asymptomatic
    recovery_D=paras$recovery_D #death pathway: icu-death
    
    #case-fatality rate
    pD=paras$pD
    
    #Reff
    R0_c=paras$R0_c*p_susceptible # community transmission
    R0_s=paras$R0_s*p_susceptible_a # sexual transmission
    
    #diagnosis rate
    dr=paras$dr
    
    #surveillance efforts
    ct_max=paras$ct_max #number of contact tracing
    ct_connect=paras$ct_connect #average number of contacts per infected
    notif=paras$notif #notification period (i.e., from diagnosis to isolation)
    c_iso=paras$c_iso  #reduced capability of infecting others
    t_quar=paras$t_quar #contact tracing duration
    
    #importation on each day
    imp=paras$imp
    
  }
  
  #preparing the model parameters------------
  {
    gamma_I0=1/recovery_I0
    gamma_II=1/(recovery_I-notif)
    gamma_D=1/recovery_D
    
    beta_c=R0_c/N/recovery_I
    beta_s=R0_s/(N*p_susceptible_a0)/recovery_I
    alpha=1/incub
    alpha_q=1/(incub-t_quar)
    
    phi_I=1/notif
    phi_Q=1/t_quar
  }
  
  #initial value, assuming population either susceptible or exposed--------------
  {
    init <-c(S1=N*p_susceptible_a0, S2=N*(1-p_susceptible_a0),
             E1=0, E2=0, Q1=0,Q2=0,Q11=0,Q12=0, QI1=0, QI2=0,
             I1_d=0, I2_d=0, I1_u=0, I2_u=0,
             I1_iso=0, I2_iso=0, I1_q=0, I2_q=0, I1_s=0, I2_s=0, 
             RI_q=0, RI_u=0, RI_qu=0, RD=0)
  }
  
  
  #modelling------------
  {
    parameters <- c(beta_s = beta_s, beta_c=beta_c, alpha=alpha,
                    gamma_I0=gamma_I0, gamma_II = gamma_II,
                    gamma_D=gamma_D,
                    phi_I=phi_I,phi_Q=phi_Q,
                    pD=pD, N=N)
    if(is.null(paras$high_init)|is.null(paras$low_init)){
      init['E1']=0
    }else{
      init['E1']=paras$high_init; init['S1']=init['S1']-paras$high_init
      init['E2']=paras$low_init; init['S2']=init['S2']-paras$low_init
    }
    #with importation
    init1=init
    out2=c()
    case_ct=c() #number of case to contact trace
    quarantined=0 #number of quarantined individuals at each time point
    for(t in 1:time)
    {
      #importation
      init1['E1']=init1['E1']+imp[t]*p_susceptible_a0; init1['S1']=init1['S1']-imp[t]*p_susceptible_a0
      init1['E2']=init1['E2']+imp[t]*(1-p_susceptible_a0); init1['S2']=init1['S2']-imp[t]*(1-p_susceptible_a0)
      
      #quarantine close contacts (local)
      if(t>1){
        q_a=min(min(wild_case,ct_max)*ct_connect*p_susceptible_a0,init1['E1'])
        q_l=min(min(wild_case,ct_max)*ct_connect*(1-p_susceptible_a0),init1['E2'])
        init1['E1']=init1['E1']-q_a; init1['Q1']=init1['Q1']+q_a
        init1['E2']=init1['E2']-q_l; init1['Q2']=init1['Q2']+q_l
        quarantined=c(quarantined,q_a+q_l)
      }
      prev_case=sum(init1[c('Q1','Q2','Q11','Q12','QI1','QI2','I1_iso','I2_iso',
                            'I1_q','I2_q','I1_s','I2_s','RI_q','RI_qu','RD')])
      out <- ode(y = init1, times = seq(0, 1), func = seir, parms = parameters)
      out1=data.frame(out)%>%filter(time==1)%>%select(-time)%>%as.matrix()
      names(out1)=names(init)
      init1=out1
      out2=rbind(out2,out1)
      #derive number of new, unquarantined local infection
      curr_case=sum(out1[c('Q1','Q2','Q11','Q12','QI1','QI2','I1_iso','I2_iso',
                           'I1_q','I2_q','I1_s','I2_s','RI_q','RI_qu','RD')])
      wild_case=max(curr_case-prev_case,0)
      case_ct=c(case_ct,wild_case)
    }
    out2=cbind(1:time, out2)%>%data.frame()
    names(out2)[1]='time'
    
  }
  
  #calculations----------------
  {
    #active
    active_i=rowSums(out2[,c('I1_d','I2_d','QI1','QI2',
                             'I1_u','I2_u','I1_s','I2_s',
                             'I1_iso','I2_iso','I1_q','I2_q')])
    active_q=rowSums(out2[,c('Q1','Q2','Q11','Q12')]) #active quarantined but not infected
    active_qi=rowSums(out2[,c('QI1','QI2','I1_q','I2_q')]) #active quarantined infected (suspected cases)
    active_iso=rowSums(out2[,c('I1_iso','I2_iso','I1_s','I2_s')])#exclude those who are about to be notified
    
    #infections (including asymptomatic)
    cum_inf=rowSums(out2[,c('I1_d','I2_d','I1_u','I2_u','QI1','QI2',
                            'I1_iso','I2_iso','I1_q','I2_q','I1_s','I2_s',
                            'RI_q','RI_u','RI_qu','RD')])
    new_inf=cum_inf-c(0,cum_inf[-time])
    
    #reported cases (after notification)
    cum_case=rowSums(out2[,c('I1_iso','I2_iso','I1_s','I2_s','RI_q','RD')])
    new_case=cum_case-c(0,cum_case[-time])
    
    #diagnosed during quarantine
    cum_qi=cumsum(quarantined)-rowSums(out2[,c('Q1','Q2','Q11','Q12','QI1','QI2','I1_q','I2_q','RI_qu')])
    new_qi=cum_qi-c(0,cum_qi[-time])
    #infectious during quarantine (before notification)
    cum_qi0=cumsum(quarantined)-rowSums(out2[,c('Q1','Q2','Q11','Q12')])
    new_qi0=cum_qi0-c(0,cum_qi0[-time])
    
    #new isolated (including quarantined and later diagnosed)
    cum_iso_all=rowSums(out2[,c('I1_iso','I2_iso','I1_s','I2_s','RI_q','RD')])
    new_iso_all=cum_iso_all-c(0,cum_iso_all[-time])
    
    #cumulative deaths
    cum_d=out2$RD
  }
  
  #output------------
  {
    #peak infection, peak time, prevalence (local), peak active quarantine/isolation/infection/hospitalization/icu, deaths
    summary_stat=c(peak_case=max(new_case), peak_inf=max(new_inf),
                   cum_inf=cum_inf[time], prev=cum_inf[time]/N*100,
                   #prev_affected=cum_inf[time]/N/p_susceptible*100,
                   t0_cum=which.max(cum_case>=10), #localized outbreak starting point
                   t0_new=which.max(new_case>=10), #localized outbreak starting point
                   # cum_qso_pd=sum(active_q+active_qi+active_iso),
                   # max_quarantie=max(active_q+active_qi), max_iso=max(active_iso),
                   # max_qso=max(active_q+active_qi+active_iso),
                   max_inf=max(active_i), all_deaths=cum_d[time], 
                   max_inf_t = which(max(new_inf)==new_inf)/30) #in month
    full_stat=data.frame(time=1:time, cum_inf=cum_inf,
                         new_case=new_case, cum_case=cum_case,
                         cum_prev=cum_inf/N*100,
                         # new_case_ct=case_ct, active_quarantine=active_q,
                         # active_q=active_q+active_qi,active_iso=active_iso,
                         # active_qso=active_q+active_qi+active_iso,
                         active_i=active_i, cum_deaths=cum_d)
    list(summary_stat=summary_stat,full_stat=full_stat)
  }
}

#all-or-nothing vaccine model
quarantine_sim_a=function(paras)
{
  #SEIR model----------------
  seir <- function(t, state, parameters){
    with(as.list(c(state, parameters)),{
      dS1 <- -(beta_c*(I1_u+I2_u+I1_d+I2_d+c_iso*(QI1+QI2+I1_iso+I2_iso+I1_q+I2_q+I1_s+I2_s))+
                 beta_s*(I1_u+I1_d+c_iso*(I1_iso+I1_q+I1_s+QI1)))*S1
      dS2 <- -beta_c*(I1_u+I2_u+I1_d+I2_d+c_iso*(QI1+QI2+I1_iso+I2_iso+I1_q+I2_q+I1_s+I2_s))*S2
      
      dE1 <- (beta_c*(I1_u+I2_u+I1_d+I2_d+c_iso*(QI1+QI2+I1_iso+I2_iso+I1_q+I2_q+I1_s+I2_s))+
                beta_s*(I1_u+I1_d+c_iso*(I1_iso+I1_q+I1_s+QI1)))*S1-alpha*E1 #
      dE2 <- beta_c*(I1_u+I2_u+I1_d+I2_d+c_iso*(QI1+QI2+I1_iso+I2_iso+I1_q+I2_q+I1_s+I2_s))*S2-alpha*E2 #
      
      #isolation due to contact tracing (assuming they are exposed)
      dQ1 <- -phi_Q*Q1
      dQ2 <- -phi_Q*Q2
      #remaining incubation
      dQ11 <- phi_Q*Q1-alpha_q*Q11
      dQ12 <- phi_Q*Q2-alpha_q*Q12
      
      #infected after quarantine
      dQI1 <- alpha_q*Q11*dr-phi_I*QI1
      dQI2 <- alpha_q*Q12*dr-phi_I*QI2
      dI1_q <- alpha_q*Q11*(1-dr)-gamma_I0*I1_q
      dI2_q <- alpha_q*Q12*(1-dr)-gamma_I0*I2_q
      
      #distinguishing infections (undiagnosed vs diagnosed)
      #to be diagnosed
      dI1_d <- alpha*E1*dr-phi_I*I1_d
      dI2_d <- alpha*E2*dr-phi_I*I2_d
      #undiagnosed: directly to recovery
      dI1_u <- alpha*E1*(1-dr)-gamma_I0*I1_u
      dI2_u <- alpha*E2*(1-dr)-gamma_I0*I2_u
      
      #from diagnosis to recovery/death
      #add two compartments to balance the hospitalization rate for quarantined ppl
      dI1_iso <- phi_I*(I1_d+QI1)*(1-pD)-gamma_II*I1_iso 
      dI2_iso <- phi_I*(I2_d+QI2)*(1-pD)-gamma_II*I2_iso
      dI1_s <- phi_I*(I1_d+QI1)*pD-gamma_D*I1_s 
      dI2_s <- phi_I*(I2_d+QI2)*pD-gamma_D*I2_s
      
      #recovered for ordinary infections
      dRI_q <- gamma_II*(I1_iso+I2_iso) #diagnosed
      dRI_u <- gamma_I0*(I1_u+I2_u) #undiagnosed
      dRI_qu <- gamma_I0*(I1_q+I2_q) #mild quarantined
      
      #death pathway
      dRD <- gamma_D*(I1_s+I2_s)
      
      list(c(dS1, dS2, dE1, dE2,
             dQ1, dQ2, dQ11, dQ12, dQI1, dQI2,
             dI1_d, dI2_d, dI1_u, dI2_u,
             dI1_iso, dI2_iso, dI1_q, dI2_q, dI1_s, dI2_s,
             dRI_q, dRI_u, dRI_qu,dRD))
    })
  }
  
  #parameters---------------
  {
    N=paras$N
    time=paras$time
    
    #population involved
    p_susceptible=paras$p_susceptible #people not immuned
    p_active=paras$p_active #highly sexually active proportion among suscepted population
    p_susceptible_a=paras$p_susceptible_a # % not immuned among highly sexually active population
    p_susceptible_a0=paras$p_susceptible_a0 # % highly sexually active population
    p_susceptible1=paras$p_susceptible1 # highly sexually active susceptible population
    p_susceptible2=p_susceptible-p_susceptible1
    
    #transmissibility parameter
    incub=paras$incub #incubation period
    recovery_I=paras$recovery_I
    recovery_I0=paras$recovery_I0
    recovery_D=paras$recovery_D #death pathway: icu-death
    
    #case-fatality rate
    pD=paras$pD
    
    #Reff
    R0_c=paras$R0_c # community transmission
    R0_s=paras$R0_s # sexual transmission
    
    #diagnosis rate
    dr=paras$dr
    
    #surveillance efforts
    ct_max=paras$ct_max #number of contact tracing
    ct_connect=paras$ct_connect #average number of contacts per infected
    notif=paras$notif #notification period (i.e., from diagnosis to isolation)
    c_iso=paras$c_iso  #reduced capability of infecting others
    t_quar=paras$t_quar #contact tracing duration
    
    #importation on each day
    imp=paras$imp
    
  }
  
  #preparing the model parameters------------
  {
    gamma_I0=1/recovery_I0
    gamma_II=1/(recovery_I-notif)
    gamma_D=1/recovery_D
    
    beta_c=R0_c/N/recovery_I
    beta_s=R0_s/(N*p_susceptible_a0)/recovery_I
    alpha=1/incub
    alpha_q=1/(incub-t_quar)
    
    phi_I=1/notif
    phi_Q=1/t_quar
  }
  
  #initial value, assuming population either susceptible or exposed--------------
  {
    init <-c(S1=N*p_susceptible1, S2=N*p_susceptible2,
             E1=0, E2=0, Q1=0,Q2=0,Q11=0,Q12=0, QI1=0, QI2=0,
             I1_d=0, I2_d=0, I1_u=0, I2_u=0,
             I1_iso=0, I2_iso=0, I1_q=0, I2_q=0, I1_s=0, I2_s=0, 
             RI_q=0, RI_u=0, RI_qu=0, RD=0)
  }
  
  
  #modelling------------
  {
    parameters <- c(beta_s = beta_s, beta_c=beta_c, alpha=alpha,
                    gamma_I0=gamma_I0, gamma_II = gamma_II,
                    gamma_D=gamma_D,
                    phi_I=phi_I,phi_Q=phi_Q,
                    pD=pD, N=N)
    if(is.null(paras$high_init)|is.null(paras$low_init)){
      init['E1']=0
    }else{
      init['E1']=paras$high_init; init['S1']=init['S1']-paras$high_init
      init['E2']=paras$low_init; init['S2']=init['S2']-paras$low_init
    }
    #with importation
    init1=init
    out2=c()
    case_ct=c() #number of case to contact trace
    quarantined=0 #number of quarantined individuals at each time point
    for(t in 1:time)
    {
      #importation
      init1['E1']=init1['E1']+imp[t]*p_active; init1['S1']=init1['S1']-imp[t]*p_active
      init1['E2']=init1['E2']+imp[t]*(1-p_active); init1['S2']=init1['S2']-imp[t]*(1-p_active)
      
      #quarantine close contacts (local)
      if(t>1){
        q_a=min(min(wild_case,ct_max)*ct_connect*p_susceptible_a0,init1['E1'])
        q_l=min(min(wild_case,ct_max)*ct_connect*(1-p_susceptible_a0),init1['E2'])
        init1['E1']=init1['E1']-q_a; init1['Q1']=init1['Q1']+q_a
        init1['E2']=init1['E2']-q_l; init1['Q2']=init1['Q2']+q_l
        quarantined=c(quarantined,q_a+q_l)
      }
      prev_case=sum(init1[c('Q1','Q2','Q11','Q12','QI1','QI2','I1_iso','I2_iso',
                            'I1_q','I2_q','I1_s','I2_s','RI_q','RI_qu','RD')])
      out <- ode(y = init1, times = seq(0, 1), func = seir, parms = parameters)
      out1=data.frame(out)%>%filter(time==1)%>%select(-time)%>%as.matrix()
      names(out1)=names(init)
      init1=out1
      out2=rbind(out2,out1)
      #derive number of new, unquarantined local infection
      curr_case=sum(out1[c('Q1','Q2','Q11','Q12','QI1','QI2','I1_iso','I2_iso',
                           'I1_q','I2_q','I1_s','I2_s','RI_q','RI_qu','RD')])
      wild_case=max(curr_case-prev_case,0)
      case_ct=c(case_ct,wild_case)
    }
    out2=cbind(1:time, out2)%>%data.frame()
    names(out2)[1]='time'
    
  }
  
  #calculations----------------
  {
    #active infectious
    active_i=rowSums(out2[,c('I1_d','I2_d','QI1','QI2',
                             'I1_u','I2_u','I1_s','I2_s',
                             'I1_iso','I2_iso','I1_q','I2_q')])
    
    #infections (including asymptomatic)
    cum_inf=rowSums(out2[,c('I1_d','I2_d','I1_u','I2_u','QI1','QI2',
                            'I1_iso','I2_iso','I1_q','I2_q','I1_s','I2_s',
                            'RI_q','RI_u','RI_qu','RD')])
    new_inf=cum_inf-c(0,cum_inf[-time])
    
    #reported cases (after notification)
    cum_case=rowSums(out2[,c('I1_iso','I2_iso','I1_s','I2_s','RI_q','RD')])
    new_case=cum_case-c(0,cum_case[-time])
    
    #cumulative deaths
    cum_d=out2$RD
  }
  
  #output------------
  {
    #peak infection, peak time, prevalence (local), peak active quarantine/isolation/infection/hospitalization/icu, deaths
    summary_stat=c(peak_case=max(new_case), peak_inf=max(new_inf),
                   cum_inf=cum_inf[time], prev=cum_inf[time]/N*100,
                   t0_cum=which.max(cum_case>=10), #localized outbreak starting point
                   t0_new=which.max(new_case>=10), #localized outbreak starting point
                   max_inf=max(active_i), all_deaths=cum_d[time], 
                   max_inf_t = which(max(new_inf)==new_inf)/30) #in month
    full_stat=data.frame(time=1:time, cum_inf=cum_inf,
                         new_case=new_case, cum_case=cum_case,
                         cum_prev=cum_inf/N*100,
                         active_i=active_i, cum_deaths=cum_d)
    list(summary_stat=summary_stat,full_stat=full_stat)
  }
}