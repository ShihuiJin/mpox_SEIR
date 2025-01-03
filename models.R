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
      dQ1 <- -alpha_q*Q1
      dQ2 <- -alpha_q*Q2
      
      #infected after quarantine
      #alpha_q: alpha_q^(-1)=alpha^(-1)-phi_Q^(-1)
      dQI1 <- alpha_q*Q1*dr-phi_I*QI1
      dQI2 <- alpha_q*Q2*dr-phi_I*QI2
      dI1_q <- alpha_q*Q1*(1-dr)-gamma_I0*I1_q
      dI2_q <- alpha_q*Q2*(1-dr)-gamma_I0*I2_q
      
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
      dRI1_q <- gamma_II*I1_iso #diagnosed
      dRI1_u <- gamma_I0*I1_u #undiagnosed
      dRI1_qu <- gamma_I0*I1_q #mild quarantined
      dRI2_q <- gamma_II*I2_iso #diagnosed
      dRI2_u <- gamma_I0*I2_u #undiagnosed
      dRI2_qu <- gamma_I0*I2_q #mild quarantined
      
      #death pathway
      dRD1 <- gamma_D*I1_s
      dRD2 <- gamma_D*I2_s
      
      
      list(c(dS1, dS2, dE1, dE2,
             dQ1, dQ2, dQI1, dQI2,
             dI1_d, dI2_d, dI1_u, dI2_u,
             dI1_iso, dI2_iso, dI1_q, dI2_q, dI1_s, dI2_s,
             dRI1_q, dRI1_u, dRI1_qu,dRD1,dRI2_q, dRI2_u, dRI2_qu,dRD2))
    })
  }
  
  #parameters---------------
  {
    N=paras$N
    time=paras$time
    
    #population involved
    p_susceptible=paras$p_susceptible #people not immuned
    p_susceptible_a=paras$p_susceptible_a # % not immuned among highly sexually active population
    p_susceptible_a0=paras$p_susceptible_a0 # % highly sexually active population

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
    #phi_Q=1/t_quar
  }
  
  #initial value, assuming population either susceptible or exposed--------------
  {
    init <-c(S1=N*p_susceptible_a0, S2=N*(1-p_susceptible_a0),
             E1=0, E2=0, Q1=0,Q2=0, QI1=0, QI2=0,
             I1_d=0, I2_d=0, I1_u=0, I2_u=0,
             I1_iso=0, I2_iso=0, I1_q=0, I2_q=0, I1_s=0, I2_s=0, 
             RI1_q=0, RI1_u=0, RI1_qu=0, RD1=0,
             RI2_q=0, RI2_u=0, RI2_qu=0, RD2=0)
  }
  
  
  #modelling------------
  {
    parameters <- c(beta_s = beta_s, beta_c=beta_c, alpha=alpha,
                    gamma_I0=gamma_I0, gamma_II = gamma_II,
                    gamma_D=gamma_D,
                    phi_I=phi_I,#phi_Q=phi_Q,
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
    quarantined=rep(0,time) #number of quarantined individuals at each time point
    for(t in 1:time)
    {
      #importation
      init1['E1']=init1['E1']+imp[t]*p_susceptible_a0; init1['S1']=init1['S1']-imp[t]*p_susceptible_a0
      init1['E2']=init1['E2']+imp[t]*(1-p_susceptible_a0); init1['S2']=init1['S2']-imp[t]*(1-p_susceptible_a0)
      
      #quarantine close contacts (local)
      if(t>3 & ct_max>0){
        q_a0=sum(init1[c('E1','I1_d','I1_u')])
        q_l0=sum(init1[c('E2','I2_d','I2_u')])
        q0=min(sum(case_ct[t-3,]),ct_max)
        q_a1=min(q0*case_ct[t-3,1]/sum(case_ct[t-3,]),case_ct[t-3,1])
        q_l1=min(q0-q_a1,case_ct[t-3,2])
        q_a=min(q_a1*ct_connect,q_a0)
        q_l=min(q_l1*ct_connect,q_l0)
        init1['Q1']=init1['Q1']+q_a/q_a0*init1['E1']; init1['E1']=init1['E1']*(1-q_a/q_a0)
        init1['Q2']=init1['Q2']+q_l/q_l0*init1['E2']; init1['E2']=init1['E2']*(1-q_l/q_l0)
        init1['QI1']=init1['QI1']+q_a/q_a0*init1['I1_d']; init1['I1_d']=init1['I1_d']*(1-q_a/q_a0)
        init1['QI2']=init1['QI2']+q_l/q_l0*init1['I2_d']; init1['I2_d']=init1['I2_d']*(1-q_l/q_l0)
        init1['I1_q']=init1['I1_q']+q_a/q_a0*init1['I1_u']; init1['I1_u']=init1['I1_u']*(1-q_a/q_a0)
        init1['I2_q']=init1['I2_q']+q_l/q_l0*init1['I2_u']; init1['I2_u']=init1['I2_u']*(1-q_l/q_l0)
        quarantined[t]=q_a+q_l
      }
      prev_case=c(sum(init1[c('Q1','QI1','I1_iso','I1_q','I1_s','RI1_q','RI1_qu','RD1')]),
                  sum(init1[c('Q2','QI2','I2_iso','I2_q','I2_s','RI2_q','RI2_qu','RD2')]))
      out <- ode(y = init1, times = seq(0, 1), func = seir, parms = parameters)
      out1=data.frame(out)%>%filter(time==1)%>%select(-time)%>%as.matrix()
      names(out1)=names(init)
      init1=out1
      out2=rbind(out2,out1)
      #derive number of new, unquarantined local infection
      curr_case=c(sum(out1[c('Q1','QI1','I1_iso','I1_q','I1_s','RI1_q','RI1_qu','RD1')]),
                  sum(out1[c('Q2','QI2','I2_iso','I2_q','I2_s','RI2_q','RI2_qu','RD2')]))
      
      wild_case=lapply(1:2, function(i){
        max(curr_case[i]-prev_case[i],0)
      })%>%unlist()
      case_ct=rbind(case_ct,wild_case)
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
    
    #infections (including asymptomatic)
    cum_inf=rowSums(out2[,c('I1_d','I2_d','I1_u','I2_u','QI1','QI2',
                            'I1_iso','I2_iso','I1_q','I2_q','I1_s','I2_s',
                            'RI1_q','RI1_u','RI1_qu','RD1',
                            'RI2_q','RI2_u','RI2_qu','RD2')])
    new_inf=cum_inf-c(0,cum_inf[-time])
    
    #reported cases (after notification)
    cum_case=rowSums(out2[,c('I1_iso','I2_iso','I1_s','I2_s','RI1_q','RD1','RI2_q','RD2')])
    new_case=cum_case-c(0,cum_case[-time])
    
    #cumulative deaths
    cum_d=rowSums(out2[,c('RD1','RD2')])
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
    full_stat=data.frame(time=1:time, cum_inf=cum_inf,new_inf=new_inf,
                         new_case=new_case, cum_case=cum_case,
                         cum_prev=cum_inf/N*100,
                         
                         active_i=active_i, cum_deaths=cum_d)
    list(summary_stat=summary_stat,full_stat=full_stat)
  }
}

#all-or-nothing vaccine model (quarantine not considered)
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
      dQ1 <- -alpha_q*Q1
      dQ2 <- -alpha_q*Q2
      
      #infected after quarantine
      dQI1 <- alpha_q*Q1*dr-phi_I*QI1
      dQI2 <- alpha_q*Q2*dr-phi_I*QI2
      dI1_q <- alpha_q*Q1*(1-dr)-gamma_I0*I1_q
      dI2_q <- alpha_q*Q2*(1-dr)-gamma_I0*I2_q
      
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
      dRI1_q <- gamma_II*I1_iso #diagnosed
      dRI1_u <- gamma_I0*I1_u #undiagnosed
      dRI1_qu <- gamma_I0*I1_q #mild quarantined
      dRI2_q <- gamma_II*I2_iso #diagnosed
      dRI2_u <- gamma_I0*I2_u #undiagnosed
      dRI2_qu <- gamma_I0*I2_q #mild quarantined
      
      #death pathway
      dRD1 <- gamma_D*I1_s
      dRD2 <- gamma_D*I2_s
      
      
      list(c(dS1, dS2, dE1, dE2,
             dQ1, dQ2, dQI1, dQI2,
             dI1_d, dI2_d, dI1_u, dI2_u,
             dI1_iso, dI2_iso, dI1_q, dI2_q, dI1_s, dI2_s,
             dRI1_q, dRI1_u, dRI1_qu,dRD1,dRI2_q, dRI2_u, dRI2_qu,dRD2))
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
    recovery_I0=paras$recovery_I0 #recovery period for asymptomatic
    recovery_D=paras$recovery_D #death pathway: icu-death
    
    #case-fatality rate
    pD=paras$pD
    
    #Reff
    R0_c=paras$R0_c # community transmission
    R0_s=paras$R0_s# sexual transmission
    
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
  }
  
  #initial value, assuming population either susceptible or exposed--------------
  {
    init <-c(S1=N*p_susceptible1, S2=N*(1-p_susceptible2),
             E1=0, E2=0, Q1=0,Q2=0, QI1=0, QI2=0,
             I1_d=0, I2_d=0, I1_u=0, I2_u=0,
             I1_iso=0, I2_iso=0, I1_q=0, I2_q=0, I1_s=0, I2_s=0, 
             RI1_q=0, RI1_u=0, RI1_qu=0, RD1=0,
             RI2_q=0, RI2_u=0, RI2_qu=0, RD2=0)
  }
  
  
  #modelling------------
  {
    parameters <- c(beta_s = beta_s, beta_c=beta_c, alpha=alpha,
                    gamma_I0=gamma_I0, gamma_II = gamma_II,
                    gamma_D=gamma_D,
                    phi_I=phi_I,
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
    quarantined=rep(0,time) #number of quarantined individuals at each time point
    for(t in 1:time)
    {
      #importation
      init1['E1']=init1['E1']+imp[t]*p_susceptible_a0; init1['S1']=init1['S1']-imp[t]*p_susceptible_a0
      init1['E2']=init1['E2']+imp[t]*(1-p_susceptible_a0); init1['S2']=init1['S2']-imp[t]*(1-p_susceptible_a0)
      
      #quarantine close contacts (local)
      if(t>3 & ct_max>0){
        q_a0=sum(init1[c('E1','I1_d','I1_u')])
        q_l0=sum(init1[c('E2','I2_d','I2_u')])
        q0=min(sum(case_ct[t-3,]),ct_max)
        q_a1=min(q0*case_ct[t-3,1]/sum(case_ct[t-3,]),case_ct[t-3,1])
        q_l1=min(q0-q_a1,case_ct[t-3,2])
        q_a=min(q_a1*ct_connect,q_a0)
        q_l=min(q_l1*ct_connect,q_l0)
        init1['Q1']=init1['Q1']+q_a/q_a0*init1['E1']; init1['E1']=init1['E1']*(1-q_a/q_a0)
        init1['Q2']=init1['Q2']+q_l/q_l0*init1['E2']; init1['E2']=init1['E2']*(1-q_l/q_l0)
        init1['QI1']=init1['QI1']+q_a/q_a0*init1['I1_d']; init1['I1_d']=init1['I1_d']*(1-q_a/q_a0)
        init1['QI2']=init1['QI2']+q_l/q_l0*init1['I2_d']; init1['I2_d']=init1['I2_d']*(1-q_l/q_l0)
        init1['I1_q']=init1['I1_q']+q_a/q_a0*init1['I1_u']; init1['I1_u']=init1['I1_u']*(1-q_a/q_a0)
        init1['I2_q']=init1['I2_q']+q_l/q_l0*init1['I2_u']; init1['I2_u']=init1['I2_u']*(1-q_l/q_l0)
        quarantined[t]=q_a+q_l
      }
      prev_case=c(sum(init1[c('Q1','QI1','I1_iso','I1_q','I1_s','RI1_q','RI1_qu','RD1')]),
                  sum(init1[c('Q2','QI2','I2_iso','I2_q','I2_s','RI2_q','RI2_qu','RD2')]))
      out <- ode(y = init1, times = seq(0, 1), func = seir, parms = parameters)
      out1=data.frame(out)%>%filter(time==1)%>%select(-time)%>%as.matrix()
      names(out1)=names(init)
      init1=out1
      out2=rbind(out2,out1)
      #derive number of new, unquarantined local infection
      curr_case=c(sum(out1[c('Q1','QI1','I1_iso','I1_q','I1_s','RI1_q','RI1_qu','RD1')]),
                  sum(out1[c('Q2','QI2','I2_iso','I2_q','I2_s','RI2_q','RI2_qu','RD2')]))
      
      wild_case=lapply(1:2, function(i){
        max(curr_case[i]-prev_case[i],0)
      })%>%unlist()
      case_ct=rbind(case_ct,wild_case)
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
    
    #infections (including asymptomatic)
    cum_inf=rowSums(out2[,c('I1_d','I2_d','I1_u','I2_u','QI1','QI2',
                            'I1_iso','I2_iso','I1_q','I2_q','I1_s','I2_s',
                            'RI1_q','RI1_u','RI1_qu','RD1',
                            'RI2_q','RI2_u','RI2_qu','RD2')])
    new_inf=cum_inf-c(0,cum_inf[-time])
    
    #reported cases (after notification)
    cum_case=rowSums(out2[,c('I1_iso','I2_iso','I1_s','I2_s','RI1_q','RD1','RI2_q','RD2')])
    new_case=cum_case-c(0,cum_case[-time])
    
    #cumulative deaths
    cum_d=rowSums(out2[,c('RD1','RD2')])
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
    full_stat=data.frame(time=1:time, cum_inf=cum_inf,new_inf=new_inf,
                         new_case=new_case, cum_case=cum_case,
                         cum_prev=cum_inf/N*100,
                         
                         active_i=active_i, cum_deaths=cum_d)
    list(summary_stat=summary_stat,full_stat=full_stat)
  }
}