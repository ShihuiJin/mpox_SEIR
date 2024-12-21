#functions utilised in simulations
#simulation for different scenarios--------------
mpox_sim=function(c,Time=1800, sens=F, risky=0.1){
  c_info=city_info%>%filter(ISOALPHA==c)
  paras$N=c_info$population
  paras$time=Time
  paras$p_susceptible=c_info$p_susceptible
  
  if(sens){
    age_dist_c=age_dist%>%filter(ISOALPHA==c)%>%pull(age_dist)
    sus_c=immunity%>%filter(ISOALPHA==c)%>%pull(mu)
    sus_c=1-sus_c*0.807/100
    arange=c(4:10) #sexually active population age range
    paras$p_susceptible1=sum(sus_c[arange]*age_dist_c[arange])*risky
    paras$p_susceptible_a=sum(sus_c[arange]*age_dist_c[arange])/sum(age_dist_c[arange])
    paras$p_susceptible_a0=sum(age_dist_c[arange])*risky
    paras$p_active=paras$p_susceptible1/paras$p_susceptible
  }else{
    paras$p_susceptible1=c_info$p_susceptible1
    paras$p_susceptible_a=c_info$p_susceptible_a
    paras$p_susceptible_a0=c_info$p_susceptible_a0
    paras$p_active=c_info$p_active
  }
  res=quarantine_sim(paras)
  res
}

#functions for output derivation-----------
#summary statistics
estim=function(x){
  c(mean(x), quantile(x, c(0.5,0.025,0.975,0.25,0.75)))
}
#doubling time matrix
doubling_matrix=function(out_full, case_thre=2^c(0:10)){
  out=as.matrix(do.call('rbind',lapply(countries, function(c){
    res=out_full[[c]]
    time0=unlist(lapply(case_thre, function(i){ #first 128 local cases
      which(res$full_stat$cum_inf>=i)[1]
    }))
    doubling=(time0-lag(time0))[-1]
    c(1/mean(1/doubling, na.rm = T), #harmonic mean-->sensitive to importation frequency
      mean(doubling, na.rm = T),time0[1],doubling) #arithmatic mean
  })))%>%round(1)
  colnames(out)=c('mean_harmonic','mean_arith',paste0('T',1:length(case_thre)-1))
  out
}

#time to the first xxx cases
#x is the number of new cases
first_time=function(x, thre=10,n.sim=1e3){
  n.t=length(x)
  mclapply(1:n.sim, function(i){
    set.seed(i)
    cum=rpois(n.t, x)%>%cumsum()
    out=which(cum>=thre)[1]
    if(length(out)==0) out=n.t
    out
  },mc.cores=n.cores)%>%unlist()%>%estim()
}
timing=function(out_full, month=T, thre=1){
  unlist(lapply(countries, function(c){
    res=out_full[[c]]
    #t0=which(res$full_stat$cum_case>=thre)[1]
    t0=first_time(res$full_stat$new_case, thre)[2]
    #t0=res$summary_stat['t0_cum']
    if(month){
      round(t0/30,1)
    }else{
      t0
    }
  }))
}

#cumulative case counts
case_count=function(out_full, months=3){
  t1=months*30
  unlist(lapply(countries, function(c){
    res=out_full[[c]]
    t2=length(out_full[[c]]$full_stat$cum_case)
    res$full_stat$cum_case[min(t1,t2)]
  }))
}

#cumulative number of potential deaths
death_count=function(out_full,months=3){
  t1=months*30
  unlist(lapply(countries, function(c){
    res=out_full[[c]]
    t2=length(out_full[[c]]$full_stat$cum_case)
    res$full_stat$cum_deaths[min(t1,t2)]
  }))
}

#summary statistics derived from raw outputs
#time to index case, cases/deaths by xx days
out_sum=function(out_full, delta_t=6, scale=F, thre=1){
  lapply(1:3, function(k){
    out_sum=as.matrix(do.call('cbind', lapply(1:length(out_full), function(i){
      if(k==1) out=timing(out_full[[i]], F, thre)
      if(k==2){
        N=rep(1,nrow(city_info))
        if(scale) N=city_info$population/1e6
        out=case_count(out_full[[i]],delta_t)/N
      }
      if(k==3) out=death_count(out_full[[i]],delta_t)
      out
    })))
  })
}

out_sum_compare=function(summ,n.iso=4,n.quar=2){
  base=summ[,1]
  lapply(1:n.quar, function(i){
    as.matrix(do.call('cbind',lapply(1:n.iso, function(j){
      1-summ[,j+(i-1)*(n.iso+1)+1]/base
    })))*100
  })
}


