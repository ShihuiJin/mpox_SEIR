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

#cumulative number of infections
inf_count=function(out_full, months=3){
  t1=months*30
  unlist(lapply(countries, function(c){
    res=out_full[[c]]
    res$full_stat$cum_inf[t1]
  }))
}

#cumulative case counts
case_count=function(out_full, months=3){
  t1=months*30
  unlist(lapply(countries, function(c){
    res=out_full[[c]]
    res$full_stat$cum_case[t1]
  }))
}

#cumulative number of potential deaths
death_count=function(out_full,months=3){
  t1=months*30
  unlist(lapply(countries, function(c){
    res=out_full[[c]]
    res$full_stat$cum_deaths[t1]
  }))
}

#summary statistics derived from raw outputs
#infection size/cases/deaths by xx days
out_sum=function(out_full, delta_t=6, scale=F){
  lapply(1:3, function(k){
    out_sum=as.matrix(do.call('cbind', lapply(1:length(out_full), function(i){
      if(k==1)
        out=inf_count(out_full[[i]],delta_t)/country_info$population*1e6
      if(k==2){
        N=rep(1,nrow(country_info))
        if(scale) N=country_info$population/1e6
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


