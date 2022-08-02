myweightedsum=function(DD,ix)
{
  return(sum(DD[ix,1],na.rm=T))
}


addrisk = function(contactsperson)
{
  contactsperson %>%
    filter(dcat>=0)%>%
    filter(!is.na(age)) %>%
    mutate(agealspac = ifelse(age<30,1,ifelse(age>=80,6,floor(30:79/10)-1))) %>%   #alspac age category
    mutate(probLFT = riskmatrix$willtest[agealspac]) %>%
    mutate(probretail = riskmatrix$limitexposure[agealspac]) %>%
    mutate(probtravel = riskmatrix$nobus[agealspac]) %>%
    mutate(probmask = riskmatrix$mask[agealspac]) %>%
    mutate(probWFH = riskmatrix$WFH[agealspac])-> contactsperson
  
  return(contactsperson)
}


degreeperson_calcR = function(degred,input,mult)
{  
  degred %>%
    mutate(rsquare1 = mult[1]*wcon*wcon) %>%
    mutate(rind1 = mult[1]*wcon*sum(wcon)) %>%
    mutate(rvac1_unvac = mult[1]*wcon,
           rvac1_dose1=mult[1]*wcon*(1-input$ve_inf),
           rvac1_dose2=mult[1]*wcon*(1-input$ve_inf2)) %>%
    mutate(rvac_unvac = mult[1]*wcon,
           rvac_d1az=mult[1]*wcon*(1-input$ve_i_az),
           rvac_d2az=mult[1]*wcon*(1-input$ve_i2_az),
           rvac_d1pf=mult[1]*wcon*(1-input$ve_i_pf),
           rvac_d2pf=mult[1]*wcon*(1-input$ve_i2_pf),
           rvac_d1mod=mult[1]*wcon*(1-input$ve_i_mod),
           rvac_d2mod=mult[1]*wcon*(1-input$ve_i2_mod),
           r_nat=mult[1]*wcon*(1-input$ve_i_nat),
           rvac_boost=mult[1]*wcon*(1-input$ve_i_boost), 
           r_vacnat=mult[1]*wcon*(1-input$ve_i_vacnat)) %>%
    mutate(cj_unvac = wcon, 
           cj_d1 = (1-input$ve_trans)*wcon, 
           cj_d2 = (1-input$ve_trans2)*wcon) %>%
    mutate(cj_unvac = wcon, cj_d1az = (1-input$ve_t_az)*wcon, 
           cj_d2az = (1-input$ve_t2_az)*wcon, 
           cj_d1pf = (1-input$ve_t_pf)*wcon, cj_d2pf = (1-input$ve_t2_pf)*wcon,
           cj_d1mod = (1-input$ve_t_mod)*wcon, cj_d2mod = (1-input$ve_t2_mod)*wcon,
           cj_nat = (1-input$ve_t_nat)*wcon,
           cj_vacnat = (1-input$ve_t_vacnat)*wcon,
           cj_boost = (1-input$ve_t_boost)*wcon) %>%
    mutate(rvsquare = (rd0*cj_unvac + natinf*cj_nat + rboost*cj_boost + vacnat*cj_vacnat + 
                         rd1az*cj_d1az + rd2az*cj_d2az + 
                         rd1pf*cj_d1pf + rd2pf*cj_d2pf+ 
                         rd1mod*cj_d1mod + rd2mod*cj_d2mod)*(rd0*rvac_unvac + rvac_d1az*rd1az + 
                                                               rvac_d2az*rd2az + rvac_d1pf*rd1pf + rvac_d2pf*rd2pf + 
                                                               rvac_d1mod*rd1mod + rvac_d2mod*rd2mod + 
                                                               r_nat*natinf + rvac_boost*rboost)       ) %>%
    #mutate(RR1 = sum(rd0*cj_unvac*(1-s1_unvac)+ rd1*cj_d1*(1-s1_d1) + rd2*cj_d2*(1-s1_d2)))%>%
    mutate(RR1 = sum(rd0*cj_unvac*(1-s1_unvac) + natinf*cj_nat*(1-s1_nat) + rboost*cj_boost*(1-s1_boost) + vacnat*cj_vacnat*(1-s1_vacnat) +
                       rd1az*cj_d1az*(1-s1_d1az) + rd2az*cj_d2az*(1-s1_d2az) + 
                       rd1pf*cj_d1pf*(1-s1_d1pf) + rd2pf*cj_d2pf*(1-s1_d2pf)+ 
                       rd1mod*cj_d1mod*(1-s1_d1mod) + rd2mod*cj_d2mod*(1-s1_d2mod)))%>%
    mutate(s1_unvac = exp(-rvac1_unvac*RR1), 
           s1_d1 = exp(-rvac1_dose1*RR1),
           s1_d2 = exp(-rvac1_dose2*RR1)) %>%
    mutate(s1_unvac = exp(-rvac_unvac*RR1), 
           s1_d1az = exp(-rvac_d1az*RR1),
           s1_d2az = exp(-rvac_d2az*RR1),
           s1_d1pf = exp(-rvac_d1pf*RR1),
           s1_d2pf = exp(-rvac_d2pf*RR1),
           s1_d1mod = exp(-rvac_d1mod*RR1),
           s1_d2mod = exp(-rvac_d2mod*RR1),
           s1_nat = exp(-r_nat*RR1),
           s1_vacnat = exp(-r_vacnat*RR1),
           s1_boost = exp(-rvac_boost*RR1)) %>%
    mutate(cases = ppp*(rd0*(1-s1_unvac) + rd1az*(1-s1_d1az) + rd2az*(1-s1_d2az) + rd1pf*(1-s1_d1pf) + rd2pf*(1-s1_d2pf) + rd1mod*(1-s1_d1mod) + rd2mod*(1-s1_d2mod) + 
                          natinf*(1-s1_nat) + rboost*(1-s1_boost) + natinf*(1-s1_vacnat))) %>%
    mutate(D1unvac = mortality[age5]*ppp*(rd0*(1-s1_unvac) + natinf*(1-s1_nat)*(1-input$ve_deathrednat))) %>%
    mutate(D1r1 = ppp*mortality[age5]*((rd1az*(1-s1_d1az)*(1-input$ve_deathredaz1)+rd1pf*(1-s1_d1pf)*(1-input$ve_deathredpf1)+rd1mod*(1-s1_d1mod)*(1-input$ve_deathredmod1)))) %>%
    mutate(D1r2 = ppp*mortality[age5]*(rd2az*(1-s1_d2az)*(1-input$ve_deathredaz2)+rd2pf*(1-s1_d2pf)*(1-input$ve_deathredpf2)+rd2mod*(1-s1_d2mod)*(1-input$ve_deathredmod2) + 
                                         r_vacnat*(1-s1_vacnat)*(1-input$ve_deathredvacnat))) %>%
    mutate(D1boost = ppp*mortality[age5]*rboost*(1-s1_boost)*(1-input$ve_deathredboost)) %>%
    #mutate(hospadmissions = hosprate[age5]*ppp*(rd0*(1-s1_unvac) + 
    #                                       rd1az*(1-s1_d1az)*(1-input$ve_deathredaz1) + rd2az*(1-s1_d2az)*(1-input$ve_deathredaz2) + 
    #                                      rd1mod*(1-s1_d1mod)*(1-input$ve_deathredmod1) + rd2mod*(1-s1_d2mod)*(1-input$ve_deathredmod2) + 
    #                                       rd1pf*(1-s1_d1pf)*(1-input$ve_deathredpf1) + rd2pf*(1-s1_d2pf)*(1-input$ve_deathredpf2) + 
    #                                       natinf*(1-s1_nat)*(1-input$ve_deathrednat) + rboost*(1-s1_boost)*(1-input$ve_deathredboost) + 
    #                                       r_vacnat*(1-s1_vacnat)*(1-input$ve_deathredvacnat)))%>%
    mutate(hospadmissions = hosprate[age5]*ppp*(rd0*(1-s1_unvac) + 
                                                  rd1az*(1-s1_d1az)*(1-input$ve_deathredaz1) + rd2az*(1-s1_d2az)*(1-input$ve_deathredaz2) + 
                                                  rd1pf*(1-s1_d1pf)*(1-input$ve_deathredpf1) + rd2pf*(1-s1_d2pf)*(1-input$ve_deathredpf2) + 
                                                  rd1mod*(1-s1_d1mod)*(1-input$ve_deathredmod1) + rd2mod*(1-s1_d2mod)*(1-input$ve_deathredmod2) + 
                                                  natinf*(1-s1_nat)*(1-input$ve_deathrednat) + rboost*(1-s1_boost)*(1-input$ve_deathredboost) + 
                                                  r_vacnat*(1-s1_vacnat)*(1-input$ve_deathredvacnat))) %>%
    mutate(deaths = mortality[age5]*ppp*(rd0*(1-s1_unvac) + 
                                           rd1az*(1-s1_d1az)*(1-input$ve_deathredaz1) + rd2az*(1-s1_d2az)*(1-input$ve_deathredaz2) + 
                                           rd1pf*(1-s1_d1pf)*(1-input$ve_deathredpf1) + rd2pf*(1-s1_d2pf)*(1-input$ve_deathredpf2) + 
                                           rd1mod*(1-s1_d1mod)*(1-input$ve_deathredmod1) + rd2mod*(1-s1_d2mod)*(1-input$ve_deathredmod2) + 
                                           natinf*(1-s1_nat)*(1-input$ve_deathrednat) + rboost*(1-s1_boost)*(1-input$ve_deathredboost) + 
                                           r_vacnat*(1-s1_vacnat)*(1-input$ve_deathredvacnat))) -> degred
  
  
  
  return(degred)
}


degreeperson = function(contactsperson,input,mult)
{
  SAR = runif(length(contactsperson$PID),min=input$SARlo,max=input$SARhi)
  SARhome = runif(length(contactsperson$PID),min=input$SARhomelo,max=input$SARhomehi)
  contactsperson %>%
    filter(dcat>=0)%>%
    filter(!is.na(age)) %>%
    mutate(entrynum=row_number())%>%
    mutate(mySAR=SAR[entrynum], mySARhome = SARhome[entrynum]) %>%
    mutate(trans=ifelse(age<=11 | Home==1 | wearsmask==0,1,(1-input$covidsec))) %>%
    mutate(infect = ifelse(age<=11,input$infchild*trans,trans)) %>%
    mutate(Num1 = ifelse(contacttrace==TRUE | LFT==TRUE,(1-input$CTF1)*Numbers*(1-comply)*infect*mySAR,Numbers*(1-comply)*infect*mySAR)) %>%
    mutate(Num2 = ifelse(contacttrace==TRUE | LFT==TRUE,(1-input$CTF1)*Numbers*(1-comply)*infect*mySARhome,Numbers*(1-comply)*infect*mySARhome)) %>%
    mutate(homenum = Num2*Home, worknum=Num1*WorkSchool, 
           othernum=Num1*OtherLeisure+Num1*Travel)%>%
    group_by(PID)%>%
    summarise(degree=sum(homenum)+sum(worknum)+sum(othernum),
              totdcat=sum(dcat),
              wcon=sum(homenum*dcat2)+sum(worknum*dcat2)+sum(othernum*dcat2),age=mean(age),
              ageweight=mean(ageweight))%>%
    mutate(rsquare1 = wcon*wcon,aweight=ageweight/sum(ageweight))-> degred
  
  degred %>%
    filter(!is.na(age)) %>%
    mutate(s1_unvac=0.5, s1_d1=0.5, s1_d2=0.5) %>%
    mutate(s1_d1az=0.5, s1_d2az=0.5, s1_d1pf=0.5, s1_d2pf=0.5, s1_d1mod=0.5, s1_d2mod=0.5) %>%
    mutate(s1_nat=0.5, s1_boost=0.5, s1_vacnat=0.5) %>%
    mutate(age5 = ifelse(age<=17,1+floor(age/5),ifelse(age<=90,2+floor(age/5),2+floor(90/5)))) %>%   #age category
    mutate(ppp = agenumbers$ppp[age5]) %>%          #number of people represented by this person
    mutate(rboost = vaccinestatus$booster[age5]) %>%           #proportion that have received a booster
    mutate(rd2az = vaccinestatus$az2[age5]) %>%                                                     #proportion that have received 2 doses AZ
    mutate(rd1az = vaccinestatus$az1[age5] ) %>%         #proportion that have received 1 dose AZ
    mutate(rd2pf = vaccinestatus$pf2[age5]) %>%                                                     #proportion that have received 2 doses Pfizer
    mutate(rd1pf = vaccinestatus$pf1[age5] ) %>%         #proportion that have received 1 dose Pfizer
    mutate(rd2mod = vaccinestatus$mod2[age5]) %>%                                                    #proportion that have received 2 doses Moderna
    mutate(rd1mod = vaccinestatus$mod1[age5] ) %>% #proportion that have received 1 dose Moderna
    mutate(rd0 = vaccinestatus$unprotected[age5]) %>%                        #proportion that are unvaccinated
    mutate(natinf = vaccinestatus$natimm[age5]) %>%                                          #proportion with natural immunity only
    mutate(vacnat = vaccinestatus$vacnat[age5]) -> degred
  #degred = degreeperson_calcR(degred,input,mult)
  
  return(degred)
  
}

degreeperson0 = function(contactsperson,input)
{
  SAR = runif(length(contactsperson$PID),min=input$SARlo,max=input$SARhi)
  SARhome = runif(length(contactsperson$PID),min=input$SARhomelo,max=input$SARhomehi)
  contactsperson %>%
    filter(dcat>=0)%>%
    filter(!is.na(age)) %>%
    mutate(entrynum=row_number())%>%
    mutate(mySAR=SAR[entrynum], mySARhome = SARhome[entrynum]) %>%
    mutate(agealspac = ifelse(age<30,1,ifelse(age>=80,6,floor(30:79/10)-1))) %>%   #alspac age category
    mutate(probLFT = riskmatrix$willtest[agealspac]) %>%
    mutate(probretail = riskmatrix$limitexposure[agealspac]) %>%
    mutate(probtravel = riskmatrix$nobus[agealspac]) %>%
    mutate(probmask = riskmatrix$mask[agealspac]) %>%
    mutate(probWFH = riskmatrix$WFH[agealspac]) %>%
    mutate(trans=1) %>% #SAR in households for delta
    mutate(infect = ifelse(age<=11,input$infchild*trans,trans)) %>%
    #mutate(infect = ifelse(age<=11,input$infchild,ifelse(age<=18,max(1,2*input$infchild),1))) %>%
    mutate(Num1 = ifelse(contacttrace==TRUE,(Numbers)*(1-comply)*infect*mySAR,Numbers*(1-comply)*infect*mySAR)) %>%
    mutate(Num2 = ifelse(contacttrace==TRUE,(Numbers)*(1-comply)*infect*mySARhome,Numbers*(1-comply)*infect*mySARhome)) %>%
    mutate(homenum = Num2*Home, worknum=Num1*WorkSchool, 
           othernum=Num1*OtherLeisure+Num1*Travel)%>%
    group_by(PID)%>%
    summarise(degree=sum(homenum)+sum(worknum)+sum(othernum),
              totdcat=sum(dcat),
              wcon=sum(homenum*dcat2)+sum(worknum*dcat2)+sum(othernum*dcat2),age=mean(age),
              ageweight=mean(ageweight))%>%
    mutate(rsquare = wcon*wcon,aweight=ageweight/sum(ageweight))-> degred
  
  return(degred)
}

