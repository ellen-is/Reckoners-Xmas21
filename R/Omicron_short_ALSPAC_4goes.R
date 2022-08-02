library(tidyverse)
library(readxl)
library(boot)
library(fields)
library(cowplot)
library(metR)
library(patchwork)
library(colorspace)
library(ggtext)

RUNMODEL=1
if(RUNMODEL==1)
{

  source("R/parameters.R")
  source("R/loaddata.R")
  source("R/reckonerfunctions.R")

mycols2 = rainbow(3,start=3/6,end=5/6,alpha=0.5)
mycols3 = rainbow(3,start=0/6,end=2/6,alpha=0.5)


baseR=c(1,5,3)
degall=degreeperson0(contactsperson,input)
SCSdata = contactsperson

mult=rep(0,3)
mult[1]=baseR[1]/sum(degall$rsquare)
mult[2]=baseR[2]/sum(degall$rsquare)
mult[3]=baseR[3]/sum(degall$rsquare)

sum(degall$rsquare)*mult[1]  
#check works with boot
myboot=boot(degall[,c("rsquare","age")],statistic=myweightedsum,R=NR) 
mult[2]*mean(myboot$t)
mult[3]*mean(myboot$t)
mult[1]*mean(myboot$t)
#perfect! 


vaccinestatus=vaccinestatus_new

#find the contact over 18
ixO18=which(contactsperson$age>18)
#find the school contacts
ixU18=which(contactsperson$age<=18 & contactsperson$WorkSchool==1)
#find leisure contacts for children
ixU18_O=which(contactsperson$age<=18 & contactsperson$WorkSchool==0)
#find home contacts
ixH = contactsperson$Home==1
ixT = contactsperson$Travel==1
ixW = contactsperson$Work==1
ixShops = contactsperson$Other==1
numreps=1
cc=0; 
x=1




##########RUN FOUR SCENARIOS
numscenarios = 4

escape=c(4,4,4,4)
escape_severity=c(4,4,1,1)
severity=c(0.8,0.6,0.8,0.6)

variant=2 #1=delta, 2=omicron
outdf = data.frame(scenario=rep(1:numscenarios,each=NR*numreps),simnum=rep(1:numreps,numscenarios*NR),
                   bootnum=rep(1:NR,numscenarios*numreps),Reff=0,Infections=0,Cases=0,Deaths=0,
                   Reffrisk=0,InfectionsRisk=0,Casesrisk=0,Deathsrisk=0)


numreps=10

i=1

for(SC in 1:numscenarios)
{
  mortality = severity[SC]*mortality2
  for(r1 in 1:numreps)
  {
    input = input2
    input[susix]=input2[susix]^escape[SC]
    input[infix]=input2[infix]^escape[SC]
    input[deathix]=input2[deathix]^escape_severity[SC]
    
    input$SAR = SAR_notHH[variant]
    input$SARhome = SAR_HH[variant]
    input$SARlo = SAR_notHH1[variant]
    input$SARhomelo = SAR_HH1[variant]
    input$SARhi = SAR_notHH2[variant]
    input$SARhomehi = SAR_HH2[variant]
  
    contactsperson = SCSdata
    degall=degreeperson0(contactsperson,input)
    mult[1]=R0/sum(degall$rsquare)
    
    contactsperson$contacttrace = 0*contactsperson$dcat
    contactsperson$LFT = 0*contactsperson$dcat
    contactsperson$wearsmask = 0*contactsperson$dcat
    contactsperson$reduceretail = 0*contactsperson$dcat
    contactsperson$reducetravel = 0*contactsperson$dcat
    contactsperson$WFH = 0*contactsperson$dcat  
    #compliance flag: 0 means contact takes place; 1 means it doesn't
    contactsperson$comply = 0*contactsperson$dcat
    
    ##contact tracing flag based on probability of symptoms
    myrands = runif(length(contactsperson$hassymp),min=0,max=1)  #not sure if I should draw a random number each time
    contactsperson$contacttrace = myrands < contactsperson$hassymp
    
    contactsperson$comply[ixU18]=0 #school contacts do/do not take place
    contactsperson$comply[ixH]=0 #home contacts do take place
    
    #mortality = mortality2*(1+escape
  degred=degreeperson(contactsperson,input,mult)
  degred = degreeperson_calcR(degred,input,mult)
  #degred
  for(r2 in 1:20){degred = degreeperson_calcR(degred,input,mult)}
  print(sum(degred$rvsquare))
  print(sum(degred$deaths))
  
  myboot=boot(cbind(degred[,"rvsquare"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
  outdf$Reff[outdf$scenario==SC & outdf$simnum==r1] = myboot$t
  
  myboot=boot(cbind(degred[,"deaths"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
  outdf$Deaths[outdf$scenario==SC & outdf$simnum==r1] = myboot$t
  
  myboot=boot(cbind(degred[,"hospadmissions"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
  #myboot=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
  outdf$Cases[outdf$scenario==SC & outdf$simnum==r1] = myboot$t
  
  myboot=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
  outdf$Infections[outdf$scenario==SC & outdf$simnum==r1] = myboot$t
  head(outdf)
  
  #WITH RISK MITIGATIONS 
  contactsperson = addrisk(contactsperson)
  contactsperson$LFT = myrands < contactsperson$probLFT*LFTsensitivity
  contactsperson$wearsmask = myrands < contactsperson$probmask
  
  contactsperson$reduceretail = myrands < contactsperson$probretail #this is 1 if retail is reduced and 0 if not
  contactsperson$reducetravel = myrands < contactsperson$probtravel
  contactsperson$WFH = myrands < contactsperson$probWFH
  
  contactsperson$comply[contactsperson$WFH*ixW==1] = 1
  contactsperson$comply[contactsperson$reduceretail*ixShops==1] = 1
  contactsperson$comply[contactsperson$reducetravel*ixT==1] = 1
  #contactsperson$comply[ixU18]=1 #school contacts do/do not take place
  
  degred=degreeperson(contactsperson,input,mult)
  degred = degreeperson_calcR(degred,input,mult)
  for(r2 in 1:20){degred = degreeperson_calcR(degred,input,mult)}
  print(sum(degred$rvsquare))
  print(sum(degred$deaths))
  
  myboot2=boot(cbind(degred[,"rvsquare"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
  outdf$Reffrisk[outdf$scenario==SC & outdf$simnum==r1] = myboot2$t

  myboot3=boot(cbind(degred[,"deaths"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
  outdf$Deathsrisk[outdf$scenario==SC & outdf$simnum==r1] = myboot3$t
  
  myboot2=boot(cbind(degred[,"hospadmissions"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
  #myboot2=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
  outdf$Casesrisk[outdf$scenario==SC & outdf$simnum==r1] = myboot2$t
  
  myboot2=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
  outdf$Infectionsrisk[outdf$scenario==SC & outdf$simnum==r1] = myboot2$t
  
    i=i+1
}
}




##########2222222222222222222
#####changing susceptibility and infectiousness
Rvalues=seq(7,11,1)
transadvantage=c(2)
escape=c(4,4)
SAR_notHH=c(0.03,0.087)
SAR_HH=c(0.103,0.158)

SAR_notHH1=c(0.028,0.075)
SAR_HH1=c(0.101,0.143)
SAR_notHH2=c(0.032,0.10)
SAR_HH2=c(0.105,0.175)
#transadvantage=c(1)
LFTsensitivity=0.5
numreps=10
outdf2 = data.frame(tau=rep(transadvantage,each=NR*numreps),simnum=rep(1:numreps,length(transadvantage)*NR),
                   bootnum=rep(1:NR,length(transadvantage)*numreps),Reff=0,Infections=0,Cases=0,Deaths=0,Reffrisk=0,Infectionsrisk=0,Casesrisk=0,Deathsrisk=0)

i=1
susix = seq(5,26,2)
infix = seq(6,26,2)
deathix = c(27:37)
R0=7
mortality = 0.6*mortality2
mortality = 0.5*mortality2
for(tau in transadvantage)
{
  for(r1 in 1:numreps)
  {
    input = input2
    input[susix]=input2[susix]^escape[tau]
    input[infix]=input2[infix]^escape[tau]
    input[deathix]=input2[deathix]^escape[tau]
    #input[susix]=0
    #input[infix]=0
    
    input$SAR = SAR_notHH[tau]
    input$SARhome = SAR_HH[tau]
    input$SARlo = SAR_notHH1[tau]
    input$SARhomelo = SAR_HH1[tau]
    input$SARhi = SAR_notHH2[tau]
    input$SARhomehi = SAR_HH2[tau]
    
    contactsperson = SCSdata
    degall=degreeperson0(contactsperson,input)
    mult[1]=R0/sum(degall$rsquare)
    
    contactsperson$contacttrace = 0*contactsperson$dcat
    contactsperson$LFT = 0*contactsperson$dcat
    contactsperson$wearsmask = 0*contactsperson$dcat
    contactsperson$reduceretail = 0*contactsperson$dcat
    contactsperson$reducetravel = 0*contactsperson$dcat
    contactsperson$WFH = 0*contactsperson$dcat  
    #compliance flag: 0 means contact takes place; 1 means it doesn't
    contactsperson$comply = 0*contactsperson$dcat
    
    ##contact tracing flag based on probability of symptoms
    myrands = runif(length(contactsperson$hassymp),min=0,max=1)  #not sure if I should draw a random number each time
    contactsperson$contacttrace = myrands < contactsperson$hassymp
    
    contactsperson$comply[ixU18]=0 #school contacts do/do not take place
    contactsperson$comply[ixH]=0 #home contacts do take place
    
    #mortality = mortality2*(1+escape
    degred=degreeperson(contactsperson,input,mult)
    degred = degreeperson_calcR(degred,input,mult)
    #degred
    for(r2 in 1:20){degred = degreeperson_calcR(degred,input,mult)}
    print(sum(degred$rvsquare))
    print(sum(degred$deaths))
    
    myboot=boot(cbind(degred[,"rvsquare"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf2$Reff[outdf2$tau==tau & outdf2$simnum==r1] = myboot$t
    
    myboot=boot(cbind(degred[,"deaths"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf2$Deaths[outdf2$tau==tau& outdf2$simnum==r1] = myboot$t
    
    myboot=boot(cbind(degred[,"hospadmissions"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    #myboot=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf2$Cases[outdf2$tau==tau& outdf2$simnum==r1] = myboot$t
    
    myboot=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf2$Infections[outdf2$tau==tau& outdf2$simnum==r1] = myboot$t
    
    #WITH RISK MITIGATIONS 
    contactsperson = addrisk(contactsperson)
    contactsperson$LFT = myrands < contactsperson$probLFT*LFTsensitivity
    contactsperson$wearsmask = myrands < contactsperson$probmask
    
    contactsperson$reduceretail = myrands < contactsperson$probretail #this is 1 if retail is reduced and 0 if not
    contactsperson$reducetravel = myrands < contactsperson$probtravel
    contactsperson$WFH = myrands < contactsperson$probWFH
    
    contactsperson$comply[contactsperson$WFH*ixW==1] = 1
    contactsperson$comply[contactsperson$reduceretail*ixShops==1] = 1
    contactsperson$comply[contactsperson$reducetravel*ixT==1] = 1
    #contactsperson$comply[ixU18]=1 #school contacts do/do not take place
    
    degred=degreeperson(contactsperson,input,mult)
    degred = degreeperson_calcR(degred,input,mult)
    for(r2 in 1:20){degred = degreeperson_calcR(degred,input,mult)}
    print(sum(degred$rvsquare))
    print(sum(degred$deaths))
    
    myboot2=boot(cbind(degred[,"rvsquare"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf2$Reffrisk[outdf2$tau==tau & outdf2$simnum==r1] = myboot2$t
    
    myboot3=boot(cbind(degred[,"deaths"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf2$Deathsrisk[outdf2$tau==tau & outdf2$simnum==r1] = myboot3$t
    
    myboot4=boot(cbind(degred[,"hospadmissions"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    #myboot2=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf2$Casesrisk[outdf2$tau==tau & outdf2$simnum==r1] = myboot4$t
    
    myboot2=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf2$Infectionsrisk[outdf2$tau==tau & outdf2$simnum==r1] = myboot2$t
    
    i=i+1
  }
}



##########3333333333333333333333
#####changing susceptibility and infectiousness
Rvalues=seq(7,11,1)
transadvantage=c(2)
escape=c(4,4)
SAR_notHH=c(0.03,0.087)
SAR_HH=c(0.103,0.158)

SAR_notHH1=c(0.028,0.075)
SAR_HH1=c(0.101,0.143)
SAR_notHH2=c(0.032,0.10)
SAR_HH2=c(0.105,0.175)
#transadvantage=c(1)
LFTsensitivity=0.5
numreps=10
outdf3 = data.frame(tau=rep(transadvantage,each=NR*numreps),simnum=rep(1:numreps,length(transadvantage)*NR),
                   bootnum=rep(1:NR,length(transadvantage)*numreps),Reff=0,Infections=0,Cases=0,Deaths=0,
                   Reffrisk=0,Infectionsrisk=0,Casesrisk=0,Deathsrisk=0)

i=1
susix = seq(5,26,2)
infix = seq(6,26,2)
deathix = c(27:37)
R0=7
mortality = 0.8*mortality2
for(tau in transadvantage)
{
  for(r1 in 1:numreps)
  {
    input = input2
    input[susix]=input2[susix]^escape[tau]
    input[infix]=input2[infix]^escape[tau]
    #input[deathix]=input2[deathix]^escape[tau]
    #input[susix]=0
    #input[infix]=0
    
    input$SAR = SAR_notHH[tau]
    input$SARhome = SAR_HH[tau]
    input$SARlo = SAR_notHH1[tau]
    input$SARhomelo = SAR_HH1[tau]
    input$SARhi = SAR_notHH2[tau]
    input$SARhomehi = SAR_HH2[tau]
    
    contactsperson = SCSdata
    degall=degreeperson0(contactsperson,input)
    mult[1]=R0/sum(degall$rsquare)
    
    contactsperson$contacttrace = 0*contactsperson$dcat
    contactsperson$LFT = 0*contactsperson$dcat
    contactsperson$wearsmask = 0*contactsperson$dcat
    contactsperson$reduceretail = 0*contactsperson$dcat
    contactsperson$reducetravel = 0*contactsperson$dcat
    contactsperson$WFH = 0*contactsperson$dcat  
    #compliance flag: 0 means contact takes place; 1 means it doesn't
    contactsperson$comply = 0*contactsperson$dcat
    
    ##contact tracing flag based on probability of symptoms
    myrands = runif(length(contactsperson$hassymp),min=0,max=1)  #not sure if I should draw a random number each time
    contactsperson$contacttrace = myrands < contactsperson$hassymp
    
    contactsperson$comply[ixU18]=0 #school contacts do/do not take place
    contactsperson$comply[ixH]=0 #home contacts do take place
    
    #mortality = mortality2*(1+escape
    degred=degreeperson(contactsperson,input,mult)
    degred = degreeperson_calcR(degred,input,mult)
    #degred
    for(r2 in 1:20){degred = degreeperson_calcR(degred,input,mult)}
    print(sum(degred$rvsquare))
    print(sum(degred$deaths))
    
    myboot=boot(cbind(degred[,"rvsquare"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf3$Reff[outdf3$tau==tau & outdf3$simnum==r1] = myboot$t
    
    myboot=boot(cbind(degred[,"deaths"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf3$Deaths[outdf3$tau==tau& outdf3$simnum==r1] = myboot$t
    
    myboot=boot(cbind(degred[,"hospadmissions"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    #myboot=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf3$Cases[outdf3$tau==tau& outdf3$simnum==r1] = myboot$t
    
    myboot=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf3$Infections[outdf3$tau==tau& outdf3$simnum==r1] = myboot$t
    
    
    #WITH RISK MITIGATIONS 
    contactsperson = addrisk(contactsperson)
    contactsperson$LFT = myrands < contactsperson$probLFT*LFTsensitivity
    contactsperson$wearsmask = myrands < contactsperson$probmask
    
    contactsperson$reduceretail = myrands < contactsperson$probretail #this is 1 if retail is reduced and 0 if not
    contactsperson$reducetravel = myrands < contactsperson$probtravel
    contactsperson$WFH = myrands < contactsperson$probWFH
    
    contactsperson$comply[contactsperson$WFH*ixW==1] = 1
    contactsperson$comply[contactsperson$reduceretail*ixShops==1] = 1
    contactsperson$comply[contactsperson$reducetravel*ixT==1] = 1
    #contactsperson$comply[ixU18]=1 #school contacts do/do not take place
    
    degred=degreeperson(contactsperson,input,mult)
    degred = degreeperson_calcR(degred,input,mult)
    for(r2 in 1:20){degred = degreeperson_calcR(degred,input,mult)}
    print(sum(degred$rvsquare))
    print(sum(degred$deaths))
    
    myboot2=boot(cbind(degred[,"rvsquare"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf3$Reffrisk[outdf3$tau==tau & outdf3$simnum==r1] = myboot2$t
    
    myboot3=boot(cbind(degred[,"deaths"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf3$Deathsrisk[outdf3$tau==tau & outdf3$simnum==r1] = myboot3$t
    
    myboot2=boot(cbind(degred[,"hospadmissions"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    #myboot2=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf3$Casesrisk[outdf3$tau==tau & outdf3$simnum==r1] = myboot2$t
    
    myboot2=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf3$Infectionsrisk[outdf3$tau==tau & outdf3$simnum==r1] = myboot2$t
    
    i=i+1
  }
}

##########44444444444444444444444
#####changing susceptibility and infectiousness
Rvalues=seq(7,11,1)
transadvantage=c(2)
escape=c(4,4)
SAR_notHH=c(0.03,0.087)
SAR_HH=c(0.103,0.158)

SAR_notHH1=c(0.028,0.075)
SAR_HH1=c(0.101,0.143)
SAR_notHH2=c(0.032,0.10)
SAR_HH2=c(0.105,0.175)
#transadvantage=c(1)
LFTsensitivity=0.5
numreps=10
outdf4 = data.frame(tau=rep(transadvantage,each=NR*numreps),simnum=rep(1:numreps,length(transadvantage)*NR),
                    bootnum=rep(1:NR,length(transadvantage)*numreps),Reff=0,Infections=0,Cases=0,Deaths=0,
                    Reffrisk=0,Infectionsrisk=0,Casesrisk=0,Deathsrisk=0)

i=1
susix = seq(5,26,2)
infix = seq(6,26,2)
deathix = c(27:37)
R0=7
mortality = 0.5*mortality2
for(tau in transadvantage)
{
  for(r1 in 1:numreps)
  {
    input = input2
    input[susix]=(input2[susix]^escape[tau]) #to account for waning?
    input[infix]=(input2[infix]^escape[tau])
    input[deathix]=input2[deathix] #should we also account for waning in VE against death?
    #input[susix]=0
    #input[infix]=0
    
    input$SAR = SAR_notHH[tau]
    input$SARhome = SAR_HH[tau]
    input$SARlo = SAR_notHH1[tau]
    input$SARhomelo = SAR_HH1[tau]
    input$SARhi = SAR_notHH2[tau]
    input$SARhomehi = SAR_HH2[tau]
    
    contactsperson = SCSdata
    degall=degreeperson0(contactsperson,input)
    mult[1]=R0/sum(degall$rsquare)
    
    contactsperson$contacttrace = 0*contactsperson$dcat
    contactsperson$LFT = 0*contactsperson$dcat
    contactsperson$wearsmask = 0*contactsperson$dcat
    contactsperson$reduceretail = 0*contactsperson$dcat
    contactsperson$reducetravel = 0*contactsperson$dcat
    contactsperson$WFH = 0*contactsperson$dcat  
    #compliance flag: 0 means contact takes place; 1 means it doesn't
    contactsperson$comply = 0*contactsperson$dcat
    
    ##contact tracing flag based on probability of symptoms
    myrands = runif(length(contactsperson$hassymp),min=0,max=1)  #not sure if I should draw a random number each time
    contactsperson$contacttrace = myrands < contactsperson$hassymp
    
    contactsperson$comply[ixU18]=0 #school contacts do/do not take place
    contactsperson$comply[ixH]=0 #home contacts do take place
    
    #mortality = mortality2*(1+escape
    degred=degreeperson(contactsperson,input,mult)
    degred = degreeperson_calcR(degred,input,mult)
    #degred
    for(r2 in 1:20){degred = degreeperson_calcR(degred,input,mult)}
    print(sum(degred$rvsquare))
    print(sum(degred$deaths))
    
    myboot=boot(cbind(degred[,"rvsquare"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf4$Reff[outdf4$tau==tau & outdf4$simnum==r1] = myboot$t
    
    myboot=boot(cbind(degred[,"deaths"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf4$Deaths[outdf4$tau==tau& outdf4$simnum==r1] = myboot$t
    
    myboot=boot(cbind(degred[,"hospadmissions"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    #myboot=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf4$Cases[outdf4$tau==tau& outdf4$simnum==r1] = myboot$t
    
    myboot=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf4$Infections[outdf4$tau==tau& outdf4$simnum==r1] = myboot$t
    
    
    #WITH RISK MITIGATIONS 
    contactsperson = addrisk(contactsperson)
    contactsperson$LFT = myrands < contactsperson$probLFT*LFTsensitivity
    contactsperson$wearsmask = myrands < contactsperson$probmask
    
    contactsperson$reduceretail = myrands < contactsperson$probretail #this is 1 if retail is reduced and 0 if not
    contactsperson$reducetravel = myrands < contactsperson$probtravel
    contactsperson$WFH = myrands < contactsperson$probWFH
    
    contactsperson$comply[contactsperson$WFH*ixW==1] = 1
    contactsperson$comply[contactsperson$reduceretail*ixShops==1] = 1
    contactsperson$comply[contactsperson$reducetravel*ixT==1] = 1
    #contactsperson$comply[ixU18]=1 #school contacts do/do not take place
    
    degred=degreeperson(contactsperson,input,mult)
    degred = degreeperson_calcR(degred,input,mult)
    for(r2 in 1:20){degred = degreeperson_calcR(degred,input,mult)}
    print(sum(degred$rvsquare))
    print(sum(degred$deaths))
    
    myboot2=boot(cbind(degred[,"rvsquare"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf4$Reffrisk[outdf4$tau==tau & outdf4$simnum==r1] = myboot2$t
    
    myboot3=boot(cbind(degred[,"deaths"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf4$Deathsrisk[outdf4$tau==tau & outdf4$simnum==r1] = myboot3$t
    
    myboot2=boot(cbind(degred[,"hospadmissions"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    #myboot2=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf4$Casesrisk[outdf4$tau==tau & outdf4$simnum==r1] = myboot2$t
    
    myboot2=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf4$Infectionsrisk[outdf4$tau==tau & outdf4$simnum==r1] = myboot2$t
    
    i=i+1
  }
}
}


#####PLOTS
xh1=min(outdf$Cases,outdf2$Casesrisk,outdf3$Casesrisk,outdf4$Casesrisk)
xh2=max(outdf$Cases,outdf2$Casesrisk,outdf3$Casesrisk,outdf4$Casesrisk)
xd1=min(outdf$Deaths,outdf2$Deathsrisk,outdf3$Deathsrisk,outdf4$Deathsrisk)
xd2=1.3*max(outdf$Deaths,outdf2$Deathsrisk,outdf3$Deathsrisk,outdf4$Deaths)


nbins=80
outdf %>% 
  select(tau,Reff,Reffrisk) %>% 
  #filter(tau%in%c(1,1.3)) %>% 
  pivot_longer(-tau) %>% 
  ggplot(aes(x=value,fill=name))+
  geom_density(alpha=0.5)+
  geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
  # facet_grid(tau~.,scales = 'free_y')+
  xlab('R Effective')+
  theme_cowplot()+
  scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
#  scale_x_continuous(labels = function(x){return(paste0("10^", x))}) +
  theme(
    axis.text.x = element_markdown(),
    legend.title = element_blank(),
    legend.position=c(0.2,0.8)
    )-> p

p

outdf %>% 
  select(tau,Deaths,Deathsrisk) %>% 
  pivot_longer(-tau) %>% 
  ggplot(aes(x=value,fill=name))+
  geom_density(alpha=0.5)+
  geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
  #  facet_grid(tau~.,scales = 'free_y')+
  xlab('Cumulative deaths')+
  theme_cowplot()+
  scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
  #  scale_x_continuous(labels = function(x){return(paste0("10^", x))}) +
  geom_vline(xintercept = 15208) +
  theme(
    axis.text.x = element_markdown(),
    legend.title = element_blank(),
    legend.position=''
  )+
  xlim(c(xd1,xd2))+
  geom_textbox(aes(x=xd2,y=Inf),label="Moderate severity, reduced VE",
               vjust=1,hjust=0,
               orientation="right-rotated",fill="lightgrey",width = grid::unit(0.8, "npc"))->p2

p2


outdf %>% 
  select(tau,Cases,Casesrisk) %>% 
  pivot_longer(-tau) %>% 
  ggplot(aes(x=value,fill=name))+
  geom_density(alpha=0.5)+
  geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
  #  facet_grid(tau~.,scales = 'free_y')+
  xlab('Cumulative hospitalisations')+
  theme_cowplot()+
  scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
  scale_x_continuous(labels = scales::comma) +
  theme(
    axis.text.x = element_markdown(),
    legend.title = element_blank(),
    legend.position='')+ xlim(c(xh1,xh2)) -> p3


p3

outdf2 %>% 
  select(tau,Reff,Reffrisk) %>% 
  #filter(tau%in%c(1,1.3)) %>% 
  pivot_longer(-tau) %>% 
  ggplot(aes(x=value,fill=name))+
  geom_density(alpha=0.5)+
  geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
  # facet_grid(tau~.,scales = 'free_y')+
  xlab('R Effective')+
  theme_cowplot()+
  scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
  scale_x_continuous(labels = scales::comma) +
  theme(
    axis.text.x = element_markdown(),
    legend.title = element_blank(),
    legend.position=''
  ) -> pb

print(pb)

outdf2 %>% 
  select(tau,Deaths,Deathsrisk) %>% 
  pivot_longer(-tau) %>% 
  ggplot(aes(x=value,fill=name))+
  geom_density(alpha=0.5)+
  geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
  #  facet_grid(tau~.,scales = 'free_y')+
  xlab('Cumulative deaths')+
  theme_cowplot()+
  scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
  scale_x_continuous(labels = scales::comma) +
  geom_vline(xintercept = 15208) +
  theme(
    axis.text.x = element_markdown(),
    legend.title = element_blank(),
    legend.position=''
  )+
  xlim(c(xd1,xd2)) +
  geom_textbox(aes(x=xd2,y=Inf),label="Low severity, reduced VE",
               vjust=1,hjust=0,
               orientation="right-rotated",fill="lightgrey",width = grid::unit(0.9, "npc"))-> pb2

pb2

outdf2 %>% 
  select(tau,Cases,Casesrisk) %>% 
  pivot_longer(-tau) %>% 
  ggplot(aes(x=value,fill=name))+
  geom_density(alpha=0.5)+
  geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
  #  facet_grid(tau~.,scales = 'free_y')+
  xlab('Cumulative hospitalisations')+
  theme_cowplot()+
  scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
  scale_x_continuous(labels = scales::comma) +
  theme(
    axis.text.x = element_markdown(),
    legend.title = element_blank(),
    legend.position='')+xlim(c(xh1,xh2)) -> pb3

pb3


outdf3 %>% 
  select(tau,Reff,Reffrisk) %>% 
  #filter(tau%in%c(1,1.3)) %>% 
  pivot_longer(-tau) %>% 
  ggplot(aes(x=value,fill=name))+
  geom_density(alpha=0.5)+
  geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
  # facet_grid(tau~.,scales = 'free_y')+
  xlab('R Effective')+
  theme_cowplot()+
  scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
  scale_x_continuous(labels = scales::comma) +
  theme(
    axis.text.x = element_markdown(),
    legend.title = element_blank(),
    legend.position=''
  ) -> pc

pc

outdf3 %>% 
  select(tau,Deaths,Deathsrisk) %>% 
  pivot_longer(-tau) %>% 
  ggplot(aes(x=value,fill=name))+
  geom_density(alpha=0.5)+
  geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
  #  facet_grid(tau~.,scales = 'free_y')+
  xlab('Cumulative deaths')+
  geom_vline(xintercept = 15208) +
  theme_cowplot()+
  scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
  scale_x_continuous(labels = scales::comma) +
  theme(
    axis.text.x = element_markdown(),
    legend.title = element_blank(),
    legend.position=''
  )+  xlim(c(xd1,xd2)) +
  geom_textbox(aes(x=xd2,y=Inf),label="Moderate severity, high VE",
               vjust=1,hjust=0,
               orientation="right-rotated",fill="lightgrey",width = grid::unit(0.9, "npc")) -> pc2

pc2 

outdf3 %>% 
  select(tau,Cases,Casesrisk) %>% 
  pivot_longer(-tau) %>% 
  ggplot(aes(x=value,fill=name))+
  geom_density(alpha=0.5)+
  geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
  #  facet_grid(tau~.,scales = 'free_y')+
  xlab('Cumulative hospitalisations')+
  theme_cowplot()+
  scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
  scale_x_continuous(labels = scales::comma) +
  theme(
    axis.text.x = element_markdown(),
    legend.title = element_blank(),
    legend.position='')+ xlim(c(xh1,xh2)) -> pc3


pc3


outdf4 %>% 
  select(tau,Reff,Reffrisk) %>% 
  #filter(tau%in%c(1,1.3)) %>% 
  pivot_longer(-tau) %>% 
  ggplot(aes(x=value,fill=name))+
  geom_density(alpha=0.5)+
  geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
  # facet_grid(tau~.,scales = 'free_y')+
  xlab('R Effective')+
  theme_cowplot()+
  scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
  scale_x_continuous(labels = scales::comma) +
  theme(
    axis.text.x = element_markdown(),
    legend.title = element_blank(),
    legend.position=''
  ) -> pd



outdf4 %>% 
  select(tau,Deaths,Deathsrisk) %>% 
  pivot_longer(-tau) %>% 
  ggplot(aes(x=value,fill=name))+
  geom_density(alpha=0.5)+
  geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
  #  facet_grid(tau~.,scales = 'free_y')+
  xlab('Cumulative deaths')+
  
  theme_cowplot()+
  scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
  scale_x_continuous(labels = scales::comma) +
  geom_vline(xintercept = 15208) +
  theme(
    axis.text.x = element_markdown(),
    legend.title = element_blank(),
    legend.position=''
  )+  xlim(c(xd1,xd2)) +
  geom_textbox(aes(x=xd2,y=Inf),label="Low severity, high VE",
               vjust=1,hjust=0,
               orientation="right-rotated",fill="lightgrey",width = grid::unit(0.9, "npc")) -> pd2


pd2

outdf4 %>% 
  select(tau,Cases,Casesrisk) %>% 
  pivot_longer(-tau) %>% 
  ggplot(aes(x=value,fill=name))+
  geom_density(alpha=0.5)+
  geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
  #  facet_grid(tau~.,scales = 'free_y')+
  xlab('Cumulative hospitalisations')+
  theme_cowplot()+
  scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
  theme(
    axis.text.x = element_markdown(),
    legend.title = element_blank(),
    legend.position='')+ 
      xlim(c(xh1,xh2))  -> pd3


pd3


print((p+p3+p2)/(pb+pb3+pb2)/(pc+pc3+pc2)/(pd+pd3+pd2))+plot_annotation(tag_levels = 'A')->pall
pall


GoldenRatio=(1+sqrt(5))/2

ggsave(pall,width = 8*GoldenRatio,height = 8,filename = paste('figs/','ALSPAC_FIGURE_4PANELS','.png',sep="")) 

