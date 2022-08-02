library(tidyverse)
library(readxl)
library(boot)
library(fields)
library(cowplot)
library(metR)
library(patchwork)
library(colorspace)
library(ggtext)


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
cc=0; 
x=1




##########RUN FOUR SCENARIOS
numscenarios = 4

escape=c(4,4,4,4)
escape_severity=c(4,4,1,2)
severity=c(0.8,0.5,0.8,0.5)

variant=2 #1=delta, 2=omicron
outdf = data.frame(scenario=rep(1:numscenarios,each=NR*numreps),simnum=rep(1:numreps,numscenarios*NR),
                   bootnum=rep(1:NR,numscenarios*numreps),Reff=0,Infections=0,Cases=0,Deaths=0,
                   Reffrisk=0,InfectionsRisk=0,Casesrisk=0,Deathsrisk=0)


numreps=100

i=1

for(SC in 1:numscenarios)
{
  mortality = severity[SC]*mortality2
  for(r1 in 1:numreps)
  {
    cat(paste(r1,"..",sep=""))
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
    
    degred=degreeperson(contactsperson,input,mult)
    degred = degreeperson_calcR(degred,input,mult)
  
    for(r2 in 1:20){degred = degreeperson_calcR(degred,input,mult)}
  
    myboot=boot(cbind(degred[,"rvsquare"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf$Reff[outdf$scenario==SC & outdf$simnum==r1] = myboot$t
  
    myboot=boot(cbind(degred[,"deaths"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf$Deaths[outdf$scenario==SC & outdf$simnum==r1] = myboot$t
  
    myboot=boot(cbind(degred[,"hospadmissions"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    #myboot=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf$Cases[outdf$scenario==SC & outdf$simnum==r1] = myboot$t
  
    myboot=boot(cbind(degred[,"cases"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    outdf$Infections[outdf$scenario==SC & outdf$simnum==r1] = myboot$t
  
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



# 
# 
# #####PLOTS
xh1=min(outdf$Casesrisk)
xh2=max(outdf$Cases)
xd1=min(outdf$Deathsrisk)
xd2=1.3*max(outdf$Deaths)

nbins=80

Rpanel = function(outdf,plotscenario)
{
  if(plotscenario==1)
  {
    outdf %>%
      select(scenario,Reff,Reffrisk) %>%
      filter(scenario%in%c(plotscenario)) %>%
      pivot_longer(-scenario) %>%
      ggplot(aes(x=value,fill=name))+
      geom_density(alpha=0.5)+
      geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
    #facet_grid(scenario~.,scales = 'free_y')+
      xlab('R Effective')+
      theme_cowplot()+
      scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
    #  scale_x_continuous(labels = function(x){return(paste0("10^", x))}) +
      theme(
        axis.text.x = element_markdown(),
        legend.title = element_blank(),
        legend.position=c(0.2,0.8)
      ) -> p
  }
  else{
    outdf %>%
      select(scenario,Reff,Reffrisk) %>%
      filter(scenario%in%c(plotscenario)) %>%
      pivot_longer(-scenario) %>%
      ggplot(aes(x=value,fill=name))+
      geom_density(alpha=0.5)+
      geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
      #facet_grid(scenario~.,scales = 'free_y')+
      xlab('R Effective')+
      theme_cowplot()+
      scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
      #  scale_x_continuous(labels = function(x){return(paste0("10^", x))}) +
      theme(
        axis.text.x = element_markdown(),
        legend.title = element_blank(),
        legend.position=''
      ) -> p
  }
  
}
hosppanel = function(outdf,plotscenario)
{
  outdf %>%
    select(scenario,Cases,Casesrisk) %>%
    filter(scenario%in%c(plotscenario)) %>%
    pivot_longer(-scenario) %>%
    ggplot(aes(x=value,fill=name))+
    geom_density(alpha=0.5)+
    geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
    xlab('Cumulative hospitalisations')+
    theme_cowplot()+
    scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
    theme(
      axis.text.x = element_markdown(),
      legend.title = element_blank(),
      legend.position='')+ 
    xlim(c(xh1,xh2))  -> pd
  return(pd)
}
Deathspanel = function(outdf,plotscenario,mylabel)
{
  outdf %>%
    select(scenario,Deaths,Deathsrisk) %>%
    filter(scenario%in%c(plotscenario)) %>%
    pivot_longer(-scenario) %>%
    ggplot(aes(x=value,fill=name))+
    geom_density(alpha=0.5)+
    geom_histogram(aes(y=..density..),bins=nbins,alpha=0.3)+
    #facet_grid(scenario~.,scales = 'free_y')+
    xlab('Cumulative deaths')+
    theme_cowplot()+
    scale_fill_discrete_qualitative(palette='Set 2',labels = c("Baseline","Reported risk reduction"))+
    geom_vline(xintercept = 15208) +
    theme(
      axis.text.x = element_markdown(),
      legend.title = element_blank(),
      legend.position=''
    )+
    xlim(c(xd1,xd2))+
    geom_textbox(aes(x=xd2,y=Inf),label=mylabel,
                 vjust=1,hjust=0,
                 orientation="right-rotated",fill="lightgrey",width = grid::unit(0.8, "npc"))->p2
  
  return(p2)
}

p1=Rpanel(outdf,1)
p2=Deathspanel(outdf,1,mylabel="Moderate severity, reduced VE")
p3=hosppanel(outdf,1)

pb=Rpanel(outdf,2)
pb2=Deathspanel(outdf,2,mylabel="Low severity, reduced VE")
pb3=hosppanel(outdf,2)

pc=Rpanel(outdf,3)
pc2=Deathspanel(outdf,3,mylabel="Moderate severity, high VE")
pc3=hosppanel(outdf,3)

pd=Rpanel(outdf,4)
pd2=Deathspanel(outdf,4,mylabel="Low severity, high VE")
pd3=hosppanel(outdf,4)

print((p1+p3+p2)/(pb+pb3+pb2)/(pc+pc3+pc2)/(pd+pd3+pd2))+plot_annotation(tag_levels = 'A')->pall
pall

GoldenRatio=(1+sqrt(5))/2

ggsave(pall,width = 8*GoldenRatio,height = 8,filename = paste('figs/','ALSPAC_FIGURE_4PANELS_2','.png',sep="")) 
 
