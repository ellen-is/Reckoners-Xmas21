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
escape_severity=c(4,4,1,1)
severity=c(0.8,0.5,0.8,0.5)
numreps=1000

variant=2 #1=delta, 2=omicron
outdf = data.frame(scenario=rep(1:numscenarios,each=NR*numreps),simnum=rep(1:numreps,numscenarios*NR),
                   bootnum=rep(1:NR,numscenarios*numreps),Reff=0,Infections=0,Cases=0,Deaths=0,
                   Reffrisk=0,InfectionsRisk=0,Casesrisk=0,Deathsrisk=0,
                   covidsec=rep(input2$covidsec,each=NR*numreps*numscenarios),
                   CTF=rep(input2$CTF1,each=NR*numreps*numscenarios),
                   LFTsens=rep(input2$covidsec,each=NR*numreps*numscenarios),
                   escape=rep(input2$CTF1,each=NR*numreps*numscenarios),
                   POI=rep("Missing",each=NR*numreps*numscenarios))




i=1

for(SC in 4)
{
  mortality = severity[SC]*mortality2
  for(r1 in 1:numreps)
  {
    cat(paste(r1,"..",sep=""))
    input = input2 #back to default parameters
    if(r1<=(numreps/4))
      {
        input$LFTsensitivity = runif(1,0.8,1.2)*input2$LFTsensitivity
        outdf$POI[outdf$scenario==SC & outdf$simnum==r1] = "LFT sensitivity"
    }
    else if (r1<=(2*numreps/4)) 
      {
        myrand=runif(1,0.8,1.2)
        input$transmissionadvantage = myrand*input2$transmissionadvantage
        input$escape_severity = myrand*input2$escape_severity
        outdf$POI[outdf$scenario==SC & outdf$simnum==r1] = "Vaccine escape"
    }
    else if (r1<=(3*numreps/4)) 
      {
        input$covidsec = runif(1,0.8,1.2)*input2$covidsec
        outdf$POI[outdf$scenario==SC & outdf$simnum==r1] = "COVID security"
    }
    else if (r1<=(4*numreps/4))   
      {
        input$CTF1 = runif(1,0.8,1.2)*input2$CTF1
        outdf$POI[outdf$scenario==SC & outdf$simnum==r1] = "Tracing effectiveness"
    }
    
    
    input[susix]=input2[susix]^input$transmissionadvantage
    input[infix]=input2[infix]^input$transmissionadvantage
    input[deathix]=input2[deathix]^escape_severity[SC]


    outdf$LFT[outdf$simnum==r1] = input$LFTsensitivity
    outdf$escape[outdf$simnum==r1] = input$transmissionadvantage
    outdf$covidsec[outdf$simnum==r1] = input$covidsec
    outdf$CTF[outdf$simnum==r1] = input$CTF1

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
    contactsperson$LFT = myrands < contactsperson$probLFT*input$LFTsensitivity
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



outdf %>%
  filter(POI!="Missing") %>%
  mutate(changeindeaths = Deathsrisk/mean(Deathsrisk)) %>%
  ggplot(aes(POI,Deathsrisk)) + 
  geom_boxplot() + 
  coord_flip() +
  ylab("Cumulative Deaths") + xlab("") -> psens

GoldenRatio=(1+sqrt(5))/2
ggsave(psens,width = 3*GoldenRatio,height = 3,filename = paste('figs/','boxplot_sensitivity','.png',sep="")) 
  




outdf %>%
  filter(POI=="LFT sensitivity") %>%
  ggplot(aes(x=LFT,y=Deathsrisk)) + 
  geom_point(alpha=0.1,color="cadetblue2")+
  ylab("Cumulative Deaths") + xlab("LFT sensitivity") -> ps1

outdf %>%
  filter(POI=="Vaccine escape") %>%
  ggplot(aes(x=escape,y=Deathsrisk)) + 
  geom_point(alpha=0.1,color="coral")+
  ylab("Cumulative Deaths") + xlab("Vaccine escape") -> ps2

ps2

outdf %>%
  filter(POI=="Tracing effectiveness") %>%
  ggplot(aes(x=CTF,y=Deathsrisk)) + 
  geom_point(alpha=0.1,color="aquamarine2")+
  ylab("Cumulative Deaths") + xlab("Tracing effectiveness") -> ps3

outdf %>%
  filter(POI=="COVID security") %>%
  ggplot(aes(x=covidsec,y=Deathsrisk)) + 
  geom_point(alpha=0.1,color="darkorchid")+
  ylab("Cumulative Deaths") + xlab("COVID security") -> ps4



print((ps1+ps2)/(ps3+ps4))+plot_annotation(tag_levels = 'A')->psensall
psensall

ggsave(psensall,width = 5*GoldenRatio,height = 5,filename = paste('figs/','scatter_sensitivity','.png',sep="")) 


outdf %>%
  filter(POI=="LFT sensitivity") %>%
  ggplot(aes(x=LFT,y=Reffrisk)) + 
  geom_point(alpha=0.1,color="cadetblue2")+
  ylab("R") + xlab("LFT sensitivity") -> ps1

outdf %>%
  filter(POI=="Vaccine escape") %>%
  ggplot(aes(x=escape,y=Reffrisk)) + 
  geom_point(alpha=0.1,color="coral")+
  ylab("R") + xlab("Vaccine escape") -> ps2

ps2

outdf %>%
  filter(POI=="Tracing effectiveness") %>%
  ggplot(aes(x=CTF,y=Reffrisk)) + 
  geom_point(alpha=0.1,color="aquamarine2")+
  ylab("R") + xlab("Tracing effectiveness") -> ps3

outdf %>%
  filter(POI=="COVID security") %>%
  ggplot(aes(x=covidsec,y=Reffrisk)) + 
  geom_point(alpha=0.1,color="darkorchid")+
  ylab("R") + xlab("COVID security") -> ps4



print((ps1+ps2)/(ps3+ps4))+plot_annotation(tag_levels = 'A')->psensall
psensall

ggsave(psensall,width = 5*GoldenRatio,height = 5,filename = paste('figs/','scatter_sensitivity_R','.png',sep="")) 



outdf %>%
  group_by(POI,LFT,escape,covidsec,CTF) %>%
  summarise(meanR=mean(Reffrisk),minR=min(Reffrisk),maxR=max(Reffrisk),meanDeaths=mean(Deathsrisk),minDeaths=min(Deathsrisk),maxDeaths=max(Deathsrisk)) %>%
  group_by(POI) %>%
  summarise(popR=mean(meanR),Rlo=min(meanR),Rhi=max(meanR),popD=mean(meanDeaths),Dlo=min(meanDeaths),Dhi=max(meanDeaths)) %>%
  mutate(perchangeR = (Rhi-Rlo)/2,perchangeD = (Dhi-Dlo)/2)
  
