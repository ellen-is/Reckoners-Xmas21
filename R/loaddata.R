load("data/riskmatrix.RData")
load("data/vaccinestatus.RData")
load(file="data/contactsperson_withkids6.Rdata")


contactsperson$age[contactsperson$age==-1] = sample(contactsperson$age[contactsperson$age>0],size=sum(contactsperson$age==-1))
contactsperson %>%
  mutate(hassymp = ifelse(age<25,0.25,0.75),comply=0,contacttrace=FALSE) -> contactsperson



contactsperson %>%
  group_by(PID) %>%
  summarise(age=mean(age))%>%
  ungroup() %>%
  mutate(age5 = ifelse(age<=17,1+floor(age/5),ifelse(age<=90,2+floor(age/5),2+floor(90/5)))) %>%
  group_by(age5) %>%
  summarise(sumage=n()) %>%
  mutate(nuk = Eng5_2021[age5]) %>%
  mutate(ppp = nuk/sumage) ->agenumbers