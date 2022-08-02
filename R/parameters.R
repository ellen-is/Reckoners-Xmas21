COMP=seq(0,1,0.1)
myinterval=0.95
NR=100
Nuk=67e6
agestats=data.frame(agegroup=c('0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49',
                               '50-54','55-59','60-64','65-69','70-74','75-79','80+'),
                    N=c(0.058,0.063,0.058,0.055,0.063,0.068,0.068,0.066,0.061,0.068,0.070,0.064,0.054,0.050,0.049,0.033,0.049),
                    IFR=c(0.0,0.0,0.0,0.0001,0.0002,0.0004,0.0006,0.001,0.001,0.002,0.005,0.008,0.017,0.027,0.043,0.062,0.096),
                    Hosp=c(0.001,0.001,0.001,0.002,0.005,0.01,0.016,0.023,0.029,0.039,0.058,0.072,0.102,0.117,0.146,0.177,0.18),
                    HFR=c(0.038,0.038,0.038,0.038,0.038,0.038,0.038,0.04,0.045,0.056,0.078,0.113,0.169,0.232,0.291,0.348,0.535)
)
IFR =c(0.0,0.0,0.0,0.0001,0.00015,0.0002,0.0004,0.0006,0.001,0.001,0.002,0.005,0.008,0.017,0.027,0.043,0.062,0.096,0.096,0.096)
Eng5_2021=c(3359,3716,3695,2010,1340,3612,3947,4074,3951,3739,3724,4085,4048,3501,2997,3016,2255,1545,966,434+117+15)*1000
hosprate = c(0.001,0.001,0.001,0.001,0.002,0.005,0.01,0.016,0.023,0.029,0.039,0.058,0.072,0.102,0.117,0.146,0.177,0.18,0.2,0.2)

mortality = 2*IFR
mortality2 = 2*IFR

popsize = Nuk*c(rep(agestats$N,each=5)/5,rep((1-sum(rep(agestats$N,each=5)/5)),12))
popsize5 = Nuk*agestats$N

input = data.frame(seroprev=0.3,infchild=0.25,covidsec=0.2,CTF1=0.2,
                   ve_inf=0.34,ve_trans=0.45,ve_inf2=0.73,ve_trans2=0.45,
                   ve_i_az=0.30,ve_t_az=0.45,ve_i2_az=0.60,ve_t2_az=0.40,
                   ve_i_pf=0.30,ve_t_pf=0.45,ve_i2_pf=0.80,ve_t2_pf=0.65,
                   ve_i_mod=0.30,ve_t_mod=0.45,ve_i2_mod=0.8,ve_t2_mod=0.65,
                   ve_i_nat=0.50,ve_t_nat=0.65,
                   ve_i_boost=0.6,ve_t_boost=0.95,
                   ve_deathred1=0.85,ve_deathred2=0.96,
                   ve_deathredaz1=0.8,ve_deathredaz2=0.94,
                   ve_deathredpf1=0.85,ve_deathredpf2=0.97,
                   ve_deathredmod1=0.85,ve_deathredmod2=0.97,
                   ve_deathrednat=0.97,ve_deathredvacnat=0.995,
                   ve_deathredboost=0.98,
                   transmissionadvantage=3,SAR=0.03,SARhome=0.103,
                   SARlo=0.028,SARhomelo=0.101,SARhi=0.032,SARhomehi=0.105,
                   maxY=4)


input2 = data.frame(seroprev=0.3,infchild=0.25,covidsec=0.2,CTF1=0.2,
                    ve_inf=0.34,ve_trans=0.45,ve_inf2=0.73,ve_trans2=0.45,
                    ve_i_az=0.30,ve_t_az=0.45,ve_i2_az=0.60,ve_t2_az=0.40,
                    ve_i_pf=0.30,ve_t_pf=0.45,ve_i2_pf=0.80,ve_t2_pf=0.65,
                    ve_i_mod=0.30,ve_t_mod=0.45,ve_i2_mod=0.8,ve_t2_mod=0.65,
                    ve_i_nat=0.50,ve_t_nat=0.65,
                    ve_i_vacnat=0.80,ve_t_vacnat=0.85,
                    ve_i_boost=0.6,ve_t_boost=0.95,
                    ve_deathred1=0.85,ve_deathred2=0.96,
                    ve_deathredaz1=0.8,ve_deathredaz2=0.96,
                    ve_deathredpf1=0.85,ve_deathredpf2=0.97,
                    ve_deathredmod1=0.85,ve_deathredmod2=0.97,
                    ve_deathrednat=0.97,ve_deathredvacnat=0.995,
                    ve_deathredboost=0.98,
                    transmissionadvantage=3,SAR=0.03,SARhome=0.103,
                    SARlo=0.028,SARhomelo=0.101,SARhi=0.032,SARhomehi=0.105,
                    maxY=4)
susix = seq(5,26,2) #indices related to susceptibility
infix = seq(6,26,2) #indices related to infectiousness
deathix = c(27:37) #indices related to mortality rates

Rvalues=seq(7,11,1)
LFTsensitivity=0.5

#Secondary attack rates, household and non household
SAR_notHH=c(0.03,0.087)
SAR_HH=c(0.103,0.158)
SAR_notHH1=c(0.028,0.075)
SAR_HH1=c(0.101,0.143)
SAR_notHH2=c(0.032,0.10)
SAR_HH2=c(0.105,0.175)

#baseline reproduction number with no immunity
R0=7