# Set up ----
set.seed(43210) 
library(arules);library(cowplot);library(DataCombine);
library(drc) # for dose response modeling
library(dplyr);library(ellipse);library(expss);library(ggExtra)
library(ggplot2);library(ggpubr);library(ggrepel);library(ggridges)
library(Hmisc);library(modelr);library(patchwork);library(plyr)
library(purrr);library(ggpubr);library(pwr);library(readxl);library(reshape2)
library(rstatix);library(rvest);library(scales)
library(segmented) # for photorepair fluence response model
library(sfsmisc) # for integration under the curve
library(sjstats);library(stats4);library(tibble);library(tidyverse)
library(truncnorm);library(viridis)

# Set working directory ----
setwd("~/Documents/1 QMRA 2020-2022/Data")

# Photorepair fluence response model ----
data_bohrerova2007 <- read_excel("bohrerova-2007.xlsx",sheet=1);head(data_bohrerova2007)
P_LI <- data_bohrerova2007$LI;P_LR <- data_bohrerova2007$LR
P_Fluence <- data_bohrerova2007$photorepair_fluence;lin.mod <- lm(P_LR~0+P_Fluence);summary(lin.mod)
predict.break <- 1500 # mJ/cm2, guess value in Bohrerova and Linden (2007)
Photorepair_model <- segmented(lin.mod,seg.Z = ~P_Fluence, psi=predict.break)
P_fluence_new <- expand.grid(P_Fluence=c(0:8000)) # creates a prediction value data frame
pm <- predict(Photorepair_model, newdata=P_fluence_new,interval="confidence",level=0.95) # predict model values
P_fluence_new$p <- pm[,1] # predicted values of log reactivation
P_fluence_new$pmin <- ifelse(pm[,2]<0,0,pm[,2]) # predicted lower CI, prevents negative values
P_fluence_new$pmax <- pm[,3] # creates CI values # predicted upper CI
P_fluence_new$dif <- P_fluence_new$pmax-P_fluence_new$p # predict ymax - predicted y
P_fluence_new$random <- runif(n=8001,min=P_fluence_new$pmin,max=P_fluence_new$pmax);head(P_fluence_new)
summary(Photorepair_model);Photorepair_model$psi # segmented break-point summary
get_P_fluence <- P_fluence_new$P_Fluence
write.table(x=P_fluence_new,"P_Fluence",sep=",")

photorepair_model <- ggplot()+
  geom_point(data=data_bohrerova2007,aes(x=P_Fluence,y=LR),pch=22,stroke=1.5,size=3)+
  geom_ribbon(data=P_fluence_new,aes(x=P_Fluence,ymin=pmin,ymax=pmax),alpha=0.2)+
  geom_line(data=P_fluence_new,aes(x=P_Fluence,y=p),alpha=0.75)+
  labs(x=expression(paste("Photorepair Fluence (mJ/",cm^2,")")), 
       y=expression(paste(Log[10]," Reactivation (LR)")))+
  scale_y_continuous(limits=c(0,4),breaks=c(0,1,2,3,4))+
  scale_x_continuous(breaks=c(0,1500,3000,4500,6000),labels=c(0,1500,3000,4500,6000),
                     limits=c(0,6000))+
  theme_classic2()+
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24))
photorepair_model
ggsave("PhotorepairFluenceResponse.png",photorepair_model,width=10,height=7)

#  EHEC E. coli Dose Response Data from QMRA Wiki (Weir, 2013) ----
iter <- 10000
exp.dr <- function(kd,dose) 1 - exp(-kd*dose)
kd <- 2.18e-4
kd.upperLim <- 5.99e-4;kd.lowerLim <- 1.2e-4; kd_r <- runif(iter,min=kd.lowerLim,max=kd.upperLim)
write.table(kd_r,"kd",sep="\t")

# Simulate morbidity (not used in current QMRA)
#mb <- runif(n=iter,min=0.05,max=0.1) # probability of illness
#hb <- runif(n=iter,min=0.03,max=0.05) # Health Burden for E. coli STEC, WHO, https://www.who.int/news-room/fact-sheets/detail/e-coli

# Load data for ASTM G-173-03 AM1.5 solar spectrum obtained from nrel.org ----
am1.5 <- expand.grid();am1.5 <- read.csv("AM1_5.csv", header=TRUE)
am1.5$rel <- am1.5$global_tilt/max(am1.5$global_tilt) # relative spectral irradiance
int.280_800 <- integrate.xy(x=am1.5$nm,fx=am1.5$global_tilt,a=280,b=800)
ITot_SSI <- sum(am1.5$global_tilt) # Total solar spectral irradiance, W/m2
int.280_800/ITot_SSI
AM1.5_Int <- approx(x=am1.5$nm,y=am1.5$global_tilt,xout=c(280:800)) #Interpolated AM1.5 Global Tilt Data
AM1.5_Global_Tilt <- data.frame(Wavelength=AM1.5_Int$x,Global_Tilt=AM1.5_Int$y)
AM1.5_Int_280_800 <- subset(AM1.5_Global_Tilt,Wavelength <=800)

# Load Photolyase Absorption Spectra Data ----
Photolyase <- read_excel("Photolyase.xlsx",sheet=1)
# Interpolate photolyase absorption spectra to obtain values at integer wavelengths
Photolyase_Int <- approx(x=Photolyase$Wavelength,y=Photolyase$Relative_Absorption_max,xout=c(250:800),na.rm=TRUE)

# Load Material Transmittance Spectra Data ----
Mat_Int <- read_excel("Transmittance_Spectra_20211103.xlsx",sheet=1)
Spectra_Master <- expand_grid(Wavelength=c(250:800))
Spectra_Master$Photolyase <- Photolyase_Int$y
Spectra_Master$PET_T <- Mat_Int$PET
Spectra_Master$PC_T <- Mat_Int$PC
Spectra_Master$PS_T <- Mat_Int$PS
Spectra_Master$PPCO_T <- Mat_Int$PPCO
Spectra_Master$Gatorade_T <- Mat_Int$Gatorade
Spectra_Master$Aquafina_T <- Mat_Int$Aquafina
Spectra_Master$SevenUp_T <- Mat_Int$SevenUp
Spectra_Master$Ice_T <- Mat_Int$Ice
Spectra_Master$Tritan_T <- Mat_Int$Tritan
Spectra_Master$PET_weighted <- Mat_Int$PET*Spectra_Master$Photolyase
Spectra_Master$PC_weighted <- Mat_Int$PC*Spectra_Master$Photolyase
Spectra_Master$PS_weighted <- Mat_Int$PS*Spectra_Master$Photolyase
Spectra_Master$PPCO_weighted <- Mat_Int$PPCO*Spectra_Master$Photolyase
Spectra_Master$Gatorade_weighted <- Mat_Int$Gatorade*Spectra_Master$Photolyase
Spectra_Master$Aquafina_weighted <- Mat_Int$Aquafina*Spectra_Master$Photolyase
Spectra_Master$SevenUp_weighted <- Mat_Int$SevenUp*Spectra_Master$Photolyase
Spectra_Master$Ice_weighted <- Mat_Int$Ice*Spectra_Master$Photolyase
Spectra_Master$Tritan_weighted <- Mat_Int$Tritan*Spectra_Master$Photolyase
Spectra_Master_Melt <- melt(data=Spectra_Master,id=c("Wavelength"),variable.name="Spectra",value.name = "Value")
Spectra_Master_280_800 <- subset(Spectra_Master,Wavelength >= 280)
Spectra_Master_280_800$Global_Tilt <- AM1.5_Int_280_800$Global_Tilt

Mat_g <- melt(Spectra_Master,id=c("Wavelength"))
Mat_g <- separate(Mat_g,col = variable, into = c("Brand","T"),sep="_")
Mat_g %>%
  ggplot()+
  geom_line(aes(x=Wavelength,y=value,col=Brand,lty=Brand))+
  xlab("Wavelength (nm)")+
  scale_x_continuous(breaks=100*c(2:8))+
  ylab("Transmittance (%)")+
  theme_classic()+
  theme(text = element_text(size=16))

# Calculate photorepair irradiance Global Tilt Irradiance ----
# Global Tilt x Weighted Material Transmittance for each material
# For example PET_P = Global_Tilt x PET_weighted
# Where "Ph_Irr" stands for photorepair irradiance
Ph_Irr <- expand_grid(Wavelength=c(280:800))
Ph_Irr$Global_Tilt <- AM1.5_Int_280_800$Global_Tilt
Ph_Irr$PET_P = Spectra_Master_280_800$Global_Tilt*Spectra_Master_280_800$PET_weighted
Ph_Irr$PC_P = Spectra_Master_280_800$Global_Tilt*Spectra_Master_280_800$PC_weighted
Ph_Irr$PS_P = Spectra_Master_280_800$Global_Tilt*Spectra_Master_280_800$PS_weighted
Ph_Irr$PPCO_P = Spectra_Master_280_800$Global_Tilt*Spectra_Master_280_800$PPCO_weighted
Ph_Irr$Gatorade_P = Spectra_Master_280_800$Global_Tilt*Spectra_Master_280_800$Gatorade_weighted
Ph_Irr$Aquafina_P = Spectra_Master_280_800$Global_Tilt*Spectra_Master_280_800$Aquafina_weighted
Ph_Irr$Ice_P = Spectra_Master_280_800$Global_Tilt*Spectra_Master_280_800$Ice_weighted
Ph_Irr$SevenUp_P = Spectra_Master_280_800$Global_Tilt*Spectra_Master_280_800$SevenUp_weighted
Ph_Irr$Tritan_P = Spectra_Master_280_800$Global_Tilt*Spectra_Master_280_800$Tritan_weighted

ggplot(data=Ph_Irr)+
  geom_line(aes(x=Wavelength,y=Global_Tilt),col="orange",size=1.2,lty=1)+
  geom_line(aes(x=Wavelength,y=PET_P),color="blue",size=1,lty=1)+
  theme_pubr()+
  ylab(expression(paste("Irradiance (µW/",m^2,"/nm)")))+
  xlab("Wavelength (nm)")+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        plot.margin = unit(c(1,1,1,1),units="cm"))+
  annotate(geom="text",x=515,y=0.45,label="Photorepair Irradiance",col="blue",size=4)+
  annotate(geom="text",x=735,y=1.6,label="Spectral Solar Irradiance",col="orange",size=4)
ggsave("Calculating Weighted Photorepair Factor.png",height=7,width=10)
  
# Calculate Weighted Photorepair Factor (WPF) for each material (m) ----
# Total Photorepair Irradiance (TPI) / Total Solar Irradiance (TSI)
TSI <- sum(am1.5$global_tilt)
TSI # W/m2
PF_PET <- integrate.xy(x=Ph_Irr$Wavelength,fx=Ph_Irr$PET_P,a=300,b=500)/TSI
PF_PC <- integrate.xy(x=Ph_Irr$Wavelength,fx=Ph_Irr$PC_P,a=300,b=500)/TSI
PF_PS <- integrate.xy(x=Ph_Irr$Wavelength,fx=Ph_Irr$PS_P,a=300,b=500)/TSI
PF_PPCO <- integrate.xy(x=Ph_Irr$Wavelength,fx=Ph_Irr$PPCO_P,a=300,b=500)/TSI
PF_Tritan <- integrate.xy(x=Ph_Irr$Wavelength,fx=Ph_Irr$Tritan_P,a=300,b=500)/TSI
PF_Gatorade <- integrate.xy(x=Ph_Irr$Wavelength,fx=Ph_Irr$Gatorade_P,a=300,b=500)/TSI/100
PF_Ice <- integrate.xy(x=Ph_Irr$Wavelength,fx=Ph_Irr$Ice_P,a=300,b=500)/TSI/100
PF_SevenUp <- integrate.xy(x=Ph_Irr$Wavelength,fx=Ph_Irr$SevenUp_P,a=300,b=500)/TSI/100
PF_Aquafina <- integrate.xy(x=Ph_Irr$Wavelength,fx=Ph_Irr$Aquafina_P,a=300,b=500)/TSI/100
PF_Values <- data.frame(Plastic = c("PET","PC","PS","PPCO","Tritan","Gatorade","Aquafina","SevenUp","Ice"),
                        PF = c(PF_PET,PF_PC,PF_PS,PF_PPCO,PF_Tritan,PF_Gatorade,PF_Aquafina,PF_SevenUp,PF_Ice))
PF_Values

# Load pveducation solar irradiation data, values in kW/m2 ----
# The conversion is 36000 mJ/cm2 per 1 kWh/m2
PV_data <- read_excel("PV_Data_20211103.xlsx",sheet="Formatted")
plot(x=PV_data$hour_day_1,y=PV_data$insolation_day_1)
plot(x=PV_data$hour_day_90,y=PV_data$insolation_day_90)
plot(x=PV_data$hour_day_180,y=PV_data$insolation_day_180)
PV_melt <- read_excel("PV_Data_20211103.xlsx",sheet="Melted")
PV_melt$Intensity <- factor(PV_melt$Intensity,levels=c("Low","Medium","High"))
gg_PV <- ggplot(data=PV_melt)+
  geom_line(aes(x=Time,y=Insolation,col=Intensity,lty=Intensity))+
  scale_linetype_manual(values=c(4,2,1))+
  theme_pubr()+
  scale_x_continuous(breaks=c(4,6,8,10,12,14,16,18,20),
                     limits=c(4,20))+
  labs(y=expression(paste("Direct Solar Irradiance (kW/",m^2,")")),
       x=expression(paste("Time (24h)")))+
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        axis.title.x = element_text(vjust=0,size=16),
        axis.title.y = element_text(vjust=0,size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))
gg_PV
ggsave("PV Data.png",width=7,height=7)

# Simulating collection times for uniform distributions ----
Times <- expand.grid(c(1:iter))
Times$tc_unif <- runif(iter,4,20)
Times$te_unif <- runif(iter,0,3.5)
Times$kd <- runif(iter,kd.lowerLim,kd.upperLim)
tc <- Times$tc_unif 
te <- Times$te_unif

# Integrating irradiance to obtain fluence ####
Fluence_High <- integrate.xy(x=PV_data$hour_day_180,fx=PV_data$insolation_day_180, 
                             a=tc, b=tc+te,use.spline=TRUE)*36000
Fluence_Medium <- integrate.xy(x=PV_data$hour_day_90,fx=PV_data$insolation_day_90, 
                               a=tc,b=tc+te, use.spline=TRUE)*36000
Fluence_Low <- integrate.xy(x=PV_data$hour_day_1,fx=PV_data$insolation_day_1, 
                            a=tc,b=tc+te,use.spline=TRUE)*36000
Fluence_Low <- ifelse(Fluence_Low<0,0,round(Fluence_Low,1))
Fluence_Medium <- ifelse(Fluence_Medium<0,0,round(Fluence_Medium,1))
Fluence_High <- ifelse(Fluence_High<0,0,round(Fluence_High,1))

Fluence <- data.frame(tc=c(tc),te=c(te),Low=c(Fluence_Low), Medium=c(Fluence_Medium),High=c(Fluence_High))
melt_F <- melt(Fluence, id=c("tc","te"))
head(melt_F)
ggplot(data=melt_F)+
  geom_point(aes(x=tc,y=te,col=value))+
  facet_wrap(~variable)+
  scale_x_continuous(name=expression(paste("Collection Time, ",t[c]," (24h)")),
                     limits=c(4,20),breaks = c(6,9,12,15,18),
                     labels=c("6","9","12","15","18"))+
  scale_y_continuous(name=expression(paste("Exposure Time, ",t[e]," (h)")),limits=c(0,3.5),
                     breaks=c(0,0.5,1,1.5,2,2.5,3,3.5))+
  guides(col=guide_legend(nrow=1,title=expression(paste("Total Direct Fluence (mJ",cm^2,")"))))+
  scale_colour_viridis_c()+
  theme_classic()+
  theme(legend.direction = "horizontal",
        legend.position = "bottom")
ggsave("Total Direct Fluence.png",width=9,height=4)

write.table(Fluence,"Fluence",sep=",")
write.table(PF_Values,"PF",sep = ",")

### Read in previously simulated data #### 2022-05-02 last time done ----
P_Fluence_Mat <- read_excel("Photorepair_Fluence_Mat.xlsx",sheet="Photorepair_Fluence")
head(P_Fluence_Mat)
P_Fluence_LR <- read_excel("Photorepair_Fluence_Mat.xlsx",sheet="LR_Sim")

## UV fluence response data ----
## alter UV dose response data here
LPUVFR <- read_excel("LP-fluence-response-ecoli-spp_20210928.xlsx")
#LPUVFR <- subset(LPUVFR, Pathogen=="Pathogenic")
#LPUVFR <- subset(LPUVFR, Study=="Sommer et al. (2000)")
#LPUVFR <- subset(LPUVFR, Pathogen=="Non-Pathogenic")
#LPUVFR <- subset(LPUVFR, Pathogen=="Non-Pathogenic"|Study=="Zimmer-Thomas et al. (2007)")

# Check if UV fluence response data can be aggregated using Kruskal Wallis test ----

# Kruskal Wallis test for pathogenic vs. non-pathogenic E. coli
kruskal_test(data=LPUVFR,LI~Pathogen)
pairwise.wilcox.test(LPUVFR$LI,LPUVFR$Pathogen,p.adjust.method="BH")
ggline(LPUVFR, x = "Pathogen", y = "LI",add = c("mean_se", "jitter"),ylab = "Log Inactivation", xlab = "Pathogen")
ggboxplot(LPUVFR, x = "Pathogen", y = "LI",color = "Pathogen",xlab = "Pathogen", ylab = "LI")

# Kruskal Wallis test for different strains of E. coli
kruskal_test(data=LPUVFR,LI~Ecoli)
pairwise.wilcox.test(LPUVFR$LI,LPUVFR$Ecoli,p.adjust.method="BH")
ggline(LPUVFR, x = "Ecoli", y = "LI",add = c("mean_se", "jitter"),ylab = "Log Inactivation", xlab = "E. coli Strain")
ggboxplot(LPUVFR, x = "Ecoli", y = "LI",color = "Pathogen",xlab = "E. coli", ylab = "LI")

# Kruskal Wallis test for different UV disinfection studies
pairwise.wilcox.test(LPUVFR$LI,LPUVFR$Study,p.adjust.method="BH")
kruskal_test(data=LPUVFR,LI~Study)
ggline(LPUVFR, x = "Study", y = "LI",add = c("mean_se", "jitter"),ylab = "Log Inactivation", xlab = "Study")
ggboxplot(LPUVFR, x = "Study", y = "LI",color = "Pathogen",xlab = "Study", ylab = "LI")

# Log logistic modeling for UV fluence response ----
CI <- 0.95; newdata <- expand.grid(dose=c(0:40))
LPUVFR_mod <- drm(LI ~ Dose, data=LPUVFR,fct=LL.4(fixed = c(NA,0,NA,NA),names = c("b","c","d","e")))
pm <- predict(LPUVFR_mod,newdata=newdata,interval="confidence",level=CI) # Predicted values
newdata$predLI <- pm[,1]
newdata$predLI_min <- ifelse(pm[,2]<0,0,pm[,2])
newdata$predLI_max <- pm[,3] # Data of CI
A <- summary(LPUVFR_mod); A
coeff <- LPUVFR_mod$coefficients; coeff
b <- coeff[1]
b_se <- A$coefficients[1,2]
c <- 0
d <- coeff[2]
d_se <- A$coefficients[2,2]
e <- coeff[3]
e_se <- A$coefficients[3,2]
LL4 <- function(b,c,d,e,dose) c+(d-c)/(1+exp(b*(log(dose)-log(e))))
dose<-expand.grid(Dose=c(0:40))
dose$LI <-LL4(b,c,d,e,dose$Dose)
LI_ClassB <- LL4(b,c,d,e,16)
print("LI for UV Dose 16")
LI_ClassB # 16 mJ/cm2 is the class B system dose and typical for secondary disinfection
LI_ClassA <- LL4(b,c,d,e,40)
print("LI for UV Dose 40")
LI_ClassA # 40 mJ/cm2 is the class A system dose and typical for primary drinking water disinfection
LI_Fluence8 <- LL4(b,c,d,e,8)
print("LI for UV Dose 8")
LI_Fluence8 # 8 mJ/cm2 is the dose used in Bohrerova and Linden 2007

iter <- 10000
LI_range <- newdata[c(9,17,41),]
LI_8 <- runif(iter,LI_range[1,3],LI_range[1,4])
LI_16 <- runif(iter,LI_range[2,3],LI_range[2,4])
LI_40 <- runif(iter,LI_range[3,3],LI_range[3,4])
#LI_8 <- runif(iter,4.104521,4.937081); LI_16 <- runif(iter,5.278155,6.563605); LI_40 <- runif(iter,5.160219,7.556821)

# POU
POU <- expand.grid(UV8=LI_8)
POU$UV16 <- LI_16
POU$UV40 <- LI_40
write.table(POU,"POU",sep="\t")

# UV fluence response plot ----
gg_UV_dose_response <- 
  ggplot(data=LPUVFR,aes(x=Dose,y=LI))+ 
  geom_ribbon(data=newdata,aes(x=dose,y=predLI,ymin=predLI_min,ymax=predLI_max),alpha=0.30,fill="gray")+
  geom_point(data=LPUVFR,aes(x=Dose,y=LI,shape=Study,col=Ecoli_Strain),size=3)+
  geom_line(data=newdata,aes(x=dose,y=predLI))+
  scale_shape_manual(values=c(0,1,2,5,6,7,9,10,15,16,17,18,19,20,22,23,24,25))+
  scale_x_continuous(limits = c(0,44), expand = c(0,0),breaks=c(0,5,10,15,20,25,30,35,40))+
  scale_y_continuous(limits=c(0,8),breaks=c(0,1,2,3,4,5,6,7,8),expand=c(0,0))+
  labs(x=expression(paste("UV Fluence (mJ/",cm^2,")")),y=expression(paste(Log[10]," Inactivation (LI)")))+
  theme_classic2()+
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        axis.title.x = element_text(vjust=0,size=20),
        axis.title.y = element_text(vjust=0,size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))+
  theme(legend.text=element_text(size=6),
        legend.margin = margin(0.25,0.25,0.25,0.25),
        legend.position="right",
        legend.justification = "center")+
  guides(shape=guide_legend(ncol=2))+
  guides(col=guide_legend(title="Strain"))
gg_UV_dose_response
ggsave("Dose Response UV.png",width=9,height=5)

# Set ingestion volume ----
Vol_ing <- 1

# set breaks for risk values
discrete_labels2 <- c("p ≤ 1E-7","1E-7 < p ≤ 1E-6",
                      "1E-6 < p ≤ 1E-5","1E-5 < p ≤ 1E-4", "1E-4 < p ≤ 1E-3",
                      "1E-3 < p ≤ 1")
discrete_breaks2 <- c(0,1e-7,1e-6,1e-5,1e-4,1e-3,1)
discrete_labels2 <- c("p ≤ 1E-6","1E-6 < p ≤ 1E-4",
                      "1E-4 < p ≤ 1E-2","1E-2 < p ≤ 1")
discrete_breaks2 <- c(0,1e-6,1e-4,1e-2,1)
colours <- c("p ≤ 1E-7"="yellow",
             "1E-7 < p ≤ 1E-6"="forestgreen",
             "1E-6 < p ≤ 1E-5"="darkseagreen",
             "1E-5 < p ≤ 1E-4"="darkcyan",
             "1E-4 < p ≤ 1E-3"="mediumblue",
             "1E-3 < p ≤ 1"="red")
Vol_ing<-1
# Dose response data frame: Variable = Initial Concentration: 100, 1000 and 10000 CFU/L ####
# Constants = Material: PET, Intensity: High, UV Dose = 16 mJ/cm2
DR_conc <- expand.grid(c(1:iter))
DR_conc$tc_unif <- P_Fluence_Mat$tc
DR_conc$te_unif <- P_Fluence_Mat$te
DR_conc$ph_fluence <- P_Fluence_Mat$High_PET
DR_conc$LR <- P_Fluence_LR$High_PET
DR_conc$LI <- LI_16
DR_conc$kd <- kd_r
#DR_conc$Cing_C1 <- log(1,base=10) - DR_conc$LI + DR_conc$LR
#DR_conc$Cing_C10 <- log(10,base=10)-DR_conc$LI+DR_conc$LR
#DR_conc$Cing_C100 <- log(100,base=10)-DR_conc$LI+DR_conc$LR
DR_conc$Cing_C1 <- log(100,base=10) - DR_conc$LI + DR_conc$LR # 100 CFU/L
DR_conc$Cing_C10 <- log(1000,base=10)-DR_conc$LI+DR_conc$LR # 1000 CFU/L
DR_conc$Cing_C100 <- log(10000,base=10)-DR_conc$LI+DR_conc$LR # 10000 CFU/L
DR_conc$DR_C1 <- exp.dr(kd_r,10^(DR_conc$Cing_C1)*Vol_ing)
DR_conc$DR_C10 <- exp.dr(kd_r,10^(DR_conc$Cing_C10)*Vol_ing)
DR_conc$DR_C100 <- exp.dr(kd_r,10^(DR_conc$Cing_C100)*Vol_ing)
DR_conc_melt <- melt(DR_conc,id=c("Var1","tc_unif","te_unif","kd","ph_fluence","LR","LI"))
DR_conc_melt <- separate(data=DR_conc_melt,col="variable",into=c("Step","Conc"))
DR_conc_melt <- subset(DR_conc_melt,Step=="DR")
DR_conc_melt$Conc <- factor(DR_conc_melt$Conc, levels=c("C1","C10","C100"),labels=c("100 CFU/L","1,000 CFU/L","10,000 CFU/L"))
DR_conc_melt$discrete <- discretize(x=DR_conc_melt$value,method="fixed", breaks=c(0,1e-4,1),labels=c("p ≤ 1E-04","p > 1E-04"))
DR_conc_melt$discrete <- factor(x=DR_conc_melt$discrete,levels=c("p ≤ 1E-04","p > 1E-04"))
DR_conc_melt$discrete2 <- discretize(x=DR_conc_melt$value,method="fixed", 
                                     breaks=discrete_breaks2,
                                     labels=discrete_labels2)

## Dose response data frame: Variable = UV Fluence: 8, 16, 40 mJ/cm2 ####
# Constants = Material: PET, Intensity: High, C0 = 100 CFU/L
DR_UV <- expand.grid(c(1:iter))
DR_UV$tc_unif <- P_Fluence_Mat$tc
DR_UV$te_unif <- P_Fluence_Mat$te
DR_UV$ph_fluence <- P_Fluence_Mat$High_PET
DR_UV$LR <- P_Fluence_LR$High_PET
DR_UV$LI16 <- LI_16
DR_UV$LI40 <- LI_40
DR_UV$LI8 <- LI_8
DR_UV$kd <- kd_r
DR_UV$Cing_UV8 <- log(100,base=10)-DR_UV$LI8+DR_UV$LR
DR_UV$Cing_UV16 <- log(100,base=10)-DR_UV$LI16+DR_UV$LR
DR_UV$Cing_UV40 <- log(100,base=10)-DR_UV$LI40+DR_UV$LR
DR_UV$DR_UV8 <- exp.dr(kd_r,10^(DR_UV$Cing_UV8)*Vol_ing)
DR_UV$DR_UV16 <- exp.dr(kd_r,10^(DR_UV$Cing_UV16)*Vol_ing)
DR_UV$DR_UV40 <- exp.dr(kd_r,10^(DR_UV$Cing_UV40)*Vol_ing)
DR_UV_melt <- melt(DR_UV,id=c("Var1","tc_unif","te_unif","LI8","LI16","LI40","ph_fluence","LR","kd"))
DR_UV_melt <- subset(DR_UV_melt,variable=="DR_UV8"|variable=="DR_UV16"|variable=="DR_UV40")
DR_UV_melt <- separate(data=DR_UV_melt,col="variable",into=c("Step","Dose"))
DR_UV_melt$Dose <- factor(DR_UV_melt$Dose,levels = c("UV8","UV16","UV40"),labels=c("8 mJ/cm2","16 mJ/cm2","40 mJ/cm2"))
head(DR_UV_melt)
DR_UV_melt$discrete <- discretize(x=DR_UV_melt$value, method="fixed", breaks=c(0,1e-4,1),labels=c("p ≤ 1E-04","p > 1E-04"))
DR_UV_melt$discrete <- factor(x=DR_UV_melt$discrete,levels=c("p ≤ 1E-04","p > 1E-04"))
DR_UV_melt$discrete2 <- discretize(x=DR_UV_melt$value,method="fixed",
                                   breaks=discrete_breaks2,
                                   labels=discrete_labels2)

## Dose response data frame: Variable = Intensity ####
# Constants = Material: PET, UV Fluence = 16 mJ/cm2, C0 = 100 CFU/L
DR_Intensity <- expand.grid(c(1:iter))
DR_Intensity$tc_unif <- P_Fluence_Mat$tc
DR_Intensity$te_unif <- P_Fluence_Mat$te
DR_Intensity$LI16 <- LI_16

DR_Intensity$ph_fluence_low <- P_Fluence_Mat$Low_PET
DR_Intensity$ph_fluence_medium <- P_Fluence_Mat$Medium_PET
DR_Intensity$ph_fluence_high <- P_Fluence_Mat$High_PET

DR_Intensity$LR_low <- P_Fluence_LR$Low_PET
DR_Intensity$LR_medium <- P_Fluence_LR$Medium_PET
DR_Intensity$LR_high <- P_Fluence_LR$High_PET
DR_Intensity$kd <- kd_r

DR_Intensity$Cing_Low <- log(100,base=10)-DR_Intensity$LI16+DR_Intensity$LR_low
DR_Intensity$Cing_Med <- log(100,base=10)-DR_Intensity$LI16+DR_Intensity$LR_medium
DR_Intensity$Cing_High <- log(100,base=10)-DR_Intensity$LI16+DR_Intensity$LR_high

DR_Intensity$DR_Low <- exp.dr(DR_Intensity$kd,10^(DR_Intensity$Cing_Low)*Vol_ing)
DR_Intensity$DR_Medium <- exp.dr(DR_Intensity$kd,10^(DR_Intensity$Cing_Med)*Vol_ing)
DR_Intensity$DR_High <- exp.dr(DR_Intensity$kd,10^(DR_Intensity$Cing_High)*Vol_ing)

DR_Intensity_melt <- melt(DR_Intensity,id=c("Var1","tc_unif","te_unif",
                                            "ph_fluence_high","ph_fluence_medium",
                                            "ph_fluence_low", "LI16","LR_medium","LR_low","LR_high",
                                            "Cing_Low","Cing_Med","Cing_High","kd"))
head(DR_Intensity_melt)
DR_Intensity_melt <- separate(data=DR_Intensity_melt,col="variable",into=c("Step","Intensity"))
DR_Intensity_melt$Intensity <- factor(x=DR_Intensity_melt$Intensity,levels=c("Low","Medium","High"))
head(DR_Intensity_melt)
hist(DR_Intensity_melt$value)
DR_Intensity_melt$discrete <- discretize(x=DR_Intensity_melt$value,method="fixed",breaks=c(0,1e-4,1),labels=c("p ≤ 1E-04","p > 1E-04"))
DR_Intensity_melt$discrete <- factor(x=DR_Intensity_melt$discrete,levels=c("p ≤ 1E-04","p > 1E-04"))
DR_Intensity_melt$discrete2 <- discretize(x=DR_Intensity_melt$value,method="fixed",
                                          breaks=discrete_breaks2,
                                          labels=discrete_labels2)

## Dose response data frame: Variable = Material ####
# Constants = UV Dose = 16 mJ/cm2, Intensity: High, C0 = 100 CFU/L
DR_Material <- expand.grid(c(1:iter))
DR_Material$tc_unif <- P_Fluence_Mat$tc
DR_Material$te_unif <- P_Fluence_Mat$te
DR_Material$LR_PET <- P_Fluence_LR$High_PET
DR_Material$LR_PS <- P_Fluence_LR$High_PS
DR_Material$LR_PC <- P_Fluence_LR$High_PC
DR_Material$LR_PPCO <- P_Fluence_LR$High_PPCO
DR_Material$LR_Tritan <- P_Fluence_LR$High_Tritan
DR_Material$LR_Aquafina <- P_Fluence_LR$High_Aquafina
DR_Material$LR_SevenUp <- P_Fluence_LR$High_SevenUp
DR_Material$LR_Gatorade <- P_Fluence_LR$High_Gatorade
DR_Material$LR_Ice <- P_Fluence_LR$High_Ice
DR_Material$kd <- kd_r
DR_Material$Cing_PET <- log(100,base=10)-LI_16+DR_Material$LR_PET
DR_Material$Cing_PC <- log(100,base=10)-LI_16+DR_Material$LR_PC
DR_Material$Cing_PS <- log(100,base=10)-LI_16+DR_Material$LR_PS
DR_Material$Cing_PPCO <- log(100,base=10)-LI_16+DR_Material$LR_PPCO
DR_Material$Cing_Tritan <- log(100,base=10)-LI_16+DR_Material$LR_Tritan
DR_Material$Cing_Aquafina <- log(100,base=10)-LI_16+DR_Material$LR_Aquafina
DR_Material$Cing_SevenUp <- log(100,base=10)-LI_16+DR_Material$LR_SevenUp
DR_Material$Cing_Gatorade <- log(100,base=10)-LI_16+DR_Material$LR_Gatorade
DR_Material$Cing_Ice <- log(100,base=10)-LI_16+DR_Material$LR_Ice

DR_Material$DR_PET <- exp.dr(kd_r,10^(DR_Material$Cing_PET)*Vol_ing)
DR_Material$DR_PS <- exp.dr(kd_r,10^(DR_Material$Cing_PS)*Vol_ing)
DR_Material$DR_PC <- exp.dr(kd_r,10^(DR_Material$Cing_PC)*Vol_ing)
DR_Material$DR_PPCO <- exp.dr(kd_r,10^(DR_Material$Cing_PPCO)*Vol_ing)
DR_Material$DR_Tritan <- exp.dr(kd_r,10^(DR_Material$Cing_Tritan)*Vol_ing)
DR_Material$DR_Gatorade <- exp.dr(kd_r,10^(DR_Material$Cing_Gatorade)*Vol_ing)
DR_Material$DR_Ice <- exp.dr(kd_r,10^(DR_Material$Cing_Ice)*Vol_ing)
DR_Material$DR_Aquafina <- exp.dr(kd_r,10^(DR_Material$Cing_Aquafina)*Vol_ing)
DR_Material$DR_SevenUp <- exp.dr(kd_r,10^(DR_Material$Cing_SevenUp)*Vol_ing)
DR_Material_melt <- melt(DR_Material,id=c("Var1","tc_unif","te_unif","kd"))
DR_Material_melt <- separate(data=DR_Material_melt,col="variable",into=c("Step","Material"))
DR_Material_melt <- subset(DR_Material_melt,Step=="DR")
DR_Material_melt$Material <- factor(x=DR_Material_melt$Material,levels=c("PET","PS","PC","PPCO","Tritan","Aquafina","Gatorade","SevenUp","Ice"))
DR_Material_melt$discrete <- discretize(x=DR_Material_melt$value,method="fixed",breaks=c(0,1e-4,1),labels=c("p ≤ 1E-04","p > 1E-04"))
DR_Material_melt$discrete <- factor(x=DR_Material_melt$discrete,levels=c("p ≤ 1E-04","p > 1E-04"))
DR_Material_melt$discrete2 <- discretize(x=DR_Material_melt$value, method="fixed",
                                         breaks=discrete_breaks2,
                                         labels=discrete_labels2)

# plots with values more discretized ----  
gg_DR_Vary_UV_2 <- 
  ggplot(data=DR_UV_melt)+
  geom_point(aes(x=tc_unif,y=te_unif,colour=log10(value)))+
  scale_color_gradient2(midpoint = -4 , low = "white", mid = "steelblue",high = "red",space="Lab")+
  geom_point(data=subset(DR_UV_melt,log10(value)>-4),aes(x=tc_unif,y=te_unif),color="red")+
  facet_grid(~Dose)+
  title(main="A) Concentration")+
  scale_x_continuous(name=expression(paste("Collection Time, ",t[c]," (24h)")),
                     limits=c(4,20),breaks = c(6,9,12,15,18),labels=c("6","9","12","15","18"))+
  scale_y_continuous(name=expression(paste("Exposure Time, ",t[e]," (h)")),limits=c(0,3.5),
                     breaks=c(0,0.5,1,1.5,2,2.5,3,3.5))+
  theme_classic(base_size = 20)+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.text = element_text(size=24),
        legend.title = element_text(size=24),
        strip.text = element_text(size=20),
        plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "horizontal",legend.position = "top")+
  guides(color=guide_legend(title=expression(paste("Risk of Infection (",Log[10],10^X,")")),override.aes = list(size = 5)))

gg_DR_Vary_Intensity_2 <- 
  ggplot(data=DR_Intensity_melt)+
  geom_point(aes(x=tc_unif,y=te_unif,colour=log10(value)))+
  scale_color_gradient2(midpoint = -4 , low = "white", mid = "steelblue",high = "red",space="Lab")+
  geom_point(data=subset(DR_Intensity_melt,log10(value)>-4),aes(x=tc_unif,y=te_unif),color="red")+
  facet_grid(~Intensity)+
  scale_x_continuous(name=expression(paste("Collection Time, ",t[c]," (24h)")),
                     limits=c(4,20),breaks = c(6,9,12,15,18),labels=c("6","9","12","15","18"))+
  scale_y_continuous(name=expression(paste("Exposure Time, ",t[e]," (h)")),limits=c(0,3.5),
                     breaks=c(0,0.5,1,1.5,2,2.5,3,3.5))+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        strip.text = element_text(size=20),
        plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "horizontal",legend.position = "top") +
  guides(color=guide_legend(title=expression(paste("Risk of Infection (",Log[10],10^X,")")),override.aes = list(size = 5)))

gg_DR_Vary_Material_2 <- 
  ggplot(data=DR_Material_melt)+
  geom_point(aes(x=tc_unif,y=te_unif,colour=log10(value)))+
  scale_color_gradient2(midpoint = -4 , low = "white", mid = "steelblue",high = "red",space="Lab")+
  geom_point(data=subset(DR_Material_melt,log10(value)>-4),aes(x=tc_unif,y=te_unif),color="red")+
  facet_wrap(~Material,nrow=3,ncol=3)+
  scale_x_continuous(name=expression(paste("Collection Time, ",t[c]," (24h)")),
                     limits=c(4,20),breaks = c(6,9,12,15,18),labels=c("6","9","12","15","18"))+
  scale_y_continuous(name=expression(paste("Exposure Time, ",t[e]," (h)")),limits=c(0,3.5),
                     breaks=c(0,0.5,1,1.5,2,2.5,3,3.5))+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        strip.text = element_text(size=20),
        plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "vertical",legend.position = "none")+
  guides(color=guide_legend(title=expression(paste("Risk of Infection (",Log[10],10^X,")")),override.aes = list(size = 5)))
  
gg_DR_Vary_Conc_2 <- 
  ggplot(data=DR_conc_melt)+
  geom_point(aes(x=tc_unif,y=te_unif,colour=log10(value)))+
  scale_color_gradient2(midpoint = -4 , low = "white", mid = "steelblue",high = "red",space="Lab")+
  geom_point(data=subset(DR_conc_melt,log10(value)>-4),aes(x=tc_unif,y=te_unif),color="red")+
  facet_grid(~Conc)+
  scale_x_continuous(name=expression(paste("Collection Time, ",t[c]," (24h)")),
                     limits=c(4,20),breaks = c(6,9,12,15,18),labels=c("6","9","12","15","18"))+
  scale_y_continuous(name=expression(paste("Exposure Time, ",t[e]," (h)")),limits=c(0,3.5),
                     breaks=c(0,0.5,1,1.5,2,2.5,3,3.5))+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.text = element_text(size=28),
        legend.title = element_text(size=28),
        strip.text = element_text(size=20),
        plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "horizontal",legend.position = "top")+
  guides(color=guide_legend(title=expression(paste("Risk of Infection (",Log[10],10^X,")")),override.aes = list(size = 8)))

gg_combine_QMRA2 <- ggarrange(gg_DR_Vary_Conc_2,
                              gg_DR_Vary_Intensity_2,
                              gg_DR_Vary_UV_2,
                              gg_DR_Vary_Material_2,
                              nrow=4,legend="top",
                              common.legend = TRUE,
                              heights=c(1,1,1,2.4),
                              labels = c("A) Variable: Concentration",
                                        "B) Variable: Solar Intensity",
                                        "C) Variable: UV Fluence",
                                        "D) Variable: Material"),
                              font.label=list(size=20),
                              hjust=-0.2,
                              vjust=1.5)

ggsave("20220517 Figure 3 Combined.png",gg_combine_QMRA2,width=12,height=24)

# Violin Plots
# Variable UV
violin_UV <- 
  ggplot(data=DR_UV_melt,aes(x=Dose,y=log10(value)))+
  geom_violin(aes(x=Dose,y=log10(value),color=Dose),size=2)+
  scale_color_manual(values=c("plum3","purple","purple4"))+
  scale_y_continuous(limits=c(-10,-2),breaks=c(-10,-8,-6,-4,-2,0))+
  geom_hline(yintercept =-4,lty="solid",color="orange",size=2)+
  stat_summary(fun=median, geom="point", size=4, color="red",shape=3) +
  stat_summary(fun=mean, geom="point", size=4, color="black",shape=15) +
  theme_pubr(base_size = 20)+
  ylab(expression(paste(Log[10]," Risk of Infection")))+
  xlab(expression(paste("LP UV Fluence (mJ/",cm^2,")")))+
  scale_x_discrete(labels=c("8","16","40"))+
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "horizontal",legend.position = "none")
violin_UV;ggsave("Violin_UV_0628.png",violin_UV,width=6,height=6)

violin_conc <- 
  ggplot(data=DR_conc_melt,aes(x=Conc,y=log10(value)))+
  geom_violin(aes(x=Conc,y=log10(value),color=Conc),size=2)+
  scale_y_continuous(limits=c(-10,0),breaks=c(-10,-8,-6,-4,-2,0))+
  scale_color_manual(values=c("red","red2","red4"))+
  ylab(expression(paste(Log[10]," Risk of Infection")))+
  stat_summary(fun=median, geom="point", size=4, color="red",shape=3) +
  stat_summary(fun=mean, geom="point", size=4, color="black",shape=15) +
  scale_x_discrete(labels=c("100","1,000","10,000"))+
  xlab("Concentration (CFU/L)")+
  geom_hline(yintercept =-4,lty="solid",color="orange",size=2)+
  theme_pubr(base_size = 20)+
  ylab(expression(paste(Log[10]," Risk of Infection")))+
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "horizontal",legend.position = "none")
violin_conc;ggsave("Violin_Conc_0628.png",violin_conc,width=6,height=6)

violin_intensity <- 
  ggplot(data=DR_Intensity_melt,aes(x=Intensity,y=log10(value)))+
  geom_violin(aes(x=Intensity,y=log10(value),color=Intensity),size=2)+
  scale_y_continuous(limits=c(-10,-2),breaks=c(-10,-8,-6,-4,-2,0))+
  scale_color_manual(values=c("yellow2","orange2","orange4"))+
  stat_summary(fun=median, geom="point", size=4, color="red",shape=3) +
  stat_summary(fun=mean, geom="point", size=4, color="black",shape=15) +
  ylab(expression(paste(Log[10]," Risk of Infection")))+
  xlab("Sunlight Intensity")+
  geom_hline(yintercept =-4,lty="solid",color="orange",size=2)+
  theme_pubr(base_size = 20)+
  ylab(expression(paste(Log[10]," Risk of Infection")))+
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "horizontal",legend.position = "none")
violin_intensity;ggsave("Violin_Intensity_0628.png",violin_intensity,width=6,height=6)


violin_material <- 
  ggplot(data=DR_Material_melt,aes(x=Material,y=log10(value)))+
  #geom_violin(aes(x=0,y=log10(value),color=Material),size=2)+
  geom_violin(aes(x=Material,y=log10(value),color=Material),size=2)+
  stat_summary(fun=median, geom="point", size=4, color="red",shape=3) +
  stat_summary(fun=mean, geom="point", size=4, color="black",shape=15) +
  scale_y_continuous(limits=c(-10,-2),breaks=c(-10,-8,-6,-4,-2,0))+
  #facet_wrap(facets=vars(Material),nrow=3,ncol=3)+
  ylab(expression(paste(Log[10]," Risk of Infection")))+
  xlab("Container Material")+
  geom_hline(yintercept =-4,lty="solid",color="orange",size=2)+
  theme_pubr(base_size = 16)+
  ylab(expression(paste(Log[10]," Risk of Infection")))+
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
#        axis.ticks.x = element_blank(),
        axis.title.x=element_text(size=20),
        axis.text.x = element_text(angle=30,hjust=1),
        legend.direction = "horizontal",legend.position = "none")
violin_material;ggsave("Violin_Material.png",violin_material,width=12,height=6)

violin_combine <- ggarrange(violin_conc,violin_intensity,violin_UV,violin_material,nrow=4,ncol=1,align=c("hv"),
                            labels = c("A) Variable: Concentration",
                                       "B) Variable: Solar Intensity",
                                       "C) Variable: UV Fluence",
                                       "D) Variable: Material"),font.label=c(size=20),
                            hjust=-0.2,vjust=1.2)
ggsave("Violin_Combined_20220629.png",height=23,width=7)

#### Violin plots with WHO 10E-6 ####

# Violin Plots
# Variable UV
violinWHO_UV <- 
  ggplot(data=DR_UV_melt)+
  geom_violin(aes(x=Dose,y=log10(value*mb*hb),color=Dose),size=2)+
  scale_y_continuous(limits=c(-10,-2),breaks=c(-10,-8,-6,-4,-2,0))+
  geom_hline(yintercept =-6,lty="solid",color="purple",size=2)+
  theme_pubr(base_size = 20)+
  ylab(expression(paste(Log[10]," Risk of Infection")))+
  xlab(expression(paste("LP UV Dose (mJ/",cm^2,")")))+
  scale_x_discrete(labels=c("8","16","40"))+
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "horizontal",legend.position = "none")
violinWHO_UV;ggsave("violinWHO_UV.png",violinWHO_UV,width=6,height=6)

violinWHO_conc <- 
  ggplot(data=DR_conc_melt)+
  geom_violin(aes(x=Conc,y=log10(value*mb*hb),color=Conc),size=2)+
  scale_y_continuous(limits=c(-10,-2),breaks=c(-10,-8,-6,-4,-2,0))+
  ylab(expression(paste(Log[10]," Risk of Infection")))+
  scale_x_discrete(labels=c("100","1,000","10,000"))+
  xlab("Initial Concentration (organisms/L)")+
  geom_hline(yintercept =-6,lty="solid",color="purple",size=2)+
  theme_pubr(base_size = 20)+
  ylab(expression(paste(Log[10]," Risk of Infection")))+
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "horizontal",legend.position = "none")
violinWHO_conc
ggsave("violinWHO_Conc.png",violinWHO_conc,width=6,height=6)


violinWHO_intensity <- 
  ggplot(data=DR_Intensity_melt)+
  geom_violin(aes(x=Intensity,y=log10(value*mb*hb),color=Intensity),size=2)+
  scale_y_continuous(limits=c(-10,-2),breaks=c(-10,-8,-6,-4,-2,0))+
  ylab(expression(paste(Log[10]," Risk of Infection")))+
  xlab("Sunlight Intensity")+
  geom_hline(yintercept =-6,lty="solid",color="purple",size=2)+
  theme_pubr(base_size = 20)+
  ylab(expression(paste(Log[10]," Risk of Infection")))+
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "horizontal",legend.position = "none")
violinWHO_intensity
ggsave("violinWHO_Intensity.png",violinWHO_intensity,width=6,height=6)


violinWHO_material <- 
  ggplot(data=DR_Material_melt)+
  #geom_violinWHO(aes(x=0,y=log10(value*mb*hb),color=Material),size=2)+
  geom_violin(aes(x=Material,y=log10(value*mb*hb),color=Material),size=2)+
  scale_y_continuous(limits=c(-10,-2),breaks=c(-10,-8,-6,-4,-2,0))+
  #facet_wrap(facets=vars(Material),nrow=3,ncol=3)+
  ylab(expression(paste(Log[10]," Risk of Infection")))+
  xlab("Container Material")+
  geom_hline(yintercept =-6,lty="solid",color="purple",size=2)+
  theme_pubr(base_size = 16)+
  ylab(expression(paste(Log[10]," Risk of Infection")))+
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.direction = "horizontal",legend.position = "none")
violinWHO_material
ggsave("violinWHO_Material.png",violinWHO_material,width=12,height=6)

# CDF plots ----
cdf_UV <- ggplot(data=DR_UV_melt)+
  stat_ecdf(aes(x=value,col=Dose,lty=Dose),geom = "step")+
  ylab("CDF")+
  scale_x_log10(name=expression(paste("Risk of Infection")),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits=c(1e-10,1e-2))+
  theme_pubr()+
  theme_pubr(base_size = 12)+
  theme(legend.title = element_blank(),
        plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "vertical",
        legend.position = c(0.2,0.8))
cdf_UV

cdf_C <- ggplot(data=DR_conc_melt)+
  stat_ecdf(aes(x=value,col=Conc,lty=Conc),geom = "step")+
  ylab("CDF")+
  scale_x_log10(name=expression(paste("Risk of Infection")),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits=c(1e-10,1e-2))+
  theme_pubr()+
  theme_pubr(base_size = 12)+
  theme(legend.title = element_blank(),
        plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "vertical",
        legend.position = c(0.2,0.8))

cdf_M <- ggplot(data=DR_Material_melt)+
  stat_ecdf(aes(x=value,col=Material,lty=Material),geom = "step")+
  xlab(expression(paste(log[10]," Risk")))+
  ylab("CDF")+
  scale_x_log10(name=expression(paste("Risk of Infection")),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits=c(1e-10,1e-2))+
  theme_pubr()+
  theme_pubr(base_size = 12)+
  theme(legend.title = element_blank(),
        plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "vertical",
        legend.position = c(0.2,0.8))

cdf_I <- ggplot(data=DR_Intensity_melt)+
  stat_ecdf(geom = "step",aes(x=value,col=Intensity,lty=Intensity))+
  xlab(expression(paste(log[10]," Risk")))+
  ylab("CDF")+
  scale_linetype_manual(values=c(1,3,4))+
  scale_x_log10(name=expression(paste("Risk of Infection")),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits=c(1e-10,1e-2))+
  theme_pubr(base_size = 12)+
  theme(legend.title = element_blank(),
        plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "vertical",
        legend.position = c(0.2,0.8))
cdf_I

cdf_combine <- ggarrange(cdf_C,cdf_I,cdf_UV,cdf_M,nrow=4,ncol=1,align=c("hv"),
                         labels = c("A) Variable: Concentration",
                                    "B) Variable: Solar Intensity",
                                    "C) Variable: UV Fluence",
                                    "D) Variable: Material"),font.label=c(size=20),
                         hjust=-0.2,vjust=1.2)

ggsave("Pi Plot Vary UV 2.png",plot=gg_DR_Vary_UV_2,width=18,height=6)
ggsave("Pi Plot Vary Intensity 2.png",plot=gg_DR_Vary_Intensity_2,width=18,height=6)
ggsave("Pi Plot Vary Material 2.png",plot=gg_DR_Vary_Material_2, width=18,height=10.5) 
ggsave("Pi Plot Vary Conc 2.png",plot=gg_DR_Vary_Conc_2, width=18,height=6) 

ggsave("Pi Boxplot-Material-Vary.png",plot=box_M, width=6,height=6) 
ggsave("Pi Boxplot-Intensity-Vary.png",plot=box_I,width=6,height=6) 
ggsave("Pi Boxplot-Conc-Vary.png",plot=box_C,width=6,height=6) 
ggsave("Pi Boxplot-Dose-Vary.png",plot=box_UV,width=6,height=6)
ggsave("Pi Boxplot-Combined.png",plot=box_combine,width=12,height=12)


ggsave("CDF-UV-Vary.png",plot=cdf_UV,width=6,height=6)
ggsave("CDF-Conc-Vary.png",plot=cdf_C,width=6,height=6)
ggsave("CDF-Material-Vary.png",plot=cdf_M,width=6,height=6)
ggsave("CDF-intensity-Vary.png",plot=cdf_I,width=6,height=6)
ggsave("20220517 Figure S6 CDF-Combined.png",plot=cdf_combine,width=6,height=18)


# Risk Assessment for Culturally Adapted Water Collection Pattern ####
# Variables: High Intensity, UV Doe = 16 mJ/cm2,
# Material = PET, Concentration = 100 CFU/L
cadp <- expand.grid(Iter=c(1:iter))
y1 = rtruncnorm(n=iter,a=7,b=Inf,mean=9,sd=1.5) # morning collection peak
y2 = rnorm(iter, 16, 1.5) # afternoon collection peak
peak1 <- 0.4 # weight of morning peak to account for greater afternoon collection (0.6)
w <- rbinom(iter, 1, peak1) # Binomial distribution for peaks 1 and 2
cadp$tc <- w*y1 + (1-w)*y2 # final collection time distribution
cadp$walkspeed <- runif(iter,3.5,4) # walking speed, km/h, Bimla et al. (2003)
dist_mean <- 0.243 # km, rural, used water source, Boone et al. (2011)
dist_sd <- 0.582 # km, rural, used water source, Boone et al. (2011)
cadp$distance <- rtruncnorm(n=iter,a=0,b=Inf,
                            mean=dist_mean,sd=dist_sd) # distance distribution, prevent negative values
cadp$te <- cadp$distance/cadp$walkspeed # walking time distribution, hours

# Perform fluence calculation and simulations
cadp$Fluence <- integrate.xy(x=PV_data$hour_day_180,fx=PV_data$insolation_day_180, 
                              a=cadp$tc, b=cadp$tc+cadp$te,use.spline=TRUE)*36000
cadp$PFluence_PET <- round(cadp$Fluence * PF_Values$PF[1],0) 
cadp$LI16 <- LI_16
cadp$kd <- kd_r
write.table(cadp,"cadp",sep = ",")
cadp <- read_excel("Final_Photorepair_Simulations.xlsx",sheet="Culturally_Adapted")
cadp$Cing_LowC <- log(100,base=10) - cadp$LI16 + cadp$LR_hypot
cadp$Cing_MedC <- log(1000,base=10) - cadp$LI16 + cadp$LR_hypot
cadp$Cing_HighC <- log(10000,base=10) - cadp$LI16 + cadp$LR_hypot
cadp$DR_HighC <- exp.dr(cadp$kd,10^(cadp$Cing_HighC)*Vol_ing)
cadp$DR_MedC <- exp.dr(cadp$kd,10^(cadp$Cing_MedC)*Vol_ing)
cadp$DR_LowC <- exp.dr(cadp$kd,10^(cadp$Cing_LowC)*Vol_ing)
cadp <- melt(cadp,id=c("Iter","tc","walkspeed","distance","te","Fluence",
                         "PFluence_PET","LI16","kd","LR_hypot"))
cadp <- separate(data=cadp,col=variable,into = c("Step","Concentration"),sep = "_") 
cadp <- subset(cadp,Step == "DR")
cadp$discrete <- discretize(x=cadp$value,method="fixed",breaks=c(0,1e-4,1),labels=c("p ≤ 1E-04","p > 1E-04"))
cadp$discrete <- factor(x=cadp$discrete,levels=c("p ≤ 1E-04","p > 1E-04"))
cadp$discrete2 <- discretize(x=cadp$value, method="fixed",
                                         breaks=discrete_breaks2,
                                         labels=discrete_labels2)
cadp$Concentration <- factor(cadp$Concentration,levels=c("LowC","MedC","HighC"),
                              labels=c("100 CFU/L","1,000 CFU/L","10,000 CFU/L"))
gg_DR_cadp <- 
  ggplot(data=cadp)+
  geom_point(aes(x=tc,y=te,colour=log10(value)))+
  scale_color_gradient2(midpoint = -4 , low = "white", mid = "steelblue",high = "red",space="Lab")+
  geom_point(data=subset(cadp,log10(value)>-4),aes(x=tc,y=te),color="red")+
  facet_grid(~Concentration)+
  scale_x_continuous(name=expression(paste("Collection Time, ",t[c]," (24h)")),
                     limits=c(6,20),breaks = c(6,9,12,15,18),labels=c("6","9","12","15","18"))+
  scale_y_continuous(name=expression(paste("Exposure Time, ",t[e]," (h)")),limits=c(0,1),
                     breaks=c(0,0.25,0.5,0.75,1))+
  theme_classic()+          
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),  
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        strip.text = element_text(size=20),
        plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "horizontal",legend.position = "top")+
  guides(color=guide_legend(title="Risk of Infection",nrow=1,
         override.aes = list(size = 5)))
ggsave("Culturally Adapted 20220603.png",gg_DR_cadp,width=10,height=6) 

cadp_violin <- ggplot(data=cadp)+
  geom_violin(aes(x=Concentration,y=log10(value),color=Concentration),size=2)+
  scale_y_continuous(name=expression(paste(Log[10]," Risk of Infection")),
                     limits=c(-10,-2),
                     breaks=c(-10,-8,-6,-4,-2,0))+
  geom_hline(yintercept = -4,lty="solid",color="orange",size=2)+
  scale_x_discrete(labels=c("100","1,000","10,000"))+
  xlab("Concentration (CFU/L)")+
  theme_pubr(base_size = 24)+
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "horizontal",legend.position = "none")
cadp_violin
ggsave("Violin_CulturallyAdapted.png",cadp_violin,width=6,height=6)

cdf_H <- ggplot(data=cadp)+
  stat_ecdf(geom = "step",aes(x=value,col=Concentration,lty=Concentration))+
  ylab("CDF")+
  theme_pubr(base_size = 16)+
  scale_x_log10(name="Risk of Infection",
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits=c(1e-10,1e-2))+
  theme(legend.title = element_blank(),
        plot.margin = unit(c(1,1,1,1),"cm"),
        legend.direction = "vertical",
        legend.position = c(0.8,0.28))
cdf_H

ggsave("20220517 Hypot Box-Plot.png",plot=box_H,width=5,height=7)
ggsave("20220517 Hypot CDF.png",plot=cdf_H,width=7,height=7)
       
gg_Hy <- ggarrange(gg_DR_Hypot,box_H,cdf_H,
                 labels=c("A)","B)","C)"))