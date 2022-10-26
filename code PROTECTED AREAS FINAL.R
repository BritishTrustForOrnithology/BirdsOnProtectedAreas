
#############################################################################

### Code for running Barnes et al. 2022
### Do conservation designations provide positive benefits for bird species and communities?
### A. E. Barnes1, J. G. Davies2, B. Martay2, P.H. Boersch-Supan1, S.J. Harris1, D.G. Noble1, J.W. Pearce-Higgins1, R.A. Robinson1,*
  
### Section 1: Species-specific models of population metric ~ cover_of_protected area
###           (a) standard models
###           (b) matched models
###           This section demonstrates the methods used but can't run without raw data
###           Raw data is available on request from BTO
###           
### Section 2: Simple summaries
###           (a) numbers of species with positive and negative responses to PAs
###           (b) mean responses across all species
###           This section can be run using the file 'Suppl File 1.xlsx'
###           from: https://doi.org/10.6084/m9.figshare.20200895 
###
### Section 3: Traits analysis
###           (a) Response to PA ~ conservation status
###           (b) Response to PA ~ ecological traits
###           This section can be run using the file 'Suppl File 1.xlsx'
###           from: https://doi.org/10.6084/m9.figshare.20200895 
#############################################################################

library(readxl)
library(mgcv)
library(MCMCglmm)
library(phytools)
library(arm)
library(MatchIt)

# set drive to read your datafile
setwd("C:/Users/admin/OneDrive - British Trust for Ornithology/Protected Areas/")


##########   Section 1a
# Structure of species-specific models identifying the relationship
# between range and population measures (occupancy, colonisation, persistence,
# abundance, abundance trend and productivity) and extent of protected area
# (all protected areas, SSSIs, SPAs and SACs).
# The range and population measures were obtained from Atlas data, BBS data and 
# CES data (all available on request from BTO), and spatially matched to land 
# cover data (Land Cover Map (LCM) 2015, 1km proportion aggregate class from Great
# Britain and Northern Ireland, Rowland et al., 2017a, b), northings, easting, 
# elevation (REF), population density (REF) and proportion cover of each type of 
# protected area. 

# Step 1: log_population_density<-rescale(log(population_density+0.01)) # rescale population density

# Step 2: run species-specific models for species i
# This example uses cover of all protected areas (AllPA) but models are also run using SSSI, SPA and SAC cover)
# LCM proportion covers are labelled: Bro (broadleaved woodland), 
# Con (conifer woodland), IG (improved grassland), SNG (semi-natural grassland),
# MHB (mountainous), SW (saltwater), FW (freshwater), Coa (coastal), 
# Bui (built-up areas and gardens)
runmodels<-"no"
if(runmodels=="yes"){
#Atlas
model_allPA_occ<-bam(occupancy~allPA+Bro+Con+IG+SNG+MHB+SW+FW+Coa+Bui+log_population_density+I(log_population_density^2)+te(easting, northing, elevation, bs=c("tp","tp","tp")), family = binomial, data = Atlas_data_species_i)
model_allPA_col<-bam(colonisation~allPA+Bro+Con+IG+SNG+MHB+SW+FW+Coa+Bui+log_population_density+I(log_population_density^2)+te(easting, northing, elevation, bs=c("tp","tp","tp")), family = binomial, data = Atlas_data_species_i)
model_allPA_per<-bam(persistence~allPA+Bro+Con+IG+SNG+MHB+SW+FW+Coa+Bui+log_population_density+I(log_population_density^2)+te(easting, northing, elevation, bs=c("tp","tp","tp")), family = binomial, data = Atlas_data_species_i)

# Abundance and Trend
model_allPA_abundance<-bam(count~rescale(year)*allPA+I(rescale(year)^2)+Bro+Con+IG+SNG+MHB+SW+FW+Coa+Bui+log_population_density+I(log_population_density^2)+te(easting, northing, elevation, bs=c("tp","tp","tp"))+s(year2, bs="re"), offset=log(Nsq), family = nb() , discrete = TRUE, data = BBS_data_species_i)

#Productivity
model_allPA_productivity<-bam(cbind(Juv_count,Ad_count)~rescale(year)*allPA+Bro+Con+IG+SNG+MHB+SW+FW+Coa+Bui+log_population_density+I(log_population_density^2)+te(easting, northing, elevation, bs=c("tp","tp","tp"))+s(year2, bs="re"), family = "binomial" , discrete = TRUE, data = CES_data_species_i)

# The models were assessed using various methods described in the paper. 
# For all models that met the inclusion criteria, the coefficients for cover 
# of protected area (or the interaction between year and PA cover for abundance 
# trends) were collated into the tables in the data folder:
# https://doi.org/10.6084/m9.figshare.20200895

#########################################################
### Section 1b Methods for matching
### The analysis above was repeated using matched samples for Atlas and BBS data

# For a given PA type:
# Using the dataset described above, identify whether squares are treatment 
# squares (PA cover>10%) or control squares (PA cover = 0%). 

data1<-subset(data, PAtype>0.1 | PAtype==0)
pa_squares<-distinct(data1, gridref1, .keep_all = TRUE) #
pa_squares<-na.omit(pa_squares) # matching doesn't like NAs in any of the variables

# Carry out the matching:
m.out.nn <- matchit(treatment_or_control ~ Bro+Con+IG+SNG+MHB+SW+FW+Coa+Bui+easting+northing+elevation+log_popdens,
                    method = "nearest", distance = "mahalanobis")
## Extracting matched dataset with distance (not for Mahalanobis), weights and subclass attached
m.data.nn <- match.data(m.out.nn)
#merging the matched data with the bird data to do the modeling
pa_matched_nn<-merge(data2, m.data.nn[,c("gridref1","weights","subclass")], by="gridref1")


# Final models:
# Atlas
model_allPA_occ<-bam(occupancy~treatment_or_control+Bro+Con+IG+SNG+MHB+SW+FW+Coa+Bui+log_population_density+I(log_population_density^2)+te(easting, northing, elevation, bs=c("tp","tp","tp")), family = binomial,  weights=weights, data = pa_matched_nn_species_i)
model_allPA_col<-bam(colonisation~treatment_or_control+Bro+Con+IG+SNG+MHB+SW+FW+Coa+Bui+log_population_density+I(log_population_density^2)+te(easting, northing, elevation, bs=c("tp","tp","tp")), family = binomial,  weights=weights, data = pa_matched_nn_species_i)
model_allPA_per<-bam(persistence~treatment_or_control+Bro+Con+IG+SNG+MHB+SW+FW+Coa+Bui+log_population_density+I(log_population_density^2)+te(easting, northing, elevation, bs=c("tp","tp","tp")), family = binomial,  weights=weights, data = pa_matched_nn_species_i)

# Abundance and Trend
model_allPA_abundance<-bam(count~rescale(year)*treatment_or_control+I(rescale(year)^2)+Bro+Con+IG+SNG+MHB+SW+FW+Coa+Bui+log_population_density+I(log_population_density^2)+te(easting, northing, elevation, bs=c("tp","tp","tp"))+s(year2, bs="re"),offset=log(Nsq), family = nb(),  weights=weights , discrete = TRUE, data = pa_matched_nn_species_i)

#Productivity
model_allPA_productivity<-bam(cbind(Juv_count,Ad_count)~rescale(year)*treatment_or_control+Bro+Con+IG+SNG+MHB+SW+FW+Coa+Bui+log_population_density+I(log_population_density^2)+te(easting, northing, elevation, bs=c("tp","tp","tp"))+s(year2, bs="re"), family = "binomial" , discrete = TRUE, data = CES_data_species_i)

}
#######################################################################################################

##########   Section 2
# Summaries for model coefficients

# Section 2a: Table 1 & Figure 1b:
# numbers of species with negative and positive responses

###########

metric<-1
# Create tables for numbers and percent of species with significant and non-sig
# positive and negative responses to protected areas:
numberstab<-matrix(NA,5,16)
percenttab<-matrix(NA,5,16)
# Choose which population measure to summarise:
for(metric in 1:5){
  metriclist<-c("Occupancy", "Colonisation", "Persistence", "Abundance", "Trend")[metric]
  
# read in the data
  metric1<-metriclist
  if(metriclist=="Abundance" | metriclist=="Trend") metric1<-"BBS"
tab<-data.frame(read_xlsx("Suppl File 1.xlsx",sheet=metric1))
head(tab)

PAtype<-1
for(PAtype in 1:4){
  if(metric<4) PAtypelist<-c("AllPA","SSSI", "SPA", "SAC")[PAtype]
  if(metric==4) PAtypelist<-c("pa_abund","sssi_abund", "spa_abund", "sac_abund")[PAtype]
  if(metric==5) PAtypelist<-c("pa_trend","sssi_trend", "spa_trend", "sac_trend")[PAtype]
  
  if(metric<4){p_value_label<-paste0(PAtypelist,"_P")
  }else{p_value_label<-paste0(PAtypelist,"_p")  }

  numberstab[metric,((PAtype-1)*4+1):(PAtype*4)]<-as.vector(table(tab[,paste0(PAtypelist,"_coeff")]>0,as.numeric(tab[,p_value_label])<0.05))[c(4,2,3,1)]
  percenttab[metric,((PAtype-1)*4+1):(PAtype*4)]<-as.vector(prop.table(table(tab[,paste0(PAtypelist,"_coeff")]>0,as.numeric(tab[,p_value_label])<0.05)))[c(4,2,3,1)]*100
}
}

colnames(numberstab)<-c("allPApossSig", "allPAposnon-sig", "allPAnegsig", "allPAnegnon-sig", "SSSIpossig", "SSSIposnon-sig", "SSSInegsig", "SSSInegnon-sig", "SPApossig", "SPAposnon-sig", "SPAnegsig", "SPAnegnon-sig", "SACpossig", "SACposnon-sig", "SACnegsig", "SACnegnon-sig")
rownames(numberstab)<-c("Occupancy", "Colonisation", "Persistence", "Abundance", "Trend")
colnames(percenttab)<-c("allPApossSig", "allPAposnon-sig", "allPAnegsig", "allPAnegnon-sig", "SSSIpossig", "SSSIposnon-sig", "SSSInegsig", "SSSInegnon-sig", "SPApossig", "SPAposnon-sig", "SPAnegsig", "SPAnegnon-sig", "SACpossig", "SACposnon-sig", "SACnegsig", "SACnegnon-sig")
rownames(percenttab)<-c("Occupancy", "Colonisation", "Persistence", "Abundance", "Trend")

round(percenttab*100)/100

# Chi-squared tests for numbers of significant positive vs significant negative:
proptests_sig<-matrix(NA,5,4) 
chival_sig<-matrix(NA,5,4) 

measure<-4 ; count<-3
for(measure in 1:5){
  for(count in 1:4){
    try(proptests_sig[measure,count]<-(prop.test(numberstab[measure,1+(count-1)*4],numberstab[measure,1+(count-1)*4]+numberstab[measure,3+(count-1)*4])$p.value))
    try(chival_sig[measure,count]<-(prop.test(numberstab[measure,1+(count-1)*4],numberstab[measure,1+(count-1)*4]+numberstab[measure,3+(count-1)*4])$statistic))
  }}
round(proptests_sig*1000)/1000
round(chival_sig*100)/100

# Prop tests for productivity
prop.test(x=2, n=12)
prop.test(x=3, n=10)
prop.test(x=3, n=8)
prop.test(x=4, n=11)

#####

### Barplot
percenttab1<-percenttab
percenttab1[,c(3,4,7,8,11,12,15,16)]<-percenttab1[,c(3,4,7,8,11,12,15,16)]*-1


lightcols<-c("#6495ED","#FFA248", "#8FF0BC", "#9365A9",NA)
darkcols<-c("#2859B1","#E4660C","#53B480","#57296D",NA)

barplot(as.numeric(c(percenttab1[1,c(2,6,10,14)],NA,percenttab1[2,c(2,6,10,14)],NA,percenttab1[3,c(2,6,10,14)],NA,percenttab1[4,c(2,6,10,14)],NA,percenttab1[5,c(2,6,10,14)])),
        beside=T,ylim=c(-78,78),xlim=c(0,23),space=0,col=lightcols ,axisnames = F,ylab="Percent of species")

barplot(as.numeric(c(percenttab1[1,c(1,5,9,13)],NA,percenttab1[2,c(1,5,9,13)],NA,percenttab1[3,c(1,5,9,13)],NA,percenttab1[4,c(1,5,9,13)],NA,percenttab1[5,c(1,5,9,13)])),
        offset=as.numeric(c(percenttab1[1,c(2,6,10,14)],NA,percenttab1[2,c(2,6,10,14)],NA,percenttab1[3,c(2,6,10,14)],NA,percenttab1[4,c(2,6,10,14)],NA,percenttab1[5,c(2,6,10,14)])),
        add=T,space=0, axes=F, axisnames=F, col=darkcols)

barplot(as.numeric(c(percenttab1[1,c(2,6,10,14)+2],NA,percenttab1[2,c(2,6,10,14)+2],NA,percenttab1[3,c(2,6,10,14)+2],NA,percenttab1[4,c(2,6,10,14)+2],NA,percenttab1[5,c(2,6,10,14)+2])),
        add=T,ylim=c(-75,75),xlim=c(0,22),space=0,col=lightcols ,axisnames = F,ylab="Species with +ve and -ve association with protected area (%)")

barplot(as.numeric(c(percenttab1[1,c(1,5,9,13)+2],NA,percenttab1[2,c(1,5,9,13)+2],NA,percenttab1[3,c(1,5,9,13)+2],NA,percenttab1[4,c(1,5,9,13)+2],NA,percenttab1[5,c(1,5,9,13)+2])),
        offset=as.numeric(c(percenttab1[1,c(2,6,10,14)+2],NA,percenttab1[2,c(2,6,10,14)+2],NA,percenttab1[3,c(2,6,10,14)+2],NA,percenttab1[4,c(2,6,10,14)+2],NA,percenttab1[5,c(2,6,10,14)+2])),
        add=T,space=0, axes=F, axisnames=F, col=darkcols)


mtext(c("Occurrence","Colonisations","Persistence","Abundance","Trend"),line=3,side=1,at=c(2,7,12,17,22),cex=1.2)
mtext(rep(c("All PAs","SSSI","SPA","SAC",""),5),line=0,side=1,at=seq(0.5,24,1),cex=0.75,las=1)

###########################################################################################################

# Section 2b: Table 2 & Figure 1a:
# mean and 95% CIS of species' response to PAs

###########

mean_tab<-matrix(NA,7,4)
se_tab<-matrix(NA,7,4)
p_tab<-matrix(NA,7,4)
# Choose which population measure to summarise:
metric<-1
for(metric in 1:7){
  metriclist<-c("Occupancy", "Colonisation", "Persistence", "Abundance", "Trend","Productivity","Productivity_trend")[metric]
  
  # read in the data
  metric1<-metriclist
  if(metriclist=="Abundance" | metriclist=="Trend") metric1<-"BBS"
  if(metriclist=="Productivity" | metriclist=="Productivity_trend") metric1<-"CES"
  tab<-data.frame(read_xlsx("Suppl File 1.xlsx",sheet=metric1))
  
  for(i in 1:dim(tab)[2]){tab[,i]<-as.numeric(tab[,i])}
  
  PAtype<-1
  for(PAtype in 1:4){
    if(metric<4) PAtypelist<-c("AllPA_coeff","SSSI_coeff", "SPA_coeff", "SAC_coeff")[PAtype]
    if(metric==4) PAtypelist<-c("pa_abund_coeff","sssi_abund_coeff", "spa_abund_coeff", "sac_abund_coeff")[PAtype]
    if(metric==5) PAtypelist<-c("pa_trend_coeff","sssi_trend_coeff", "spa_trend_coeff", "sac_trend_coeff")[PAtype]
    if(metric==6) PAtypelist<-c("allPA","SSSI", "SPA", "SAC")[PAtype]
    if(metric==7) PAtypelist<-c("allPA_int","SSSI_int", "SPA_int", "SAC_int")[PAtype]
    
    # remove extremes:
    tab[tab[,PAtypelist]>10 & !is.na(tab[,PAtypelist]),PAtypelist]<-NA
    tab[tab[,PAtypelist]< -10 & !is.na(tab[,PAtypelist]),PAtypelist]<-NA
    
    mean_tab[metric,PAtype]<-summary(glm(tab[,PAtypelist]~1))$coefficients[1]
    se_tab[metric,PAtype]<-summary(glm(tab[,PAtypelist]~1))$coefficients[2]
    p_tab[metric,PAtype]<-summary(glm(tab[,PAtypelist]~1))$coefficients[4]
  }
  
  ## compare SPAs and SACs:
  if(metric<4) print(t.test(tab[,"SPA_coeff"],tab[,"SAC_coeff"],paired=T))
  if(metric==4) print(t.test(tab[,"spa_abund_coeff"],tab[,"sac_abund_coeff"],paired=T))
  if(metric==5) print(t.test(tab[,"spa_trend_coeff"],tab[,"sac_trend_coeff"],paired=T))
  if(metric==6) print(t.test(tab[,"SPA"],tab[,"SAC"],paired=T))
  if(metric==7) print(t.test(tab[,"SPA_int"],tab[,"SAC_int"],paired=T))
}
    
colnames(mean_tab)<-c("allPA","SSSI","SPA","SAC")
rownames(mean_tab)<-c("Occurrence","Colonisations","Persistence", "Abundance", "Trend","Productivity","Productivity_trend")
mean_tab
    
plot(mean_tab[1:5,1]~c(1:5),type="n",ylab="Effect",xlab="",xaxt="n",xlim=c(0.9,5.5),ylim=c(-0.2,0.7))
mtext(c("Occupancy", "Colonisation", "Persistence", "Abundance", "Trend"),side=1,at=c(1:5)+0.25)
for(i in 1:4){
  points(mean_tab[1:5,i]~c(1:5+(i-1)/6),pch=i)
segments(y0=mean_tab[1:5,i]+se_tab[1:5,i]*1.96,y1=mean_tab[1:5,i]-se_tab[1:5,i]*1.96,x0=c(1:5+(i-1)/6),x1=c(1:5+(i-1)/6))
}
lines(c(0,6),c(0,0),lty=2)

##################################################################################

### Section 3 Traits analysis

# choose whether to do simple version (popsizepredict<-"no", Table S5)
# or version with population size accounted for in predictions
# (popsizepredict=="yes", Table S6)
popsizepredict<-"yes"

# CHANGE THIS WHEN TRAITS ADDED TO FIGSHARE
traits_all<-read.csv("Analysis/Final_traits.csv")
row_name<-c("BoCC_Green","BoCC_Amber","BoCC_Red","a1no","a1yes","s1no","s1yes")

# Choose which population measure to summarise:
metric<-1
for(metric in 1:5){
  metriclist<-c("Occupancy", "Colonisation", "Persistence", "Abundance", "Trend","Productivity","Productivity_trend")[metric]
  
  # read in the data
  metric1<-metriclist
  if(metriclist=="Abundance" | metriclist=="Trend") metric1<-"BBS"
  if(metriclist=="Productivity" | metriclist=="Productivity_trend") metric1<-"CES"
  tab<-data.frame(read_xlsx("Suppl File 1.xlsx",sheet=metric1))
  for(count in 2:dim(tab)[2]){tab[,count]<-as.numeric(tab[,count])}
  
  traits<-traits_all[match(tab$Species,traits_all$speccode),]
  # This line gets rid of the undesigniated species:
  traits$BoCC1[traits$BoCC1=="-"]<-NA ; traits$BoCC1 <- factor(traits$BoCC1, levels = c("G","A","R"))
  
  # this sets up the results tables for number to be put in later on
  means_w<-matrix(NA,7,4) ; se_w<-matrix(NA,7,4)
  p_w<-matrix(NA,7,4) ;   p_w_i<-matrix(NA,7,4)
  lci_w<-matrix(NA,7,4) ; uci_w<-matrix(NA,7,4)
  
  
  # get just the PA type required:
  PAtype<-1
  for(PAtype in 1:4){
    if(metric<4) PAtypelist<-c("AllPA","SSSI", "SPA", "SAC")[PAtype]
    if(metric==4) PAtypelist<-c("pa_abund","sssi_abund", "spa_abund", "sac_abund")[PAtype]
    if(metric==5) PAtypelist<-c("pa_trend","sssi_trend", "spa_trend", "sac_trend")[PAtype]

    # select the correct PAtype data:
    tab2<-cbind(tab[,1],tab[,grepl(PAtypelist, colnames(tab))])
    
    tab3<-tab2[!is.na(tab2[,2]),] # this deletes the NAs
    logpsize<-traits$logpsize[!is.na(tab2[,2])]
    
    i<-1
    for(i in 1:3){ # This selects the conservation status type
      if(i==1){status<-traits$BoCC1[!is.na(tab2[,2])]}
      if(i==2){status<-traits$Annex1[!is.na(tab2[,2])]}
      if(i==3){status<-traits$Schedule1[!is.na(tab2[,2])]}
      
      # Models of response to PA ~ conservation status:
      if(popsizepredict=="no"){
      try(mod<-lm(tab3[,2]~status,weights=(1/(tab3[,3]^2))))
      # This gets estimates and SEs for species groups without popsize taken into account
      tagd_i<-predict(mod,newdata=data.frame(status=levels(factor(status))),se.fit=T)
      }
      if(popsizepredict=="yes"){
      mod<-lm(tab3[,2]~status+rescale(logpsize),weights=(1/(tab3[,3]^2)))
      # This gets estimates and SEs for species groups with popsize taken into account
      tagd_i<-predict(mod,newdata=data.frame(status=levels(factor(status)),logpsize=tapply(logpsize,status,mean,na.rm=T)),se.fit=T)
      }
      # Below puts the results into the results tables
      if(i==1){ means_w[1:3,PAtype]<-tagd_i$fit ; se_w[1:3,PAtype]<-tagd_i$se.fit
      lci_w[1:3,PAtype]<-c(tagd_i$fit-tagd_i$se.fit*1.96)
      uci_w[1:3,PAtype]<-c(tagd_i$fit+tagd_i$se.fit*1.96)
      p_w_i[2:3,PAtype]<-summary(mod)$coefficients[2:3,4]  }
      
      if(i==2){ means_w[4:5,PAtype]<-tagd_i$fit ; se_w[4:5,PAtype]<-tagd_i$se.fit
      lci_w[4:5,PAtype]<-c(tagd_i$fit-tagd_i$se.fit*1.96)
      uci_w[4:5,PAtype]<-c(tagd_i$fit+tagd_i$se.fit*1.96)
      p_w_i[5,PAtype]<-summary(mod)$coefficients[2,4] }
      if(i==3){ means_w[6:7,PAtype]<-tagd_i$fit ; se_w[6:7,PAtype]<-tagd_i$se.fit
      lci_w[6:7,PAtype]<-c(tagd_i$fit-tagd_i$se.fit*1.96)
      uci_w[6:7,PAtype]<-c(tagd_i$fit+tagd_i$se.fit*1.96)
      p_w_i[7,PAtype]<-summary(mod)$coefficients[2,4]}
    }
    
    # this gets results tables labelled and in the right shape:
    tab_A1<-matrix(paste0(round(means_w*1000)/1000," (",round(lci_w*100)/100,", ",round(uci_w*100)/100,"), ",round(p_w_i*1000)/1000),7,4) 
    tab_A1[c(1,4,6),]<-matrix(paste0(round(means_w[c(1,4,6),]*1000)/1000," (",round(lci_w[c(1,4,6),]*100)/100,", ",round(uci_w[c(1,4,6),]*100)/100,")"),3,4) 
    colnames(tab_A1)<-paste0(metriclist,c("allPA","SSSI","SPA","SAC")) ;  rownames(tab_A1)<-row_name
    tab_m1<-matrix(means_w,7,4) ; colnames(tab_m1)<-paste0(metriclist,c("allPA","SSSI","SPA","SAC")) ; rownames(tab_m1)<-row_name
    tab_se1<-matrix(se_w,7,4) ; colnames(tab_se1)<-paste0(metriclist,c("allPA","SSSI","SPA","SAC")) ; rownames(tab_se1)<-row_name
    tab_lci1<-matrix(lci_w,7,4) ; colnames(tab_lci1)<-paste0(metriclist,c("allPA","SSSI","SPA","SAC")) ; rownames(tab_lci1)<-row_name
    tab_uci1<-matrix(uci_w,7,4) ; colnames(tab_uci1)<-paste0(metriclist,c("allPA","SSSI","SPA","SAC")) ; rownames(tab_uci1)<-row_name
    tab_p_i1<-matrix(p_w_i,7,4) ; colnames(tab_p_i1)<-paste0(metriclist,c("allPA","SSSI","SPA","SAC")) ; rownames(tab_p_i1)<-row_name
    
  }
  
  # This combines results tables as you go through the loops:
  if(metric==1){tab_A<-tab_A1}else{tab_A<-cbind(tab_A,tab_A1)}
  if(metric==1){tab_m<-tab_m1}else{tab_m<-cbind(tab_m,tab_m1)}
  if(metric==1){tab_se<-tab_se1}else{tab_se<-cbind(tab_se,tab_se1)}
  if(metric==1){tab_lci<-tab_lci1}else{tab_lci<-cbind(tab_lci,tab_lci1)}
  if(metric==1){tab_uci<-tab_uci1}else{tab_uci<-cbind(tab_uci,tab_uci1)}
  if(metric==1){tab_p_i<-tab_p_i1}else{tab_p_i<-cbind(tab_p_i,tab_p_i1)}
  print(tab_A1)
}

rbind(tab_A[,c(1:4)],tab_A[,c(5:8)],tab_A[,c(9:12)],tab_A[,c(13:16)],tab_A[,c(17:20)])

rm(logpsize)
####################################################################################################

# section 3c: ecological traits analysis

# First create a phylogenetic tree from Jetz data
#tree2<-read.tree("Analysis/EricsonStage1_myTips_hooded.tre", keep.multi=TRUE) 
#avtree <- averageTree(tree2,method="symmetric.difference")
#write.tree(avtree,"Analysis/myaveragedtree.tre")
avtree<-read.tree("Analysis/myaveragedtree.tre", keep.multi=TRUE)

traits_all<-read.csv("Analysis/Final_traits.csv") # This is the traits file copied from google drive

##################################################

# this sets up the results tables for number to be put in later on
habitat_table<-array(NA, dim=c(7,4,5,4))
dimnames(habitat_table)<-list(c("wetland","upland","coastal","farmland","woodland","urban","unclassified"),
                              c("mean","lCI","uCI","P"),c("Occupancy","Colonisation","Persistence","Abundance","Trend"),(c("allPA","SSSI","SPA","SAC")))

vars_table<-array(NA, dim=c(5,4,5,4))
dimnames(vars_table)<-list(c("Log_Mass","logpchange","SSI","logpsize","STI"),
                           c("mean","lCI","uCI","P"),c("Occupancy","Colonisation","Persistence","Abundance","Trend"),(c("allPA","SSSI","SPA","SAC")))


# Choose which population measure to summarise:
metric<-1
for(metric in 1:5){
  metriclist<-c("Occupancy", "Colonisation", "Persistence", "Abundance", "Trend")[metric]
  
  # read in the data
  metric1<-metriclist
  if(metriclist=="Abundance" | metriclist=="Trend") metric1<-"BBS"
  tab<-data.frame(read_xlsx("Suppl File 1.xlsx",sheet=metric1))
  for(count in 2:dim(tab)[2]){tab[,count]<-as.numeric(tab[,count])}
  
  # make sure your results and traits tables are in the same order of species
  traits<-traits_all[match(tab$Species,traits_all$speccode),]

  # get just the PA type required:
  PAtype<-1
  for(PAtype in 1:4){
    if(metric<4) PAtypelist<-c("AllPA","SSSI", "SPA", "SAC")[PAtype]
    if(metric==4) PAtypelist<-c("pa_abund","sssi_abund", "spa_abund", "sac_abund")[PAtype]
    if(metric==5) PAtypelist<-c("pa_trend","sssi_trend", "spa_trend", "sac_trend")[PAtype]
    
    # select the correct PAtype data:
    tab2<-cbind(tab[,1],tab[,grepl(PAtypelist, colnames(tab))])
    traitstag<-cbind(tab2,traits)
    colnames(traitstag)[2:4]<-c("coeff","SE","P")

    
    # Get rid of the species that are not included in analysis:
    traitstag<-traitstag[!is.na(traitstag[,2]) & !is.na(rowSums(traitstag[,names(traitstag) %in% c("Log_Mass","logpchange","logpsize","SSI","STI")])),]

    traitstag$animal<-traitstag$tip.label # this is needed for the model
    rm(tab2)
    
    
    ### Remake the tree with only required species in:
    results <- lapply(avtree, keep.tip,  traitstag$tip.label)
    tree3 <- do.call(c, results)
    
    ### Reorder traits table to match the species order of the tree
    traitstag<-traitstag[match(tree3[[1]]$tip.label,traitstag$tip.label),]

    ## set up phylogeny for inclusion
    sptree<-makeNodeLabel(tree3[[1]], method="number", prefix="node")

    ## prior with expanded parameters for MCMCglmm 
    prior <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),G2 = list(V = 1, fix=1)), R = list(V = 1, nu = 0.002))

    model1 <- MCMCglmm(coeff~rescale(Log_Mass)+rescale(logpchange)+rescale(SSI)+rescale(logpsize)+rescale(STI)+habitat,random=~animal+idh(SE):units, family = "gaussian", prior = prior, data = traitstag, pedigree=sptree,nitt = 50000, burnin = 5000, thin = 25,verbose=F)
    
    #print(c("metric",metric,"PAtype",PAtype))
    #print(summary(model1)$solutions)      
    #plot(model1$Sol)
    #plot(model1$VCV)
    
    # habitats and other variables:
    vars_table[,1:4,metric,PAtype]<-summary(model1)$solutions[2:6,c(1:3,5)]   
    habitat_table[,1,metric,PAtype]<-c(mean(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatwetland']),
                                               mean(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatupland']),
                                               mean(model1$Sol[,'(Intercept)']),
                                               mean(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatfarmland']),
                                               mean(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatwoodland']),
                                               mean(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitaturban']),
                                               mean(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatunclassified']))
    
    habitat_table[,2:3,metric,PAtype]<-rbind(quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatwetland'], prob = c(0.025,0.975)),
                                                     quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatupland'], prob = c(0.025,0.975)),
                                                     quantile(model1$Sol[,'(Intercept)'], prob = c(0.025,0.975)),
                                                     quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatfarmland'], prob = c(0.025,0.975)),
                                                     quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatwoodland'], prob = c(0.025,0.975)),
                                                     quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitaturban'], prob = c(0.025,0.975)),
                                                     quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatunclassified'], prob = c(0.025,0.975)))
    tag<-matrix(NA,7,2)
    tag<-rbind(quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatwetland'], prob = c(0.0025,0.9975)),
               quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatupland'], prob = c(0.0025,0.9975)),
               quantile(model1$Sol[,'(Intercept)'], prob = c(0.0025,0.9975)),
               quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatfarmland'], prob = c(0.0025,0.9975)),
               quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatwoodland'], prob = c(0.0025,0.9975)),
               quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitaturban'], prob = c(0.0025,0.9975)),
               quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatunclassified'], prob = c(0.0025,0.9975)))
    tag2<-matrix(NA,7,2)
    tag2<-rbind(quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatwetland'], prob = c(0.0005,0.9995)),
                quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatupland'], prob = c(0.0005,0.9995)),
                quantile(model1$Sol[,'(Intercept)'], prob = c(0.0005,0.9995)),
                quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatfarmland'], prob = c(0.0005,0.9995)),
                quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatwoodland'], prob = c(0.0005,0.9995)),
                quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitaturban'], prob = c(0.0005,0.9995)),
                quantile(model1$Sol[,'(Intercept)'] + model1$Sol[,'habitatunclassified'], prob = c(0.0005,0.9995)))
    
    
    habitat_table[as.numeric(habitat_table[,2,metric,PAtype])/sqrt(as.numeric(habitat_table[,2,metric,PAtype])^2)+as.numeric(habitat_table[,3,metric,PAtype])/sqrt(as.numeric(habitat_table[,3,metric,PAtype])^2)==2,4,metric,PAtype]<-"*"
    habitat_table[tag[,1]/sqrt(tag[,1]^2)+tag[,2]/sqrt(tag[,2]^2)==2,4,metric,PAtype]<-"**"
    habitat_table[tag2[,1]/sqrt(tag2[,1]^2)+tag2[,2]/sqrt(tag2[,2]^2)==2,4,metric,PAtype]<-"***"
    habitat_table[,,metric,PAtype]
    
  }
  rm(tab)
  }

##############

dimnames(vars_table)
par(mfrow=c(1,1))
tab<-as.data.frame(cbind(habitat_table[,1:3,"Occupancy","allPA"],habitat_table[,1:3,"Colonisation","allPA"],habitat_table[,1:3,"Persistence","allPA"],habitat_table[,1:3,"Abundance","allPA"],habitat_table[,1:3,"Trend","allPA"]))
#colnames(tab)<-c("mean_occ_PA", "lCI_occ_PA", "uCI_occ_PA", "mean_col_PA", "lCI_col_PA", "uCI_col_PA", "mean_per_PA", "lCI_per_PA", "uCI_per_PA")

coeffs<-rbind(vars_table[c(1,4,2,3,5),1,,1],habitat_table[,1,,1])
LCI<-rbind(vars_table[c(1,4,2,3,5),2,,1],habitat_table[,2,,1])
UCI<-rbind(vars_table[c(1,4,2,3,5),3,,1],habitat_table[,3,,1])

spacing<-0.13

names<-c("Mass","Pop\nsize","Pop\nchange","SSI","STI","Wetland","Upland","Coastal","Farmland","Woodland","Urban","Unclassified")
plot(as.numeric(coeffs[,1])~c(1:12),pch=16,xlim=c(0.9,12.6),ylim=c(min(as.numeric(LCI)),max(as.numeric(UCI))),
     xlab="",ylab="Effect",xaxt="n",cex=1.1,las=2)
mtext(names,at=c(1:12)+0.2,side=1,cex=1,line=1)
mtext("Trait",side=1,cex=1.1,line=2.5)

for(i in 1:5){segments(1:12+spacing*i-spacing,as.numeric(LCI[,i]),1:12+spacing*i-spacing,as.numeric(UCI[,i]))}
lines(c(0,13),c(0,0),lty=2)

points(coeffs[,2]~c(c(1:12)+spacing),pch=17,cex=1.2)
points(coeffs[,3]~c(c(1:12)+spacing*2),pch=18,cex=1.2)
points(coeffs[,4]~c(c(1:12)+spacing*3),cex=1.2,pch=21,col=1, bg="white")
points(coeffs[,5]~c(c(1:12)+spacing*4),cex=1.2,pch=24,col=1, bg="white")


################################################################################################################
#Figure 2
#Plot PA extent coefficients from productivity models against those from abundance/occupancy/colonisation/persistence models
PA_coeffs_df_1km<-read.csv("PA_coeffs_df_1km.csv")

parmests_BBS<-read.csv("parmests with popdens - BBS.csv")
parmests_occ<-read.csv("parmests with popdens - Occupancy.csv")
parmests_col<-read.csv("parmests with popdens - Colonisation.csv")
parmests_per<-read.csv("parmests with popdens - Persistence.csv")

PA_coeffs_df_1km$BBS_code<-character(nrow(PA_coeffs_df_1km))
PA_coeffs_df_1km$BBS_SPA<-PA_coeffs_df_1km$BBS_SPA_se<-
  PA_coeffs_df_1km$BBS_allPA<-PA_coeffs_df_1km$BBS_allPA_se<-
  PA_coeffs_df_1km$BBS_SSSI<-PA_coeffs_df_1km$BBS_SSSI_se<-
  PA_coeffs_df_1km$BBS_SAC<-PA_coeffs_df_1km$BBS_SAC_se<-
  PA_coeffs_df_1km$BBSyr_SPA<-PA_coeffs_df_1km$BBSyr_SPA_se<-
  PA_coeffs_df_1km$BBSyr_allPA<-PA_coeffs_df_1km$BBSyr_allPA_se<-
  PA_coeffs_df_1km$BBSyr_SSSI<-PA_coeffs_df_1km$BBSyr_SSSI_se<-
  PA_coeffs_df_1km$BBSyr_SAC<-PA_coeffs_df_1km$BBSyr_SAC_se<-
  PA_coeffs_df_1km$occ_SPA<-PA_coeffs_df_1km$occ_SPA_se<-
  PA_coeffs_df_1km$occ_allPA<-PA_coeffs_df_1km$occ_allPA_se<-
  PA_coeffs_df_1km$occ_SSSI<-PA_coeffs_df_1km$occ_SSSI_se<-
  PA_coeffs_df_1km$occ_SAC<-PA_coeffs_df_1km$occ_SAC_se<-
  PA_coeffs_df_1km$col_SPA<-PA_coeffs_df_1km$col_SPA_se<-
  PA_coeffs_df_1km$col_allPA<-PA_coeffs_df_1km$col_allPA_se<-
  PA_coeffs_df_1km$col_SSSI<-PA_coeffs_df_1km$col_SSSI_se<-
  PA_coeffs_df_1km$col_SAC<-PA_coeffs_df_1km$col_SAC_se<-
  PA_coeffs_df_1km$per_SPA<-PA_coeffs_df_1km$per_SPA_se<-
  PA_coeffs_df_1km$per_allPA<-PA_coeffs_df_1km$per_allPA_se<-
  PA_coeffs_df_1km$per_SSSI<-PA_coeffs_df_1km$per_SSSI_se<-
  PA_coeffs_df_1km$per_SAC<-PA_coeffs_df_1km$per_SAC_se<-numeric(nrow(PA_coeffs_df_1km))

for(i in 1:nrow(PA_coeffs_df_1km)){
  PA_coeffs_df_1km$BBS_code[i]<-traits[traits$english_name==PA_coeffs_df_1km$Species[i],]$speccode
  PA_coeffs_df_1km$BBS_SPA[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$spa_abund_coeff
  PA_coeffs_df_1km$BBS_SPA_se[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$spa_abund_se
  PA_coeffs_df_1km$BBS_allPA[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$pa_abund_coeff
  PA_coeffs_df_1km$BBS_allPA_se[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$pa_abund_se
  PA_coeffs_df_1km$BBS_SSSI[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$sssi_abund_coeff
  PA_coeffs_df_1km$BBS_SSSI_se[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$sssi_abund_se
  PA_coeffs_df_1km$BBS_SAC[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$sac_abund_coeff
  PA_coeffs_df_1km$BBS_SAC_se[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$sac_abund_se
  PA_coeffs_df_1km$BBSyr_SPA[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$spa_trend_coeff
  PA_coeffs_df_1km$BBSyr_SPA_se[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$spa_trend_se
  PA_coeffs_df_1km$BBSyr_allPA[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$pa_trend_coeff
  PA_coeffs_df_1km$BBSyr_allPA_se[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$pa_trend_se
  PA_coeffs_df_1km$BBSyr_SSSI[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$sssi_trend_coeff
  PA_coeffs_df_1km$BBSyr_SSSI_se[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$sssi_trend_se
  PA_coeffs_df_1km$BBSyr_SAC[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$sac_trend_coeff
  PA_coeffs_df_1km$BBSyr_SAC_se[i]<-parmests_BBS[parmests_BBS$Species==PA_coeffs_df_1km$BBS_code[i],]$sac_trend_se
  PA_coeffs_df_1km$occ_SPA[i]<-parmests_occ[parmests_occ$Species==PA_coeffs_df_1km$BBS_code[i],]$SPA_coeff
  PA_coeffs_df_1km$occ_SPA_se[i]<-parmests_occ[parmests_occ$Species==PA_coeffs_df_1km$BBS_code[i],]$SPA_se
  PA_coeffs_df_1km$occ_allPA[i]<-parmests_occ[parmests_occ$Species==PA_coeffs_df_1km$BBS_code[i],]$AllPA_coeff
  PA_coeffs_df_1km$occ_allPA_se[i]<-parmests_occ[parmests_occ$Species==PA_coeffs_df_1km$BBS_code[i],]$AllPA_se
  PA_coeffs_df_1km$occ_SSSI[i]<-parmests_occ[parmests_occ$Species==PA_coeffs_df_1km$BBS_code[i],]$SSSI_coeff
  PA_coeffs_df_1km$occ_SSSI_se[i]<-parmests_occ[parmests_occ$Species==PA_coeffs_df_1km$BBS_code[i],]$SSSI_se
  PA_coeffs_df_1km$occ_SAC[i]<-parmests_occ[parmests_occ$Species==PA_coeffs_df_1km$BBS_code[i],]$SAC_coeff
  PA_coeffs_df_1km$occ_SAC_se[i]<-parmests_occ[parmests_occ$Species==PA_coeffs_df_1km$BBS_code[i],]$SAC_se
  PA_coeffs_df_1km$col_SPA[i]<-parmests_col[parmests_col$Species==PA_coeffs_df_1km$BBS_code[i],]$SPA_coeff
  PA_coeffs_df_1km$col_SPA_se[i]<-parmests_col[parmests_col$Species==PA_coeffs_df_1km$BBS_code[i],]$SPA_se
  PA_coeffs_df_1km$col_allPA[i]<-parmests_col[parmests_col$Species==PA_coeffs_df_1km$BBS_code[i],]$AllPA_coeff
  PA_coeffs_df_1km$col_allPA_se[i]<-parmests_col[parmests_col$Species==PA_coeffs_df_1km$BBS_code[i],]$AllPA_se
  PA_coeffs_df_1km$col_SSSI[i]<-parmests_col[parmests_col$Species==PA_coeffs_df_1km$BBS_code[i],]$SSSI_coeff
  PA_coeffs_df_1km$col_SSSI_se[i]<-parmests_col[parmests_col$Species==PA_coeffs_df_1km$BBS_code[i],]$SSSI_se
  PA_coeffs_df_1km$col_SAC[i]<-parmests_col[parmests_col$Species==PA_coeffs_df_1km$BBS_code[i],]$SAC_coeff
  PA_coeffs_df_1km$col_SAC_se[i]<-parmests_col[parmests_col$Species==PA_coeffs_df_1km$BBS_code[i],]$SAC_se
  if(length(parmests_per[parmests_per$Species==PA_coeffs_df_1km$BBS_code[i],]$SPA_coeff)==0) {
    PA_coeffs_df_1km$per_SPA[i]<-NA
    PA_coeffs_df_1km$per_SPA_se[i]<-NA
    PA_coeffs_df_1km$per_allPA[i]<-NA
    PA_coeffs_df_1km$per_allPA_se[i]<-NA
    PA_coeffs_df_1km$per_SSSI[i]<-NA
    PA_coeffs_df_1km$per_SSSI_se[i]<-NA
    PA_coeffs_df_1km$per_SAC[i]<-NA
    PA_coeffs_df_1km$per_SAC_se[i]<-NA
  } else{
    PA_coeffs_df_1km$per_SPA[i]<-parmests_per[parmests_per$Species==PA_coeffs_df_1km$BBS_code[i],]$SPA_coeff
    PA_coeffs_df_1km$per_SPA_se[i]<-parmests_per[parmests_per$Species==PA_coeffs_df_1km$BBS_code[i],]$SPA_se
    PA_coeffs_df_1km$per_allPA[i]<-parmests_per[parmests_per$Species==PA_coeffs_df_1km$BBS_code[i],]$AllPA_coeff
    PA_coeffs_df_1km$per_allPA_se[i]<-parmests_per[parmests_per$Species==PA_coeffs_df_1km$BBS_code[i],]$AllPA_se
    PA_coeffs_df_1km$per_SSSI[i]<-parmests_per[parmests_per$Species==PA_coeffs_df_1km$BBS_code[i],]$SSSI_coeff
    PA_coeffs_df_1km$per_SSSI_se[i]<-parmests_per[parmests_per$Species==PA_coeffs_df_1km$BBS_code[i],]$SSSI_se
    PA_coeffs_df_1km$per_SAC[i]<-parmests_per[parmests_per$Species==PA_coeffs_df_1km$BBS_code[i],]$SAC_coeff
    PA_coeffs_df_1km$per_SAC_se[i]<-parmests_per[parmests_per$Species==PA_coeffs_df_1km$BBS_code[i],]$SAC_se
  }
}

plot_df<-data.frame(matrix(nrow=5,ncol=2))
names(plot_df)<-c("response_var","ylab_text")
plot_df$response_var<-c("BBS","BBSyr","occ","col","per")
plot_df$ylab_text<-c("abundance","abundance trend","occupancy","colonisations","persistence")

#PA term first
for(i in 1:nrow(plot_df)){
  for(j in 1:length(PA_type)){
    par(mar=c(4.1, 4.1, 1.1, 1.1))
    temp_form<-as.formula(paste0(plot_df$response_var[i],"_",PA_type[j],"~",PA_type[j]))
    y_var_se<-PA_coeffs_df_1km[,paste0(plot_df$response_var[i],"_",PA_type[j],"_se")]
    y_var<-PA_coeffs_df_1km[,paste0(plot_df$response_var[i],"_",PA_type[j])]
    x_var_se<-PA_coeffs_df_1km[,paste0(PA_type[j],"_se")]
    x_var<-PA_coeffs_df_1km[,PA_type[j]]
    plot(temp_form,
         PA_coeffs_df_1km,
         bty="n",
         pch=21,
         xlim=c(min(x_var-x_var_se),max(x_var+x_var_se)),
         ylim=c(min(y_var-y_var_se,na.rm=T),max(y_var+y_var_se,na.rm=T)),
         bg="seagreen3",
         las=1,
         xlab=paste0(PA_type[j]," coefficient (productivity)"),
         ylab=paste0(PA_type[j]," coefficient (",plot_df$ylab_text[i],")"))
    arrows(x_var-x_var_se,y_var,x_var+x_var_se,y_var,length=0,col="gray50")
    arrows(x_var,y_var-y_var_se,x_var,y_var+y_var_se,length=0,col="gray50")
    m1<-lm(y_var~x_var,PA_coeffs_df_1km)
    newx<-seq(min(x_var-x_var_se),max(x_var+x_var_se),by=0.01)
    ps<-predict(m1, newdata=data.frame(x_var=newx), interval="confidence", level = 0.95)
    lines(newx,ps[,1],lty=1)
    lines(newx,ps[,2],lty=2)
    lines(newx,ps[,3],lty=2)
  }
}

#PA*time term next
for(i in 1:nrow(plot_df)){
  for(j in 1:length(PA_type)){
    par(mar=c(4.1, 4.1, 1.1, 1.1))
    temp_form<-as.formula(paste0(plot_df$response_var[i],"_",PA_type[j],"~",PA_type[j],"_int"))
    y_var_se<-PA_coeffs_df_1km[,paste0(plot_df$response_var[i],"_",PA_type[j],"_se")]
    y_var<-PA_coeffs_df_1km[,paste0(plot_df$response_var[i],"_",PA_type[j])]
    x_var_se<-PA_coeffs_df_1km[,paste0(PA_type[j],"_int_se")]
    x_var<-PA_coeffs_df_1km[,paste0(PA_type[j],"_int")]
    plot(temp_form,
         PA_coeffs_df_1km,
         bty="n",
         pch=21,
         xlim=c(min(x_var-x_var_se),max(x_var+x_var_se)),
         ylim=c(min(y_var-y_var_se,na.rm=T),max(y_var+y_var_se,na.rm=T)),
         bg="seagreen3",
         las=1,
         xlab=paste0(PA_type[j]," * time coefficient (productivity)"),
         ylab=paste0(PA_type[j]," coefficient (",plot_df$ylab_text[i],")"))
    arrows(x_var-x_var_se,y_var,x_var+x_var_se,y_var,length=0,col="gray50")
    arrows(x_var,y_var-y_var_se,x_var,y_var+y_var_se,length=0,col="gray50")
    m1<-lm(y_var~x_var,PA_coeffs_df_1km)
    newx<-seq(min(x_var-x_var_se),max(x_var+x_var_se),by=0.01)
    ps<-predict(m1, newdata=data.frame(x_var=newx), interval="confidence", level = 0.95)
    lines(newx,ps[,1],lty=1)
    lines(newx,ps[,2],lty=2)
    lines(newx,ps[,3],lty=2)
  }
}

