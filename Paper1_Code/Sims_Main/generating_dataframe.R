library(readr)
library(ggpubr)
setwd("/nesi/nobackup/uoa04039/holobiont_first_paper_sim_reps/sims_28Aug24/sims_complete_rerun_31oct//")
rds_files<-list.files(pattern="*RDS")
rds_files<-rds_files[grepl("analysed_replicate",rds_files)]
rds_files<-rds_files[grepl("2nov",rds_files)]

rds_files<-lapply(rds_files,read_rds)
for(i in 1:length(rds_files)){
  rds_files[[i]]$Replicate<-paste0("replicate_",i)
}

general_tm_1_10_data<-do.call("rbind",rds_files)


general_tm_1_10_data$MeanGroupingColumn<-paste(general_tm_1_10_data$Replicate,general_tm_1_10_data$EnvType,
                                               paste0("1gen",general_tm_1_10_data$Gen),sep=".")


general_tm_1_10_data$MeanGroupingColumn<-paste(general_tm_1_10_data$Gen,general_tm_1_10_data$VertCon,
                                               general_tm_1_10_data$EnvType,general_tm_1_10_data$NGens)
library(dplyr)
general_tm_1_10_data_gen1<-general_tm_1_10_data %>% filter(NGens == "1gen")
general_tm_1_10_data_gen1$NewGrouping_Column<-gsub(" 1gen","",general_tm_1_10_data_gen1$MeanGroupingColumn)
general_tm_1_10_data_gen1$NewGrouping_Column<-paste(general_tm_1_10_data_gen1$NewGrouping_Column,general_tm_1_10_data_gen1$Replicate)
general_tm_1_10_data$cleaned_strings <- gsub("\\S*gen\\S*", "", general_tm_1_10_data$MeanGroupingColumn)
general_tm_1_10_data$cleaned_strings <- trimws(general_tm_1_10_data$cleaned_strings)
general_tm_1_10_data$cleaned_strings<-paste(general_tm_1_10_data$cleaned_strings,general_tm_1_10_data$Replicate)

general_tm_1_10_data$gen1_microbefitness<-general_tm_1_10_data_gen1$HostMicrobeFitness[match(general_tm_1_10_data$cleaned_strings,general_tm_1_10_data_gen1$NewGrouping_Column)]
general_tm_1_10_data$gen1_hostfitness<-general_tm_1_10_data_gen1$HostFitness[match(general_tm_1_10_data$cleaned_strings,general_tm_1_10_data_gen1$NewGrouping_Column)]

general_tm_1_10_data$gen1_microbefitvar<-general_tm_1_10_data_gen1$hostmicrobe_fitvar[match(general_tm_1_10_data$cleaned_strings,general_tm_1_10_data_gen1$NewGrouping_Column)]
general_tm_1_10_data$gen1_hostfitvar<-general_tm_1_10_data_gen1$host_fitvar[match(general_tm_1_10_data$cleaned_strings,general_tm_1_10_data_gen1$NewGrouping_Column)]


general_tm_1_10_data<-general_tm_1_10_data %>% mutate(NGens=factor(NGens,levels=c("1gen","2gen","3gen","4gen","5gen","6gen","7gen","8gen","9gen",
                                                                                  "10gen","20gen","50gen","100gen","200gen")))

general_tm_1_10_data$LogFold_ngen1_Host<-(general_tm_1_10_data$HostFitness/general_tm_1_10_data$gen1_hostfitness)
general_tm_1_10_data$LogFold_ngen1_Microbe<-(general_tm_1_10_data$HostMicrobeFitness/general_tm_1_10_data$gen1_microbefitness)

general_tm_1_10_data$hostfitvar_ratio<-(general_tm_1_10_data$host_fitvar/general_tm_1_10_data$gen1_hostfitvar)
general_tm_1_10_data$microbefitvar_ratio<-(general_tm_1_10_data$hostmicrobe_fitvar/general_tm_1_10_data$gen1_microbefitvar)


general_tm_1_10_data_varfits<-general_tm_1_10_data %>%
  filter(Gen==1500)
general_tm_1_10_data_varfits$VertCon_Num<-as.numeric(general_tm_1_10_data_varfits$VertCon)
general_tm_1_10_data_varfits$NGens_Num<-as.numeric(gsub("[^0-9]","",general_tm_1_10_data_varfits$NGens))
general_tm_1_10_data_varfits %>%
  filter(NGens %in% c("200gen", "1gen")) %>%
  ggplot() +
  aes(x = VertCon, y = microbefitvar_ratio, fill = NGens,col=NGens) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_minimal() +
  facet_wrap(vars(EnvType),scales="free_y")


general_tm_1_10_data_varfits %>%
  # filter(NGens %in% "1gen") %>%
  ggplot() +
  aes(x = VertCon, y = hostmicrobe_fitvar,group=NGens) +
  geom_boxplot(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(EnvType), nrow=1)

N=20
random_levels <- sample(unique(general_tm_1_10_data$Replicate), N)
geomean<-function(x){
  exp(mean(log(x)))
}
#general_tm_1_10_data$diversity<-ifelse(is.na(general_tm_1_10_data$diversity),0,general_tm_1_10_data$diversity)
general_tm_1_10_data_means <- general_tm_1_10_data %>%
  mutate(diversity = ifelse(is.na(diversity) | is.nan(diversity), 0, diversity)) %>%
  filter(Replicate %in% random_levels) %>%
  # filter(EnvType =="lowac_incvar") %>%
  #filter(VertCon == 0 ) %>%
  group_by(EnvType,VertCon,NGens,Gen) %>%
  summarize(HostFitness.mean = mean(HostFitness),
            HostMicrobeFitness.mean = mean(HostMicrobeFitness),
            diversity.mean = mean(diversity),
            mean_reproductive_output.mean = mean(mean_reproductive_output),
            EnvMicrobeFitness.mean = mean(EnvMicrobeFitness),
            LogFoldChangeMicrobe.mean = mean(LogFoldChangeMicrobe),
            EnvCond.mean = mean(EnvCond),
            LogFoldChangeHost.mean = mean(LogFoldChangeHost),
            AdjRSq.mean = mean(adj_rsq),
            Coeff.mean = mean(host_microbe_coeff),
            LogFold_ngen1_Host.mean = mean(LogFold_ngen1_Host),
            LogFold_ngen1_Microbe.mean = mean(LogFold_ngen1_Microbe),
            BDiv.mean = mean(host_bdiv),
            FitVar_Host.mean=mean(host_fitvar),
            FitVar_Microbe.mean=mean(hostmicrobe_fitvar),
            mean_microbetrait=mean(mean_microbe_traitval),
            diversity.var=var(diversity),
            PhenoVar=var(mean_microbe_traitval)
  )



library(dplyr)

general_tm_1_10_data_means<-general_tm_1_10_data_means %>% mutate(NGens=factor(NGens,levels=c("1gen","2gen","3gen","4gen","5gen","6gen","7gen","8gen","9gen",
                                                                                              "10gen","20gen","50gen","100gen","200gen")))
general_tm_1_10_data_means$GenLab<-gsub("gen","",general_tm_1_10_data_means$NGens)
general_tm_1_10_data_means<-general_tm_1_10_data_means %>% 
  mutate(GenLab=factor(GenLab,levels=c("1","2","3","4","5","6","7","8","9","10","20","50","100","200")))


write.csv(general_tm_1_10_data_means,"general_tm_1_10_data_means.csv")
write_rds(general_tm_1_10_data_means,"general_tm_1_10_data_means.RDS")

