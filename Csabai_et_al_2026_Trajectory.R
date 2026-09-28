# FIRST STEPS

# Checking & installing missing packages
packages <- c("berryFunctions", "col2hex2col", "colorspace", "dplyr","ecotraj","emmeans", "forcats", "ggplot2", "ggpubr", "lmerTest","plyr", "rstatix","scales", "stringr", "strex", "vegan", "viridis")
install.packages(setdiff(packages, rownames(installed.packages())))
remove(packages)

# Load packages
library(berryFunctions)
library(col2hex2col)
library(colorspace)
library(dplyr)
library(ecotraj)
library(emmeans)
library(forcats)
library(ggplot2)
library(ggpubr)
library(lmerTest)
library(rstatix)
library(scales)
library(stringr)
library(strex)
library(vegan)
library(viridis)

# Start
setwd(getwd())

# Load data
data(data)
data(data_hydrology)
data(DD_365_FINAL)
DD_365_FINAL$DD_365 <- round(DD_365_FINAL$DD_365,0)


## -- BASIC STEPS & SETTINGS -- ##

# Create hydrology table
rownames(data_hydrology) <- data_hydrology$SITE
data_hydrology <- data_hydrology[,-1]
colnames(data_hydrology) <- str_sub(colnames(data_hydrology),2)
HYDR_TABLE <- data.frame("SITE"=as.numeric(as.character(data.frame(stack(lapply(data_hydrology, as.character)), "YEAR_SEASON"=rep(rownames(data_hydrology), times=ncol(data_hydrology))) %>% pull(ind))), 
                         "YEAR_SEASON"=data.frame(stack(lapply(data_hydrology, as.character)), "YEAR_SEASON"=rep(rownames(data_hydrology), times=ncol(data_hydrology))) %>% pull(YEAR_SEASON), 
                         "HYDR_STATE"=data.frame(stack(lapply(data_hydrology, as.character)), "YEAR_SEASON"=rep(rownames(data_hydrology), times=ncol(data_hydrology))) %>% pull(values))

# Create time series frame
TIME_SERIES_FRAME <- data.frame("YEAR"=rep(c(min(data$YEAR), rep(seq(from=min(data$YEAR)+1, to=max(data$YEAR)-1), each=4), rep(max(data$YEAR), times=2)), times=length(unique(data$SITE))), 
                                "SEASON"=rep(rep(unique(data$SEASON), times=length(unique(c(min(data$YEAR), rep(seq(from=min(data$YEAR)+1, to=max(data$YEAR)-1), each=4), rep(max(data$YEAR), times=2)))))[1:length(c(min(data$YEAR), rep(seq(from=min(data$YEAR)+1, to=max(data$YEAR)-1), each=4), rep(max(data$YEAR), times=2)))], times=length(unique(data$SITE))))
TIME_SERIES_FRAME$YEAR_SEASON <- as.factor(paste(TIME_SERIES_FRAME$YEAR, TIME_SERIES_FRAME$SEASON, sep="_"))
TIME_SERIES_FRAME$SITE <- rep(unique(data$SITE), each=length(unique(TIME_SERIES_FRAME$YEAR_SEASON)))
TIME_SERIES_FRAME <- left_join(TIME_SERIES_FRAME[,c(1,2,3,4)], data[, c(1,4,1)], by=c("SITE"="SITE", "YEAR_SEASON"="YEAR_SEASON"))
colnames(TIME_SERIES_FRAME)[ncol(TIME_SERIES_FRAME)] <- "SITE_sampled"
TIME_SERIES_FRAME$CAMPAIGN_Nr <- rep(c(1:length(levels(TIME_SERIES_FRAME$YEAR_SEASON))), times=length(unique(data$SITE)))
TIME_SERIES_FRAME$GLOBAL_Order <- c(1:nrow(TIME_SERIES_FRAME))
TIME_SERIES_FRAME <- left_join(TIME_SERIES_FRAME, HYDR_TABLE, by=join_by(SITE==SITE, YEAR_SEASON==YEAR_SEASON))
TIME_SERIES_FRAME$color <- ifelse(is.na(TIME_SERIES_FRAME$HYDR_STATE),ifelse(is.na(TIME_SERIES_FRAME$SITE_sampled),NA,as.character(TIME_SERIES_FRAME$SEASON)),ifelse(TIME_SERIES_FRAME$HYDR_STATE=="DRY"|TIME_SERIES_FRAME$HYDR_STATE=="PLS","DRY/PLS",as.character(TIME_SERIES_FRAME$SEASON)))
TIME_SERIES_FRAME$tile <- ifelse(is.na(TIME_SERIES_FRAME$color), NA, TIME_SERIES_FRAME$SITE)

campaign_full_list <- factor(unique(TIME_SERIES_FRAME$YEAR_SEASON), levels=c(unique(as.character(TIME_SERIES_FRAME$YEAR_SEASON))))
data <- cbind(data, "CAMPAIGN_Nr"=match(data$YEAR_SEASON, campaign_full_list))
data <- data %>% relocate(CAMPAIGN_Nr, .after=YEAR_SEASON)
data <- data %>% arrange(SITE, CAMPAIGN_Nr)

# Check the number of trajectory points per site
site_list <- unique(data$SITE)
campaign_list <- factor(unique(data$YEAR_SEASON), levels=c(unique(as.character(data$YEAR_SEASON))))

# Set basic values
seasons <- c("SP", "SU", "AU", "WI")
seasons_colour_list <- setNames(c("darkgreen", "orange", "brown", "darkblue", "darkgrey"), c(seasons, "DRY/PLS"))
samp_ID_ncol <- 6 # First n columns identifying samples
spec_ncol <- ncol(data)-samp_ID_ncol # Number of columns with abundance data of taxa (following samp_ID_ncol)
spec_list <- colnames(data)[(samp_ID_ncol+1):ncol(data)]
dist_metric <- "bray"

## REPLACE DATA TABLE (to remove dry and pools hydrological states)
data_new <- as.data.frame(as.data.frame(data %>% left_join(TIME_SERIES_FRAME %>% select(SITE, CAMPAIGN_Nr, color), join_by(SITE==SITE, CAMPAIGN_Nr==CAMPAIGN_Nr)) %>% filter(color %in% c("SP", "SU", "AU", "WI"))) %>% arrange(SITE, CAMPAIGN_Nr) %>% group_by(SITE) %>% mutate(SITE_ORDER=1:n()) %>% ungroup %>% relocate(SITE_ORDER, .after = CAMPAIGN_Nr) %>% select(-color))
data_new <- data_new[,-c(as.numeric(which(colSums(data_new[, (samp_ID_ncol+1):ncol(data_new)])==0))+samp_ID_ncol)]
data <- data_new

# Define stages
stages_start_list <- c(1,4,12,16) # CAMPAIGN_Nrs rigth after SPs
stages_end_list <- c(3,11,15,19) # CAMPAIGN_Nrs of SPs

# Check stage assigment start-end numbers
if(length(stages_start_list)!=length(stages_end_list)){
  stop("The number of start and end points of stages are not equal")}

# Create stage interval table
stage_assignment_table <- data.frame("Interval"=paste("int_", seq(1:length(stages_start_list)), sep=""), 
                                     "start"=stages_start_list, 
                                     "end"=stages_end_list)


## -- FULL DATA SET -- ##

# FULL TRAJECTORIES FOR THE FULL DATA SET
# available only for sites with at least 2 samples
# tr_directionality only available for sites with at least 3 samples
FULL_DATASET_FULL_TRAJECTORIES <- list()
sites <- as.data.frame(data %>% filter(SITE %in% as.numeric(as.data.frame(data %>% group_by(SITE) %>% summarise(n_CAMP=length(CAMPAIGN_Nr))) %>% filter(n_CAMP>1) %>% pull(SITE))))$SITE
surveys <- as.data.frame(data %>% filter(SITE %in% as.numeric(as.data.frame(data %>% group_by(SITE) %>% summarise(n_CAMP=length(CAMPAIGN_Nr))) %>% filter(n_CAMP>1) %>% pull(SITE))))$CAMPAIGN_Nr
sites_dir <- as.data.frame(data %>% filter(SITE %in% as.numeric(as.data.frame(data %>% group_by(SITE) %>% summarise(n_CAMP=length(CAMPAIGN_Nr))) %>% filter(n_CAMP>2) %>% pull(SITE))))$SITE
surveys_dir <- as.data.frame(data %>% filter(SITE %in% as.numeric(as.data.frame(data %>% group_by(SITE) %>% summarise(n_CAMP=length(CAMPAIGN_Nr))) %>% filter(n_CAMP>2) %>% pull(SITE))))$SITE
dev.new()
tr_def <- defineTrajectories(vegan::vegdist(as.data.frame(data %>% filter(SITE %in% as.numeric(as.data.frame(data %>% group_by(SITE) %>% summarise(n_CAMP=length(CAMPAIGN_Nr))) %>% filter(n_CAMP>1) %>% pull(SITE))))[,c((samp_ID_ncol+1):ncol(data))], method=dist_metric), sites, surveys)
tr_dir <- defineTrajectories(vegan::vegdist(as.data.frame(data %>% filter(SITE %in% as.numeric(as.data.frame(data %>% group_by(SITE) %>% summarise(n_CAMP=length(CAMPAIGN_Nr))) %>% filter(n_CAMP>2) %>% pull(SITE))))[,c((samp_ID_ncol+1):ncol(data))], method=dist_metric), sites_dir, surveys_dir)
tr <- trajectoryPCoA(tr_def)
dev.off()
FULL_DATASET_FULL_TRAJECTORIES <- append(FULL_DATASET_FULL_TRAJECTORIES, list(tr$eig))
FULL_DATASET_FULL_TRAJECTORIES <- append(FULL_DATASET_FULL_TRAJECTORIES, list(data.frame("SITE"=sites, "CAMPAIGN_Nr"=surveys, tr$points)))
FULL_DATASET_FULL_TRAJECTORIES <- append(FULL_DATASET_FULL_TRAJECTORIES, list(data.frame("SITE"=unique(sites),trajectoryLengths(tr_def)) %>% as.data.frame(row.names = 1:nrow(.))))
FULL_DATASET_FULL_TRAJECTORIES <- append(FULL_DATASET_FULL_TRAJECTORIES, list(data.frame("SITE"=unique(sites_dir), "tr_dir"=as.numeric(trajectoryDirectionality(tr_dir)))))
remove(tr, sites_dir, surveys_dir)
segment_dist_IF <- data.frame(segmentDistances(tr_def)$Dinifin)
rownames(segment_dist_IF) <- str_replace_all(rownames(segment_dist_IF), "\\[", "_")
rownames(segment_dist_IF) <- str_replace_all(rownames(segment_dist_IF), "-", "_")
rownames(segment_dist_IF) <- paste("I", str_sub(rownames(segment_dist_IF),1,-2), sep="_")
colnames(segment_dist_IF) <- paste("F", str_sub(str_replace_all(colnames(segment_dist_IF), "\\.", "_"), 2, -2), sep="_")
FULL_DATASET_FULL_TRAJECTORIES <- append(FULL_DATASET_FULL_TRAJECTORIES, list(data.frame(segment_dist_IF)))
remove(segment_dist_IF)
segment_dist_FF <- data.frame(as.matrix(segmentDistances(tr_def)$Dfin))
rownames(segment_dist_FF) <- str_replace_all(rownames(segment_dist_FF), "\\[", "_")
rownames(segment_dist_FF) <- str_replace_all(rownames(segment_dist_FF), "-", "_")
rownames(segment_dist_FF) <- paste("F", str_sub(rownames(segment_dist_FF),1,-2), sep="_")
colnames(segment_dist_FF) <- paste("F", str_sub(str_replace_all(colnames(segment_dist_FF), "\\.", "_"), 2, -2), sep="_")
FULL_DATASET_FULL_TRAJECTORIES <- append(FULL_DATASET_FULL_TRAJECTORIES, list(data.frame(segment_dist_FF)))
remove(segment_dist_FF)
names(FULL_DATASET_FULL_TRAJECTORIES) <- c("tr_eig", "tr_points", "tr_lengths", "tr_directionality", "segment_dist_IF", "segment_dist_FF")

# Create groups of sampling sites

grouping_table <- data.frame("SITE"=as.numeric(names(cutree(hclust(dist(xtabs(DD_365 ~ SITE + CAMPAIGN_Nr, as.data.frame(DD_365_FINAL %>% filter(SITE %in% as.numeric(TIME_SERIES_FRAME %>% filter(!is.na(color), color!="DRY/PLS", SEASON=="SP") %>% group_by(SITE) %>% summarise(n_SP=length(SEASON)) %>% filter(n_SP>3, SITE %in% unique(DD_365_FINAL$SITE)) %>% pull(SITE))) %>% left_join(TIME_SERIES_FRAME[,c(3,4,6)], join_by(SITE==SITE, YEAR_SEASON==YEAR_SEASON))))), method="ward.D2"), k=3))), 
                             "group"=as.numeric(cutree(hclust(dist(xtabs(DD_365 ~ SITE + CAMPAIGN_Nr, as.data.frame(DD_365_FINAL %>% filter(SITE %in% as.numeric(TIME_SERIES_FRAME %>% filter(!is.na(color), color!="DRY/PLS", SEASON=="SP") %>% group_by(SITE) %>% summarise(n_SP=length(SEASON)) %>% filter(n_SP>3, SITE %in% unique(DD_365_FINAL$SITE)) %>% pull(SITE))) %>% left_join(TIME_SERIES_FRAME[,c(3,4,6)], join_by(SITE==SITE, YEAR_SEASON==YEAR_SEASON))))), method="ward.D2"), k=3)))
grouping_table <- grouping_table %>% left_join(as.data.frame(DD_365_FINAL %>% left_join(grouping_table, join_by(SITE==SITE)) %>% filter(!is.na(group)) %>% group_by(group) %>% summarise(mean_DD365=mean(DD_365)) %>% arrange(mean_DD365) %>% mutate(dry_group=paste("DD365", strrep("+", c(seq(0,length(unique(grouping_table$group))-1))), sep="")))[,c(1,3)], join_by(group==group))
grouping_table <- grouping_table %>% left_join(as.data.frame(TIME_SERIES_FRAME %>% filter(SITE %in% c(grouping_table$SITE), color %in% c("SP", "SU", "AU", "WI")) %>% group_by(SITE) %>% summarise(n_SP=sum(ifelse(color=="SP",1,0)), n_SU=sum(ifelse(color=="SU",1,0)), n_AU=sum(ifelse(color=="AU",1,0)), n_WI=sum(ifelse(color=="WI",1,0)))) %>% arrange(n_SU) %>% mutate("SU_AU_WI"=n_SU+n_AU+n_WI) %>% arrange(desc(SU_AU_WI)) %>% mutate(exp_group=ifelse(SU_AU_WI==12,"Perennial", ifelse(SU_AU_WI>9,"Summer dry","More seasons dry"))) %>% select(SITE, exp_group), join_by(SITE==SITE))

# Correction for SITE 9
grouping_table$exp_group[which(grouping_table$SITE==9)] <- "Summer dry"

# Set factor levels
grouping_table$dry_group <- factor(grouping_table$dry_group, levels=c("DD365", "DD365+", "DD365++"))
grouping_table$exp_group <- factor(grouping_table$exp_group, levels=c("Perennial", "Summer dry", "More seasons dry"))

# ENVELOPES

# Create a list of sites and their campaign list

envelope_list <- NULL

S1_list <- NULL
S1_list <- append(S1_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,12,13,14,15,16,17,18,19), 
  c(1,1,1,1,2,2,2,2,3,3,3,3)))
names(S1_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S2_list <- NULL
S2_list <- append(S2_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,12,13,14,15,16,17,18,19), 
  c(1,1,1,1,2,2,2,2,3,3,3,3)))
names(S2_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S3_list <- NULL
S3_list <- append(S3_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,12,13,14,15,16,17,18,19), 
  c(1,1,1,1,2,2,2,2,3,3,3,3)))
names(S3_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S4_list <- NULL
S4_list <- append(S4_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,12,13,14,15,16,17,18,19), 
  c(1,1,1,1,2,2,2,2,3,3,3,3)))
names(S4_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S5_list <- NULL
S5_list <- append(S5_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,13,14,15,17,18,19), 
  c(1,1,1,1,2,2,2,3,3,3)))
names(S5_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S6_list <- NULL
S6_list <- append(S6_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,14,15,18,19), 
  c(1,1,1,1,2,2,3,3)))
names(S6_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S8_list <- NULL
S8_list <- append(S8_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,12,13,14,15,16,17,18,19), 
  c(1,1,1,1,2,2,2,2,3,3,3,3)))
names(S8_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S9_list <- NULL
S9_list <- append(S9_list,list(
  c(1,2,3,4), 
  c(9,10,11,13,14,15,17,18,19), 
  c(1,1,1,2,2,2,3,3,3)))
names(S9_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S10_list <- NULL
S10_list <- append(S10_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,12,13,14,15,17,18,19), 
  c(1,1,1,1,2,2,2,2,3,3,3)))
names(S10_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S11_list <- NULL
S11_list <- append(S11_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,12,13,14,15,16,17,18,19), 
  c(1,1,1,1,2,2,2,2,3,3,3,3)))
names(S11_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S12_list <- NULL
S12_list <- append(S12_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,12,13,14,15,16,17,18,19), 
  c(1,1,1,1,2,2,2,2,3,3,3,3)))
names(S12_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S13_list <- NULL
S13_list <- append(S13_list,list(
  c(1,2,3), 
  c(9,10,11,14,15,18,19), 
  c(1,1,1,2,2,3,3)))
names(S13_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S14_list <- NULL
S14_list <- append(S14_list,list(
  c(1,2,3,4), 
  c(9,10,11,13,14,15,16,17,18,19), 
  c(1,1,1,2,2,2,3,3,3,3)))
names(S14_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S15_list <- NULL
S15_list <- append(S15_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,12,13,14,15,16,17,18,19), 
  c(1,1,1,1,2,2,2,2,3,3,3,3)))
names(S15_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S16_list <- NULL
S16_list <- append(S16_list,list(
  c(1,2,3,4), 
  c(10,11,14,15,18,19), 
  c(1,1,2,2,3,3)))
names(S16_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S18_list <- NULL
S18_list <- append(S18_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,13,14,15,17,18,19), 
  c(1,1,1,1,2,2,2,3,3,3)))
names(S18_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S19_list <- NULL
S19_list <- append(S19_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,12,13,14,15,16,17,18,19), 
  c(1,1,1,1,2,2,2,2,3,3,3,3)))
names(S19_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S20_list <- NULL
S20_list <- append(S20_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,14,15,18,19), 
  c(1,1,1,1,2,2,3,3)))
names(S20_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S22_list <- NULL
S22_list <- append(S22_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,14,15,17,18,19), 
  c(1,1,1,1,2,2,3,3,3)))
names(S22_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S23_list <- NULL
S23_list <- append(S23_list,list(
  c(1,2,3,4), 
  c(9,10,11,14,15,17,18,19), 
  c(1,1,1,2,2,3,3,3)))
names(S23_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S24_list <- NULL
S24_list <- append(S24_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,13,14,15,17,18,19), 
  c(1,1,1,1,2,2,2,3,3,3)))
names(S24_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S30_list <- NULL
S30_list <- append(S30_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,12,13,14,15,17,18,19), 
  c(1,1,1,1,2,2,2,2,3,3,3)))
names(S30_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S32_list <- NULL
S32_list <- append(S32_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,12,13,14,15,17,18,19), 
  c(1,1,1,1,2,2,2,2,3,3,3)))
names(S32_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S35_list <- NULL
S35_list <- append(S35_list,list(
  c(1,2,3,4), 
  c(9,10,11,12,13,14,15,17,18,19), 
  c(1,1,1,2,2,2,2,3,3,3)))
names(S35_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S36_list <- NULL
S36_list <- append(S36_list,list(
  c(1,2,3,4), 
  c(10,11,14,15,18,19), 
  c(1,1,2,2,3,3)))
names(S36_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S37_list <- NULL
S37_list <- append(S37_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,12,13,14,15,17,18,19), 
  c(1,1,1,1,2,2,2,2,3,3,3)))
names(S37_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S38_list <- NULL
S38_list <- append(S38_list,list(
  c(1,2,3,4), 
  c(10,11,14,15,18,19), 
  c(1,1,2,2,3,3)))
names(S38_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S39_list <- NULL
S39_list <- append(S39_list,list(
  c(1,2,3,4), 
  c(10,11,14,15,18,19), 
  c(1,1,2,2,3,3)))
names(S39_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S40_list <- NULL
S40_list <- append(S40_list,list(
  c(2,3), 
  c(10,11,14,15,18,19), 
  c(1,1,2,2,3,3)))
names(S40_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S42_list <- NULL
S42_list <- append(S42_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,12,13,14,15,17,18,19), 
  c(1,1,1,1,2,2,2,2,3,3,3)))
names(S42_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S44_list <- NULL
S44_list <- append(S44_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,13,14,15,17,18,19), 
  c(1,1,1,1,2,2,2,3,3,3)))
names(S44_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S45_list <- NULL
S45_list <- append(S45_list,list(
  c(1,2,3,4), 
  c(8,9,10,11,13,14,15,17,18,19), 
  c(1,1,1,1,2,2,2,3,3,3)))
names(S45_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

S49_list <- NULL
S49_list <- append(S49_list,list(
  c(2,3), 
  c(10,11,14,15,18,19), 
  c(1,1,2,2,3,3)))
names(S49_list) <- c("selected_envelope", "selected_to_compare", "sequence_to_compare")

envelope_list <- append(envelope_list, list(
  S1_list, S2_list, S3_list, S4_list, S5_list, S6_list, S8_list, S9_list, S10_list, S11_list, S12_list, S13_list, S14_list, S15_list, S16_list, S18_list, S19_list, S20_list, S22_list, S23_list, S24_list, S30_list, S32_list, S35_list, S36_list, S37_list, S38_list, S39_list, S40_list, S42_list, S44_list, S45_list, S49_list
))
names(envelope_list) <- c("1", "2", "3", "4", "5", "6", "8", "9", "10", "11", "12", "13", "14", "15", "16", "18", "19", "20", "22", "23", "24", "30", "32", "35", "36", "37", "38", "39", "40", "42", "44", "45", "49")
remove(S1_list, S2_list, S3_list, S4_list, S5_list, S6_list, S8_list, S9_list, S10_list, S11_list, S12_list, S13_list, S14_list, S15_list, S16_list, S18_list, S19_list, S20_list, S22_list, S23_list, S24_list, S30_list, S32_list, S35_list, S36_list, S37_list, S38_list, S39_list, S40_list, S42_list, S44_list, S45_list, S49_list)

ecotraj_Q_FINAL_TABLE <- data.frame()

for(i in as.numeric(names(envelope_list))){
  selected_site <- i
  selected_envelope <- envelope_list[[which(names(envelope_list)==i)]]$selected_envelope
  selected_to_compare <- envelope_list[[which(names(envelope_list)==i)]]$selected_to_compare
  sequence_to_compare <- envelope_list[[which(names(envelope_list)==i)]]$sequence_to_compare
  
  data_raw <- data[which(data$SITE==selected_site),]
  
  # Check if selected campaigns are available
  if(any(is.na(match(c(selected_envelope, selected_to_compare), data_raw$CAMPAIGN_Nr)))){
    stop("Samples from one or more selected campaign(s) cannot be found in the dataset")}
  
  # Check if envelope and compare shares campaigns
  if(length(intersect(selected_envelope, selected_to_compare))>0){
    stop("There are identical samples in envelope and those to be compared to it")}
  
  # Check if the sequence and list of samples to compare are equal in length
  if(!length(selected_to_compare)==length(sequence_to_compare)){
    stop("The list and sequence of samples to compare are not equal in length")}
  
  data_raw <- data_raw[,-c(as.numeric(which(colSums(data_raw[-c(1:samp_ID_ncol)])==0))+samp_ID_ncol)]
  data_raw <- data_raw[match(c(selected_envelope, selected_to_compare), data_raw$CAMPAIGN_Nr),]
  data_raw <- cbind("NAME"=NA, data_raw)
  
  # Name and select envelope samples and samples to compare
  data_raw$NAME[match(selected_envelope, data_raw$CAMPAIGN_Nr)] <- paste("E_", seq(from=1, to=length(match(selected_envelope, data_raw$CAMPAIGN_Nr))), sep="")
  data_raw$NAME[match(data_raw$CAMPAIGN_Nr[-match(selected_envelope, data_raw$CAMPAIGN_Nr)], data_raw$CAMPAIGN_Nr)] <- paste("C", sequence_to_compare, "_", seq(from=1, to=length(match(data_raw$CAMPAIGN_Nr[-match(selected_envelope, data_raw$CAMPAIGN_Nr)], data_raw$CAMPAIGN_Nr))), sep="")
  data_raw$NAME <- paste(str_before_first(data_raw$NAME, "_"), "_", data.frame("NAME"=str_before_first(data_raw$NAME, "_")) %>% group_by(NAME) %>% mutate(counter=row_number(NAME)) %>% pull(counter), sep="")
  rownames(data_raw) <- data_raw$NAME
  envelope <- rownames(data_raw)[grep("E_", rownames(data_raw))]
  
  spec_dist <- vegan::vegdist(data_raw[,-c(1:(samp_ID_ncol+1))], method=dist_metric)
  
  ecotraj_Q_table <- data.frame("squared_dist"=compareToStateEnvelope(spec_dist, envelope, distances_to_envelope = TRUE)[,3], "ecotraj_Q"=compareToStateEnvelope(spec_dist, envelope, distances_to_envelope = TRUE)[,4], "group"=factor(str_before_first(data_raw$NAME, "_"), levels=c(unique(str_before_first(data_raw$NAME, "_")))), "NAME"=data_raw$NAME)
  
  ecotraj_Q_FINAL_TABLE <- rbind(ecotraj_Q_FINAL_TABLE, data.frame("SITE"=TIME_SERIES_FRAME[which(TIME_SERIES_FRAME$SITE==selected_site),][match(seq(from=min((c(selected_envelope, selected_to_compare))), to=max((c(selected_envelope, selected_to_compare)))), TIME_SERIES_FRAME[which(TIME_SERIES_FRAME$SITE==selected_site),]$CAMPAIGN_Nr),]$SITE, 
                                                                   "CAMPAIGN_Nr"=TIME_SERIES_FRAME[which(TIME_SERIES_FRAME$SITE==selected_site),][match(seq(from=min((c(selected_envelope, selected_to_compare))), to=max((c(selected_envelope, selected_to_compare)))), TIME_SERIES_FRAME[which(TIME_SERIES_FRAME$SITE==selected_site),]$CAMPAIGN_Nr),]$CAMPAIGN_Nr, 
                                                                   "HDR_STATE"=TIME_SERIES_FRAME[which(TIME_SERIES_FRAME$SITE==selected_site),][match(seq(from=min((c(selected_envelope, selected_to_compare))), to=max((c(selected_envelope, selected_to_compare)))), TIME_SERIES_FRAME[which(TIME_SERIES_FRAME$SITE==selected_site),]$CAMPAIGN_Nr),]$HYDR_STATE) %>% 
                                   left_join(data_raw[,c(2,6,1)], join_by(SITE==SITE, CAMPAIGN_Nr==CAMPAIGN_Nr)) %>%
                                   left_join(ecotraj_Q_table, join_by(NAME==NAME)) %>%
                                   left_join(as.data.frame(aggregate(data.frame("group"=ecotraj_Q_table$group, "mean_Q"=ecotraj_Q_table$ecotraj_Q), mean_Q~group, FUN=mean)), join_by(group==group)) %>%
                                   left_join(as.data.frame(aggregate(data.frame("group"=ecotraj_Q_table$group, "mean_dist"=ecotraj_Q_table$squared_dist), mean_dist~group, FUN=mean)), join_by(group==group)) %>%  
                                   mutate("E_start"=min(selected_envelope), "E_end"=max(selected_envelope)))
}

remove(data_raw, spec_dist, envelope, selected_site, selected_envelope, selected_to_compare, sequence_to_compare, sites, surveys, ecotraj_Q_table)

ecotraj_Q_FINAL_TABLE <- ecotraj_Q_FINAL_TABLE %>% left_join(grouping_table[,c(1,3)], join_by(SITE==SITE))

## FIG 7 ##

fig_7_plotlist <- list()

# Boxplot of mean squared distances from envelopes in consecutive periods within How long? groups
p <- ggplot(na.omit(ecotraj_Q_FINAL_TABLE), aes(x=factor(dry_group, levels=c("DD365", "DD365+", "DD365++")), y=squared_dist, fill=dry_group, alpha=group)) + 
  scale_fill_viridis(option="D", direction=1, discrete=TRUE) + 
  scale_y_continuous(limits=c(0,1.05), breaks=seq(0,1, by=0.2)) +
  scale_x_discrete(labels=c("DD365\nQuasi-perennial", "DD365+\nShort-term drying", "DD365++\nLong-term drying")) +
  scale_alpha_discrete(range=c(0,1), labels=c("E", "C1", "C2", "C3"),
                       guide=guide_legend(override.aes = list(fill="black", alpha=c(0,0.2,0.4,0.9)))) +
  geom_boxplot(linewidth=0.8, median.linewidth=1.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(size=13, color="black")) +
  xlab("") +
  ylab("Squared distance from state envelope") +
  labs(title="", 
       subtitle = "E: envelope, Cx: consecutive periods after envelope",
       fill="Duration of drying       ", alpha="Envelope/Period")
fig_7_plotlist <- append(fig_7_plotlist, list(p))
names(fig_7_plotlist)[length(fig_7_plotlist)] <- "7E"

# Boxplot of mean squared distances from envelopes in consecutive periods within When? groups
p <- ggplot(na.omit(ecotraj_Q_FINAL_TABLE %>% left_join(grouping_table[,c(1,4)], join_by(SITE==SITE)) %>% mutate(exp_group=factor(exp_group, levels=c("Perennial", "Summer dry", "More seasons dry")))), aes(x=factor(exp_group, levels=c("Perennial", "Summer dry", "More seasons dry")), y=squared_dist, fill=exp_group, alpha=group)) + 
  scale_fill_viridis(option="E", direction=1, discrete=TRUE) + 
  scale_y_continuous(limits=c(0,1.05), breaks=seq(0,1, by=0.2)) +
  scale_x_discrete(labels=c("Perennial", "Summer-dry", "More-seasons-dry")) +
  scale_alpha_discrete(range=c(0,1), labels=c("E", "C1", "C2", "C3"),
                       guide=guide_legend(override.aes = list(fill="black", alpha=c(0,0.2,0.4,0.9)))) +
  geom_boxplot(linewidth=0.8, median.linewidth=1.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(size=13, color="black")) +
  xlab("") +
  ylab("Squared distance from state envelope") +
  labs(title="", 
       subtitle = "E: envelope, Cx: consecutive periods after envelope",
       fill="Seasonal occurrence \nof drying", alpha="Envelope/Period")
fig_7_plotlist <- append(fig_7_plotlist, list(p))
names(fig_7_plotlist)[length(fig_7_plotlist)] <- "7A"

# Compare mean squared distances from envelopes in consecutive periods within When? groups
mixed_model_data <- as.data.frame(ecotraj_Q_FINAL_TABLE %>% left_join(grouping_table[,c(1,4)], join_by(SITE==SITE)) %>% mutate(exp_group=factor(exp_group, levels=c("Perennial", "Summer dry", "More seasons dry"))) %>% filter(grepl("C", group)))

res.anova <- as.data.frame(anova(lmer(squared_dist ~ exp_group + group + (1|SITE), mixed_model_data)))
names(res.anova) <- gsub(" ", "_", names(res.anova))
res.anova <- res.anova %>% rename(p="Pr(>F)") %>% mutate(DenDF=round(DenDF,2), F_value=round(F_value, 2), p=ifelse(p<0.001, paste("<0.001"), round(p,3)))

pwc <- as.data.frame(contrast(emmeans(lmer(squared_dist ~ exp_group + group + (1|SITE), mixed_model_data), ~ exp_group), method="pairwise", adjust="tukey"))
pwc <- pwc %>% mutate(group1=str_before_first(contrast," "), 
                      group2=str_after_first(contrast, " - "),
                      p.adj=p.value,
                      p.adj.signif=symnum(p.adj, corr=FALSE, na=FALSE, cutpoints=c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols=c("****", "***", "**", "*", "ns")),
                      y.position=as.numeric(mixed_model_data %>% tukey_hsd(squared_dist ~ exp_group) %>% add_xy_position(x="exp_group") %>% pull(y.position)),
                      xmin=as.numeric(mixed_model_data %>% tukey_hsd(squared_dist ~ exp_group) %>% add_xy_position(x="exp_group") %>% pull(xmin)),
                      xmax=as.numeric(mixed_model_data %>% tukey_hsd(squared_dist ~ exp_group) %>% add_xy_position(x="exp_group") %>% pull(xmax)))

p <- ggplot(mixed_model_data, aes(x=exp_group, y=squared_dist)) + 
  geom_boxplot(aes(fill=exp_group), linewidth=0.8, median.linewidth=1.5, show.legend = FALSE) +
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(na.omit(mixed_model_data$squared_dist))+0.05, by=0.06, length.out=length(y.position))), aes(y=y.position, x=xmin, xend=xmax)) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(na.omit(mixed_model_data$squared_dist))+0.05, by=0.06, length.out=length(y.position))), aes(x=xmin, y=y.position, yend=y.position-0.01)) +
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(na.omit(mixed_model_data$squared_dist))+0.05, by=0.06, length.out=length(y.position))), aes(x=xmax, y=y.position, yend=y.position-0.01)) +
  geom_text(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(na.omit(mixed_model_data)$squared_dist)+0.05, by=0.06, length.out=length(y.position))), aes(label=as.character(p.adj.signif), y=y.position, x=(xmin+xmax)/2,), position=position_nudge(y=0.02)) +
  scale_y_continuous(limits=c(0,1.05), breaks=seq(0,1, by=0.2)) +
  scale_x_discrete(labels=c("Perennial", "Summer-dry", "More-seasons-dry")) +
  scale_fill_viridis(option="E", discrete=TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(margin=margin(0,0,2,0)),
        plot.subtitle = element_text(margin=margin(0,0,5,0)),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(size=12, color="black")) +
  xlab("") +
  ylab("Squared distance from state envelope") +
  labs(title="",
       subtitle=paste("Mixed model ANOVA, ", "F(",res.anova$NumDF[1], ",", res.anova$DenDF[1], ")","=",res.anova$F_value[1], ", ", "p=", res.anova$p[1], sep=""))
fig_7_plotlist <- append(fig_7_plotlist, list(p))
names(fig_7_plotlist)[length(fig_7_plotlist)] <- "7B"

# Compare mean squared distances from envelopes in consecutive periods within How long? groups
res.anova <- as.data.frame(anova(lmer(squared_dist ~ dry_group + group + (1|SITE), mixed_model_data)))
names(res.anova) <- gsub(" ", "_", names(res.anova))
res.anova <- res.anova %>% rename(p="Pr(>F)") %>% mutate(DenDF=round(DenDF,2), F_value=round(F_value, 2), p=ifelse(p<0.001, paste("<0.001"), round(p,3)))

pwc <- as.data.frame(contrast(emmeans(lmer(squared_dist ~ dry_group + group + (1|SITE), mixed_model_data), ~ dry_group), method="pairwise", adjust="tukey"))
pwc <- pwc %>% mutate(group1=str_before_first(contrast," "), 
                      group2=str_after_first(contrast, " - "),
                      p.adj=p.value,
                      p.adj.signif=symnum(p.adj, corr=FALSE, na=FALSE, cutpoints=c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols=c("****", "***", "**", "*", "ns")),
                      y.position=as.numeric(mixed_model_data %>% tukey_hsd(squared_dist ~ dry_group) %>% add_xy_position(x="dry_group") %>% pull(y.position)),
                      xmin=as.numeric(mixed_model_data %>% tukey_hsd(squared_dist ~ dry_group) %>% add_xy_position(x="dry_group") %>% pull(xmin)),
                      xmax=as.numeric(mixed_model_data %>% tukey_hsd(squared_dist ~ dry_group) %>% add_xy_position(x="dry_group") %>% pull(xmax)))

p <- ggplot(mixed_model_data, aes(x=dry_group, y=squared_dist)) + 
  geom_boxplot(aes(fill=dry_group), linewidth=0.8, median.linewidth=1.5, show.legend = FALSE) +
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(na.omit(mixed_model_data$squared_dist))+0.05, by=0.06, length.out=length(y.position))), aes(y=y.position, x=xmin, xend=xmax)) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(na.omit(mixed_model_data$squared_dist))+0.05, by=0.06, length.out=length(y.position))), aes(x=xmin, y=y.position, yend=y.position-0.01)) +
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(na.omit(mixed_model_data$squared_dist))+0.05, by=0.06, length.out=length(y.position))), aes(x=xmax, y=y.position, yend=y.position-0.01)) +
  geom_text(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(na.omit(mixed_model_data)$squared_dist)+0.05, by=0.06, length.out=length(y.position))), aes(label=as.character(p.adj.signif), y=y.position, x=(xmin+xmax)/2,), position=position_nudge(y=0.02)) +
  scale_y_continuous(limits=c(0,1.05), breaks=seq(0,1, by=0.2)) +
  scale_x_discrete(labels=c("DD365\nQuasi-perennial", "DD365+\nShort-term drying", "DD365++\nLong-term drying")) +
  scale_fill_viridis(option="D", discrete=TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(margin=margin(0,0,2,0)),
        plot.subtitle = element_text(margin=margin(0,0,5,0)),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(size=12, color="black")) +
  xlab("") +
  ylab("Squared distance from state envelope") +
  labs(title="",
       subtitle=paste("Mixed model ANOVA, ", "F(",res.anova$NumDF[1], ",", res.anova$DenDF[1], ")","=",res.anova$F_value[1], ", ", "p=", res.anova$p[1], sep=""))
fig_7_plotlist <- append(fig_7_plotlist, list(p))
names(fig_7_plotlist)[length(fig_7_plotlist)] <- "7F"

# Comparison of site trajectories and trajectory envelopes instead of point envelopes (How long?)
data_raw <- data[which(!is.na(match(data$SITE, names(envelope_list)))),c(1,5,c(7:ncol(data)))]
data_raw <- data_raw %>% left_join(grouping_table[,c(1,3)], join_by(SITE==SITE))
data_raw <- data_raw %>% relocate(dry_group, .after=SITE)
data_raw$dry_group <- factor(data_raw$dry_group, levels=c("DD365", "DD365+", "DD365++"))
data_raw <- data_raw[,c(1:3,as.numeric(which(colSums(data_raw[, c(4:ncol(data_raw))])>0))+3)]

envelope_results <- as.data.frame(compareToTrajectoryEnvelope(vegdist(data_raw[,c(4:ncol(data_raw))], method=dist_metric), sites=data_raw$SITE, envelope=unique(data_raw$SITE[which(data_raw$dry_group=="DD365")]), surveys=data_raw$CAMPAIGN_Nr, distances_to_envelope=TRUE))
envelope_results <- envelope_results %>% mutate(order=match(envelope_results$Site,c(unique(data_raw$SITE[which(data_raw$dry_group=="DD365")]), unique(data_raw$SITE[which(data_raw$dry_group=="DD365+")]), unique(data_raw$SITE[which(data_raw$dry_group=="DD365++")])))) %>% arrange(order) %>% select(1,2,3,4)
colnames(envelope_results)[1] <- "SITE"

envelope_results <- envelope_results %>% left_join(grouping_table[,c(1,3)], join_by(SITE==SITE))
envelope_results$dry_group <- factor(envelope_results$dry_group, levels=c("DD365", "DD365+", "DD365++"))

# Comapre squared distances from envelope (DRY_GROUP)
res.anova <- envelope_results %>% anova_test(SquaredDist ~ dry_group)
pwc <- tukey_hsd(envelope_results, SquaredDist ~ dry_group)
pwc <- pwc %>% add_xy_position(x = "dry_group")

pwc$xmin <- as.numeric(sapply(pwc$group1, switch, "DD365"=1, "DD365+"=2, "DD365++"=3))
pwc$xmax <- as.numeric(sapply(pwc$group2, switch, "DD365"=1, "DD365+"=2, "DD365++"=3))

p <- ggplot(envelope_results, aes(x=dry_group, y=SquaredDist)) + 
  geom_boxplot(aes(fill=dry_group), linewidth=0.8, median.linewidth=1.5, show.legend = FALSE) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=0.95+0.05, by=0.04, length.out=length(y.position))), aes(y=y.position, x=xmin, xend=xmax)) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=0.95+0.05, by=0.04, length.out=length(y.position))), aes(x=xmin, y=y.position, yend=y.position-0.01)) +
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=0.95+0.05, by=0.04, length.out=length(y.position))), aes(x=xmax, y=y.position, yend=y.position-0.01)) +
  geom_text(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=0.95+0.05, by=0.04, length.out=length(y.position))), aes(label=p.adj.signif, y=y.position, x=(xmin+xmax)/2,), position=position_nudge(y=0.02)) +
  scale_y_continuous(limits=c(0,1.05), breaks=seq(0,1, by=0.2)) +
  scale_x_discrete(labels=c("DD365\nQuasi-perennial\n(as Envelope)", "DD365+\nShort-term drying", "DD365++\nLong-term drying")) +
  scale_fill_viridis(option="D", discrete=TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(margin=margin(0,0,2,0)),
        plot.subtitle = element_text(margin=margin(0,0,0,0)),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(size=12, color="black")) +
  xlab("") +
  ylab("Squared dist. from trajectory envelope") +
  labs(title="",
       subtitle=get_test_label(res.anova, detailed = TRUE))
fig_7_plotlist <- append(fig_7_plotlist, list(p))
names(fig_7_plotlist)[length(fig_7_plotlist)] <- "7H"

# Roseplot of trajectory envelople distances (dry_group)
p <- ggplot(envelope_results, aes(x=reorder(factor(SITE), SquaredDist), y=SquaredDist, fill=dry_group)) +
  geom_col(data=envelope_results %>% filter(dry_group=="DD365"), color="darkgrey") +
  geom_col(data=envelope_results %>% filter(dry_group=="DD365+"), color="darkgrey") +
  geom_col(data=envelope_results %>% filter(dry_group=="DD365++"), color="darkgrey") +
  scale_fill_manual(values=c(viridis(3)), labels=c("DD365 Quasi-perennial (as Envelope)", "DD365+ Short-term drying", "DD365++ Long-term drying")) +
  coord_polar() +
  theme_bw() +
  theme(axis.title.y = element_text(size=12),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12)) +
  labs(title="", x="", y="Squared dist. from trajectory envelope", fill="Duration of drying")
fig_7_plotlist <- append(fig_7_plotlist, list(p))
names(fig_7_plotlist)[length(fig_7_plotlist)] <- "7G"

# Comparison of site trajectories and trajectory envelopes instead of point envelopes (When?)
data_raw <- data[which(!is.na(match(data$SITE, names(envelope_list)))),c(1,5,c(7:ncol(data)))]
data_raw <- data_raw %>% left_join(grouping_table[,c(1,4)], join_by(SITE==SITE))
data_raw <- data_raw %>% relocate(exp_group, .after=SITE)
data_raw$exp_group <- factor(data_raw$exp_group, levels=c("Perennial", "Summer dry", "More seasons dry"))
data_raw <- data_raw[,c(1:3,as.numeric(which(colSums(data_raw[, c(4:ncol(data_raw))])>0))+3)]

envelope_results <- as.data.frame(compareToTrajectoryEnvelope(vegdist(data_raw[,c(4:ncol(data_raw))], method=dist_metric), sites=data_raw$SITE, envelope=unique(data_raw$SITE[which(data_raw$exp_group=="Perennial")]), surveys=data_raw$CAMPAIGN_Nr, distances_to_envelope=TRUE))
envelope_results <- envelope_results %>% mutate(order=match(envelope_results$Site,c(unique(data_raw$SITE[which(data_raw$exp_group=="Perennial")]), unique(data_raw$SITE[which(data_raw$exp_group=="Summer dry")]), unique(data_raw$SITE[which(data_raw$exp_group=="More seasons dry")])))) %>% arrange(order) %>% select(1,2,3,4)
colnames(envelope_results)[1] <- "SITE"

envelope_results <- envelope_results %>% left_join(grouping_table[,c(1,4)], join_by(SITE==SITE))
envelope_results$exp_group <- factor(envelope_results$exp_group, levels=c("Perennial", "Summer dry", "More seasons dry"))

# Comapre squared distances from envelope (EXP_GROUP)
res.anova <- envelope_results %>% anova_test(SquaredDist ~ exp_group)
pwc <- tukey_hsd(envelope_results, SquaredDist ~ exp_group)
pwc <- pwc %>% add_xy_position(x = "exp_group")

pwc$xmin <- as.numeric(sapply(pwc$group1, switch, "Perennial"=1, "Summer dry"=2, "More seasons dry"=3))
pwc$xmax <- as.numeric(sapply(pwc$group2, switch, "Perennial"=1, "Summer dry"=2, "More seasons dry"=3))

p <- ggplot(envelope_results, aes(x=exp_group, y=SquaredDist)) + 
  geom_boxplot(aes(fill=exp_group), linewidth=0.8, median.linewidth=1.5, show.legend = FALSE) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=0.85+0.05, by=0.06, length.out=length(y.position))), aes(y=y.position, x=xmin, xend=xmax)) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=0.85+0.05, by=0.06, length.out=length(y.position))), aes(x=xmin, y=y.position, yend=y.position-0.01)) +
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=0.85+0.05, by=0.06, length.out=length(y.position))), aes(x=xmax, y=y.position, yend=y.position-0.01)) +
  geom_text(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=0.85+0.05, by=0.06, length.out=length(y.position))), aes(label=p.adj.signif, y=y.position, x=(xmin+xmax)/2,), position=position_nudge(y=0.02)) +
  scale_y_continuous(limits=c(0,1.05), breaks=seq(0,1, by=0.2)) +
  scale_x_discrete(labels=c("Perennial\n(as Envelope)", "Summer-dry", "More-seasons-dry")) +
  scale_fill_viridis(option="E", discrete=TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(margin=margin(0,0,2,0)),
        plot.subtitle = element_text(margin=margin(0,0,0,0)),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(size=12, color="black")) +
  xlab("") +
  ylab("Squared dist. from trajectory envelope") +
  labs(title="",
       subtitle=get_test_label(res.anova, detailed = TRUE))
fig_7_plotlist <- append(fig_7_plotlist, list(p))
names(fig_7_plotlist)[length(fig_7_plotlist)] <- "7D"

# Roseplot of trajectory envelople distances (When?)
p <- ggplot(envelope_results, aes(x=reorder(factor(SITE), SquaredDist), y=SquaredDist, fill=exp_group)) +
  geom_col(data=envelope_results %>% filter(exp_group=="Perennial"), color="darkgrey") +
  geom_col(data=envelope_results %>% filter(exp_group=="Summer dry"), color="darkgrey") +
  geom_col(data=envelope_results %>% filter(exp_group=="More seasons dry"), color="darkgrey") +
  scale_fill_manual(values=c(viridis(3, option="E")), labels=c("Perennial (as Envelope)", "Summer-dry", "More-seasons-dry                                ")) +
  coord_polar() +
  theme_bw() +
  theme(axis.title.y = element_text(size=12),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12)) +
  labs(title="", x="", y="Squared dist. from trajectory envelope", fill="Seasonal occurrence of drying")
fig_7_plotlist <- append(fig_7_plotlist, list(p))
names(fig_7_plotlist)[length(fig_7_plotlist)] <- "7C"

fig_7_plotlist <- fig_7_plotlist[order(names(fig_7_plotlist))]
ggarrange(plotlist=fig_7_plotlist, nrow=4, ncol=2, labels=c("A", "B", "C", "D", "E", "F", "G", "H"), font.label = list(size=25))
#ggsave("REV_FINAL_Fig_7.png", height = 16, width=16, bg="white")

# SPRING-SPRING LENGTHS, DISTANCES, AVERAGE SEGMENT LEGHTS AND DISTANCES FROM INTERVAL CENTROIDS 
# only available for sites with 4 SEASON="SP" samplings

SP_SP_dist_table <- data.frame(data %>% filter(SITE %in% unique(grouping_table$SITE)) %>% select(SITE, YEAR_SEASON, CAMPAIGN_Nr, SEASON) %>% filter(SEASON=="SP"))
SP_SP_dist_table <- SP_SP_dist_table %>% 
  mutate(from=c(NA,as.numeric(data %>% filter(SITE %in% unique(grouping_table$SITE)) %>% select(SITE, CAMPAIGN_Nr, SEASON, SITE_ORDER) %>% filter(SEASON=="SP") %>% pull(CAMPAIGN_Nr))[1:(length(as.numeric(data %>% filter(SITE %in% unique(grouping_table$SITE)) %>% select(SITE, CAMPAIGN_Nr, SEASON, SITE_ORDER) %>% filter(SEASON=="SP") %>% pull(CAMPAIGN_Nr)))-1)]), to=CAMPAIGN_Nr) %>% filter(to!=3)
SP_SP_dist_table <- SP_SP_dist_table %>% left_join(data[,c(1,5,6)], join_by(SITE==SITE, from==CAMPAIGN_Nr))

SP_SP_dist_list <- NULL
for(i in c(1:nrow(SP_SP_dist_table))){
  SP_SP_dist_list <- append(SP_SP_dist_list, as.numeric(FULL_DATASET_FULL_TRAJECTORIES$segment_dist_IF[grep(paste("I_", SP_SP_dist_table[i,]$SITE, "_",SP_SP_dist_table[i,]$from,"_",  sep=""), rownames(FULL_DATASET_FULL_TRAJECTORIES$segment_dist_IF)),match(colnames(FULL_DATASET_FULL_TRAJECTORIES$segment_dist_IF)[grep(paste("F_", SP_SP_dist_table[i,]$SITE, "_", sep=""), colnames(FULL_DATASET_FULL_TRAJECTORIES$segment_dist_IF))][match(SP_SP_dist_table[i,]$to, str_after_last(colnames(FULL_DATASET_FULL_TRAJECTORIES$segment_dist_IF)[grep(paste("F_", SP_SP_dist_table[i,]$SITE, "_", sep=""), colnames(FULL_DATASET_FULL_TRAJECTORIES$segment_dist_IF))], "_"))],colnames(FULL_DATASET_FULL_TRAJECTORIES$segment_dist_IF))]))
}

SP_SP_dist_table$SP_SP_dist <- SP_SP_dist_list
SP_SP_dist_table <- SP_SP_dist_table %>% mutate(order_from=SITE_ORDER, order_to=c(as.numeric(SP_SP_dist_table %>% mutate(order_from=SITE_ORDER) %>% pull(order_from))[2:nrow(SP_SP_dist_table)]-1,0))

for (i in which(SP_SP_dist_table$order_to<SP_SP_dist_table$order_from)){
  SP_SP_dist_table$order_to[i] <- max(as.numeric(str_after_first(colnames(FULL_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))[which(!is.na(FULL_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE)))][-c(1, length(colnames(FULL_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))[which(!is.na(FULL_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE)))]))],"S")))
}

SP_SP_length_list <- NULL
SP_SP_rel_length_list <- NULL
for(i in c(1:nrow(SP_SP_dist_table))){
  SP_SP_length_list <- append(SP_SP_length_list, sum(as.numeric(as.data.frame(FULL_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))[,c(match(paste("S", SP_SP_dist_table[i,]$order_from, sep=""), colnames(FULL_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))):match(paste("S", SP_SP_dist_table[i,]$order_to, sep=""), colnames(FULL_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))))])))
  SP_SP_rel_length_list <- append(SP_SP_rel_length_list, sum(as.numeric(as.data.frame(FULL_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))[,c(match(paste("S", SP_SP_dist_table[i,]$order_from, sep=""), colnames(FULL_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))):match(paste("S", SP_SP_dist_table[i,]$order_to, sep=""), colnames(FULL_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))))]))/length(as.numeric(as.data.frame(FULL_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))[,c(match(paste("S", SP_SP_dist_table[i,]$order_from, sep=""), colnames(FULL_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))):match(paste("S", SP_SP_dist_table[i,]$order_to, sep=""), colnames(FULL_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))))]))) 
}
SP_SP_dist_table$SP_SP_length <- SP_SP_length_list
SP_SP_dist_table$SP_SP_rel_length <- SP_SP_rel_length_list

# Calculate SP-SP distance from centroid
SP_SP_DIST_CENT_TABLE <- data.frame()
for (i in unique(ecotraj_Q_FINAL_TABLE$SITE)){
  dist_cent_list <- NULL
  for (j in c(1:3)){
    dev.new()
    tr_points <- data.frame(trajectoryPCoA(defineTrajectories(vegdist(as.data.frame(data %>% filter(SITE==i, between(CAMPAIGN_Nr, stage_assignment_table$end[0+j], stage_assignment_table$end[1+j])))[,c((samp_ID_ncol+1):ncol(data))][,as.numeric(which(colSums(as.data.frame(data %>% filter(SITE==i, between(CAMPAIGN_Nr, stage_assignment_table$end[0+j], stage_assignment_table$end[1+j])))[,c((samp_ID_ncol+1):ncol(data))])>0))], method=dist_metric), sites=as.numeric(data %>% filter(SITE==i, between(CAMPAIGN_Nr, stage_assignment_table$end[0+j], stage_assignment_table$end[1+j])) %>% pull(SITE)), surveys=as.numeric(data %>% filter(SITE==i, between(CAMPAIGN_Nr, stage_assignment_table$end[0+j], stage_assignment_table$end[1+j])) %>% pull(CAMPAIGN_Nr))))$points)
    dev.off()
    for (k in c(1:ncol(tr_points))){
      tr_points <- cbind(tr_points, rep(mean(tr_points[,k]), times=nrow(tr_points)))
      colnames(tr_points)[ncol(tr_points)] <- paste("C", k, sep="")
    }
    dist_cent_list <- append(dist_cent_list, mean(sqrt(rowSums(((tr_points[,grep("X", colnames(tr_points))])-(tr_points[,grep("C", colnames(tr_points))]))^2))))
  }
  SP_SP_DIST_CENT_TABLE <- rbind(SP_SP_DIST_CENT_TABLE, data.frame("SITE"=i, "Interval"=paste("int_", c(1:3), sep=""), "dist_from_cent"=dist_cent_list))
}

SP_SP_dist_table$dist_from_int_cent <- SP_SP_DIST_CENT_TABLE$dist_from_cent

SP_SP_dist_table <- SP_SP_dist_table %>% left_join(grouping_table[,c(1,3)], join_by(SITE==SITE))

SP_SP_dist_table <- SP_SP_dist_table %>% reframe(SITE=SITE, 
                                                 YEAR_SEASON=YEAR_SEASON, 
                                                 SEASON=SEASON, 
                                                 CAMPAIGN_Nr=CAMPAIGN_Nr, 
                                                 SITE_ORDER=SITE_ORDER, 
                                                 dist_from_prev_SP=SP_SP_dist, 
                                                 length_from_prev_SP=SP_SP_length, 
                                                 rel_length_from_prev_SP=SP_SP_rel_length,
                                                 dist_from_int_cent=dist_from_int_cent,
                                                 dry_group=dry_group)

# Setting up a table to COMPARE AVERAGE SEGMENT LENGTH WITH REDUCED TRAJ LENGTHS
avg_length_full_length_comp <- data.frame("SITE"=SP_SP_dist_table$SITE, "YEAR_SEASON"=SP_SP_dist_table$YEAR_SEASON,"avg_length"=SP_SP_dist_table$rel_length_from_prev_SP)

## FIG 5 & 6 ##

fig_5_plotlist <- list()
fig_6_plotlist <- list()

# Compare SP-SP average segment lengths (How long?)
mixed_model_data <- na.omit(SP_SP_dist_table) %>% left_join(grouping_table[,c(1,4)], join_by(SITE==SITE))

res.anova <- as.data.frame(anova(lmer(rel_length_from_prev_SP ~ dry_group + YEAR_SEASON + (1|SITE), mixed_model_data)))
names(res.anova) <- gsub(" ", "_", names(res.anova))
res.anova <- res.anova %>% rename(p="Pr(>F)") %>% mutate(DenDF=round(DenDF,2), F_value=round(F_value, 2), p=ifelse(p<0.001, paste("<0.001"), round(p,3)))

pwc <- as.data.frame(contrast(emmeans(lmer(rel_length_from_prev_SP ~ dry_group + YEAR_SEASON + (1|SITE), mixed_model_data), ~ dry_group), method="pairwise", adjust="tukey"))
pwc <- pwc %>% mutate(group1=str_before_first(contrast," "), 
                      group2=str_after_first(contrast, " - "),
                      p.adj=p.value,
                      p.adj.signif=symnum(p.adj, corr=FALSE, na=FALSE, cutpoints=c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols=c("****", "***", "**", "*", "ns")),
                      y.position=as.numeric(mixed_model_data %>% tukey_hsd(rel_length_from_prev_SP ~ dry_group) %>% add_xy_position(x="dry_group") %>% pull(y.position)),
                      xmin=as.numeric(mixed_model_data %>% tukey_hsd(rel_length_from_prev_SP ~ dry_group) %>% add_xy_position(x="dry_group") %>% pull(xmin)),
                      xmax=as.numeric(mixed_model_data %>% tukey_hsd(rel_length_from_prev_SP ~ dry_group) %>% add_xy_position(x="dry_group") %>% pull(xmax)))

pwc$xmin <- as.numeric(sapply(pwc$group1, switch, "DD365"=1, "(DD365+)"=2, "(DD365++)"=3))
pwc$xmax <- as.numeric(sapply(pwc$group2, switch, "DD365"=1, "(DD365+)"=2, "(DD365++)"=3))

p <- ggplot(mixed_model_data, aes(x=factor(dry_group, levels=c("DD365", "DD365+", "DD365++")), y=rel_length_from_prev_SP)) + 
  geom_boxplot(aes(fill=dry_group), linewidth=0.8, median.linewidth=1.5, show.legend = FALSE) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$rel_length_from_prev_SP)+0.05, by=0.08, length.out=length(y.position))), aes(y=y.position, x=xmin, xend=xmax)) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$rel_length_from_prev_SP)+0.05, by=0.08, length.out=length(y.position))), aes(x=xmin, y=y.position, yend=y.position-0.01)) +
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$rel_length_from_prev_SP)+0.05, by=0.08, length.out=length(y.position))), aes(x=xmax, y=y.position, yend=y.position-0.01)) +
  geom_text(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$rel_length_from_prev_SP)+0.05, by=0.08, length.out=length(y.position))), aes(label=as.character(p.adj.signif), y=y.position, x=(xmin+xmax)/2,), position=position_nudge(y=0.02)) +
  scale_y_continuous(limits=c(0,1.2), breaks=seq(0,1,by=0.2)) +
  scale_x_discrete(labels=c("DD365\nQuasi-perennial", "DD365+\nShort-term drying", "DD365++\nLong-term drying")) +
  scale_fill_viridis(option="D", discrete=TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(margin=margin(0,0,2,0)),
        plot.subtitle = element_text(margin=margin(0,0,5,0)),
        axis.title.y = element_text(size=13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(size=13, color="black")) +
  xlab("") +
  ylab("Average segment length") +
  labs(title="", 
       subtitle=paste("Mixed model ANOVA, ", "F(",res.anova$NumDF[1], ",", res.anova$DenDF[1], ")","=",res.anova$F_value[1], ", ", "p=", res.anova$p[1], sep=""))
fig_5_plotlist <- append(fig_5_plotlist, list(p))
names(fig_5_plotlist)[length(fig_5_plotlist)] <- "5D"
remove(pwc, res.anova)

# Roseplot SP-SP average segment lengths (How long?)
p <- ggplot(SP_SP_dist_table %>% left_join(grouping_table[,c(1,3)], join_by(SITE==SITE)), aes(x=reorder(factor(SITE), rel_length_from_prev_SP), y=rel_length_from_prev_SP, alpha=forcats::fct_rev(factor(YEAR_SEASON)), fill=factor(dry_group))) + 
  geom_col(data=SP_SP_dist_table %>% left_join(grouping_table[,c(1,2)], join_by(SITE==SITE)) %>% filter(group==1), color="darkgrey") +
  geom_col(data=SP_SP_dist_table %>% left_join(grouping_table[,c(1,2)], join_by(SITE==SITE)) %>% filter(group==2), color="darkgrey") +
  geom_col(data=SP_SP_dist_table %>% left_join(grouping_table[,c(1,2)], join_by(SITE==SITE)) %>% filter(group==3), color="darkgrey") +
  scale_alpha_discrete(range=c(0.9,0.3), labels=c("2022-2023", "2021-2022", "2019-2021"),
                       guide=guide_legend(override.aes = list(fill="black", alpha=c(0.9,0.4,0.2)))) +
  scale_fill_manual(values=c(viridis(3)), labels=c("DD365 Quasi-perennial", "DD365+ Short-term drying", "DD365++ Long-term drying")) +
  coord_polar() +
  theme_bw() +
  theme(axis.title.y = element_text(size=13),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12)) +
  labs(title="", x="", y="Average segment length", fill="Duration of drying", alpha="Spring-to-spring interval")
fig_5_plotlist <- append(fig_5_plotlist, list(p))
names(fig_5_plotlist)[length(fig_5_plotlist)] <- "5C"

# Compare SP-SP average segment lengths (When?)
res.anova <- as.data.frame(anova(lmer(rel_length_from_prev_SP ~ exp_group + YEAR_SEASON + (1|SITE), mixed_model_data)))
names(res.anova) <- gsub(" ", "_", names(res.anova))
res.anova <- res.anova %>% rename(p="Pr(>F)") %>% mutate(DenDF=round(DenDF,2), F_value=round(F_value, 2), p=ifelse(p<0.001, paste("<0.001"), round(p,3)))

pwc <- as.data.frame(contrast(emmeans(lmer(rel_length_from_prev_SP ~ exp_group + YEAR_SEASON + (1|SITE), mixed_model_data), ~ exp_group), method="pairwise", adjust="tukey"))
pwc <- pwc %>% mutate(group1=str_before_first(contrast," -"), 
                      group2=str_after_first(contrast, " - "),
                      p.adj=p.value,
                      p.adj.signif=symnum(p.adj, corr=FALSE, na=FALSE, cutpoints=c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols=c("****", "***", "**", "*", "ns")),
                      y.position=as.numeric(mixed_model_data %>% tukey_hsd(rel_length_from_prev_SP ~ exp_group) %>% add_xy_position(x="exp_group") %>% pull(y.position)),
                      xmin=as.numeric(mixed_model_data %>% tukey_hsd(rel_length_from_prev_SP ~ exp_group) %>% add_xy_position(x="exp_group") %>% pull(xmin)),
                      xmax=as.numeric(mixed_model_data %>% tukey_hsd(rel_length_from_prev_SP ~ exp_group) %>% add_xy_position(x="exp_group") %>% pull(xmax)))

pwc$xmin <- as.numeric(sapply(pwc$group1, switch, "Perennial"=1, "Summer dry"=2, "More seasons dry"=3))
pwc$xmax <- as.numeric(sapply(pwc$group2, switch, "Perennial"=1, "Summer dry"=2, "More seasons dry"=3))

p <- ggplot(mixed_model_data, aes(x=factor(exp_group, levels=c("Perennial", "Summer dry", "More seasons dry")), y=rel_length_from_prev_SP)) + 
  geom_boxplot(aes(fill=exp_group), linewidth=0.8, median.linewidth=1.5, show.legend = FALSE) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$rel_length_from_prev_SP)+0.05, by=0.08, length.out=length(y.position))), aes(y=y.position, x=xmin, xend=xmax)) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$rel_length_from_prev_SP)+0.05, by=0.08, length.out=length(y.position))), aes(x=xmin, y=y.position, yend=y.position-0.01)) +
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$rel_length_from_prev_SP)+0.05, by=0.08, length.out=length(y.position))), aes(x=xmax, y=y.position, yend=y.position-0.01)) +
  geom_text(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$rel_length_from_prev_SP)+0.05, by=0.08, length.out=length(y.position))), aes(label=as.character(p.adj.signif), y=y.position, x=(xmin+xmax)/2,), position=position_nudge(y=0.02)) +
  scale_y_continuous(limits=c(0,1.2), breaks=seq(0,1,by=0.2)) +
  scale_x_discrete(labels=c("Perennial\n ", "Summer-dry\n ", "More-seasons-dry\n ")) +
  scale_fill_viridis(option="E", discrete=TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(margin=margin(0,0,2,0)),
        plot.subtitle = element_text(margin=margin(0,0,5,0)),
        axis.title.y = element_text(size=13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(size=13, color="black")) +
  xlab("") +
  ylab("Average segment length") +
  labs(title="", 
       subtitle=paste("Mixed model ANOVA, ", "F(",res.anova$NumDF[1], ",", res.anova$DenDF[1], ")","=",res.anova$F_value[1], ", ", "p=", res.anova$p[1], sep=""))
fig_5_plotlist <- append(fig_5_plotlist, list(p))
names(fig_5_plotlist)[length(fig_5_plotlist)] <- "5B"
remove(pwc, res.anova)

# Roseplot SP-SP average segment lengths (When?)
p <- ggplot(SP_SP_dist_table %>% left_join(grouping_table[,c(1,4)], join_by(SITE==SITE)), aes(x=reorder(factor(SITE), rel_length_from_prev_SP), y=rel_length_from_prev_SP, alpha=forcats::fct_rev(factor(YEAR_SEASON)), fill=factor(exp_group))) + 
  geom_col(data=SP_SP_dist_table %>% left_join(grouping_table[,c(1,4)], join_by(SITE==SITE)) %>% filter(exp_group=="Perennial"), color="darkgrey") +
  geom_col(data=SP_SP_dist_table %>% left_join(grouping_table[,c(1,4)], join_by(SITE==SITE)) %>% filter(exp_group=="Summer dry"), color="darkgrey") +
  geom_col(data=SP_SP_dist_table %>% left_join(grouping_table[,c(1,4)], join_by(SITE==SITE)) %>% filter(exp_group=="More seasons dry"), color="darkgrey") +
  scale_alpha_discrete(range=c(0.9,0.3), labels=c("2022-2023", "2021-2022", "2019-2021"),
                       guide=guide_legend(override.aes = list(fill="black", alpha=c(0.9,0.4,0.2)))) +
  scale_fill_manual(values=c(viridis(3, option="E")), labels=c("Perennial", "Summer-dry", "More-seasons-dry")) +
  coord_polar() +
  theme_bw() +
  theme(axis.title.y = element_text(size=13),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12)) +
  labs(title="", x="", y="Average segment length", fill="Seasonal occurrence\n of drying", alpha="Spring-to-spring interval       ")
fig_5_plotlist <- append(fig_5_plotlist, list(p))
names(fig_5_plotlist)[length(fig_5_plotlist)] <- "5A"

fig_5_plotlist <- fig_5_plotlist[order(names(fig_5_plotlist))]
ggarrange(plotlist=fig_5_plotlist, nrow=2, ncol=2, labels=c("A", "B", "C", "D"), font.label = list(size=25))
#ggsave("REV_FINAL_Fig_5.png", height = 8, width=16, bg="white")

# Compare distances from SP-SP interval centroid (How long?)
res.anova <- as.data.frame(anova(lmer(dist_from_int_cent ~ dry_group + YEAR_SEASON + (1|SITE), mixed_model_data)))
names(res.anova) <- gsub(" ", "_", names(res.anova))
res.anova <- res.anova %>% rename(p="Pr(>F)") %>% mutate(DenDF=round(DenDF,2), F_value=round(F_value, 2), p=ifelse(p<0.001, paste("<0.001"), round(p,3)))

pwc <- as.data.frame(contrast(emmeans(lmer(dist_from_int_cent ~ dry_group + YEAR_SEASON + (1|SITE), mixed_model_data), ~ dry_group), method="pairwise", adjust="tukey"))
pwc <- pwc %>% mutate(group1=str_before_first(contrast," -"), 
                      group2=str_after_first(contrast, " - "),
                      p.adj=p.value,
                      p.adj.signif=symnum(p.adj, corr=FALSE, na=FALSE, cutpoints=c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols=c("****", "***", "**", "*", "ns")),
                      y.position=as.numeric(mixed_model_data %>% tukey_hsd(dist_from_int_cent ~ dry_group) %>% add_xy_position(x="dry_group") %>% pull(y.position)),
                      xmin=as.numeric(mixed_model_data %>% tukey_hsd(dist_from_int_cent ~ dry_group) %>% add_xy_position(x="dry_group") %>% pull(xmin)),
                      xmax=as.numeric(mixed_model_data %>% tukey_hsd(dist_from_int_cent ~ dry_group) %>% add_xy_position(x="dry_group") %>% pull(xmax)))

p <- ggplot(mixed_model_data, aes(x=factor(dry_group, levels=c("DD365", "DD365+", "DD365++")), y=dist_from_int_cent)) + 
  geom_boxplot(aes(fill=dry_group), linewidth=0.8, median.linewidth=1.5, show.legend = FALSE) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$dist_from_int_cent)+0.05, by=0.08, length.out=length(y.position))), aes(y=y.position, x=xmin, xend=xmax)) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$dist_from_int_cent)+0.05, by=0.08, length.out=length(y.position))), aes(x=xmin, y=y.position, yend=y.position-0.01)) +
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$dist_from_int_cent)+0.05, by=0.08, length.out=length(y.position))), aes(x=xmax, y=y.position, yend=y.position-0.01)) +
  geom_text(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$dist_from_int_cent)+0.05, by=0.08, length.out=length(y.position))), aes(label=as.character(p.adj.signif), y=y.position, x=(xmin+xmax)/2,), position=position_nudge(y=0.02)) +
  scale_y_continuous(limits=c(0,0.7), breaks=seq(0,0.6, by=0.2)) +
  scale_x_discrete(labels=c("DD365\nQuasi-perennial", "DD365+\nShort-term drying", "DD365++\nLong-term drying")) +
  scale_fill_viridis(option="D", discrete=TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(margin=margin(0,0,2,0)),
        plot.subtitle = element_text(margin=margin(0,0,5,0)),
        axis.title.y = element_text(size=13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(size=13, color="black")) +
  xlab("") +
  ylab("Distance from interval centroid") +
  labs(title="", 
       subtitle=paste("Mixed model ANOVA, ", "F(",res.anova$NumDF[1], ",", res.anova$DenDF[1], ")","=",res.anova$F_value[1], ", ", "p=", res.anova$p[1], sep=""))
fig_6_plotlist <- append(fig_6_plotlist, list(p))
names(fig_6_plotlist)[length(fig_6_plotlist)] <- "6D"
remove(pwc, res.anova)

p <- ggplot(SP_SP_dist_table, aes(x=factor(dry_group, levels=c("DD365", "DD365+", "DD365++")), y=dist_from_int_cent, fill=dry_group, alpha=plyr::revalue(factor(CAMPAIGN_Nr), c("11" = "2019-2021", "15" = "2021-2022", "19"="2022-2023")))) + 
  scale_fill_viridis(option="D", direction=1, discrete=TRUE) +
  scale_x_discrete(labels=c("DD365\nQuasi-perennial", "DD365+\nShort-term drying", "DD365++\nLong-term drying")) +
  scale_y_continuous(limits=c(0,0.7), breaks=seq(0,0.6, by=0.2)) +
  scale_alpha_discrete(range=c(0.2,1), labels=c("2022-2023", "2021-2022", "2019-2021"),
                       guide=guide_legend(override.aes = list(fill="black", alpha=c(0.9,0.4,0.2)))) +
  geom_boxplot(linewidth=0.8, median.linewidth=1.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        axis.title.y = element_text(size=13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(size=13, color="black")) +
  xlab("") +
  ylab("Distance from interval centroid") +
  labs(title="", subtitle="", fill="Duration of drying", alpha="Spring-to-spring interval")
fig_6_plotlist <- append(fig_6_plotlist, list(p))
names(fig_6_plotlist)[length(fig_6_plotlist)] <- "6C"

# Compare distances from SP-SP interval centroid (When?)
res.anova <- as.data.frame(anova(lmer(dist_from_int_cent ~ exp_group + YEAR_SEASON + (1|SITE), mixed_model_data)))
names(res.anova) <- gsub(" ", "_", names(res.anova))
res.anova <- res.anova %>% rename(p="Pr(>F)") %>% mutate(DenDF=round(DenDF,2), F_value=round(F_value, 2), p=ifelse(p<0.001, paste("<0.001"), round(p,3)))

pwc <- as.data.frame(contrast(emmeans(lmer(dist_from_int_cent ~ exp_group + YEAR_SEASON + (1|SITE), mixed_model_data), ~ exp_group), method="pairwise", adjust="tukey"))
pwc <- pwc %>% mutate(group1=str_before_first(contrast," -"), 
                      group2=str_after_first(contrast, " - "),
                      p.adj=p.value,
                      p.adj.signif=symnum(p.adj, corr=FALSE, na=FALSE, cutpoints=c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols=c("****", "***", "**", "*", "ns")),
                      y.position=as.numeric(mixed_model_data %>% tukey_hsd(dist_from_int_cent ~ exp_group) %>% add_xy_position(x="exp_group") %>% pull(y.position)),
                      xmin=as.numeric(mixed_model_data %>% tukey_hsd(dist_from_int_cent ~ exp_group) %>% add_xy_position(x="exp_group") %>% pull(xmin)),
                      xmax=as.numeric(mixed_model_data %>% tukey_hsd(dist_from_int_cent ~ exp_group) %>% add_xy_position(x="exp_group") %>% pull(xmax)))

pwc$xmin <- as.numeric(sapply(pwc$group1, switch, "Perennial"=1, "Summer dry"=2, "More seasons dry"=3))
pwc$xmax <- as.numeric(sapply(pwc$group2, switch, "Perennial"=1, "Summer dry"=2, "More seasons dry"=3))

p <- ggplot(mixed_model_data, aes(x=factor(exp_group, levels=c("Perennial", "Summer dry", "More seasons dry")), y=dist_from_int_cent)) + 
  geom_boxplot(aes(fill=exp_group), linewidth=0.8, median.linewidth=1.5, show.legend = FALSE) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$dist_from_int_cent)+0.05, by=0.08, length.out=length(y.position))), aes(y=y.position, x=xmin, xend=xmax)) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$dist_from_int_cent)+0.05, by=0.08, length.out=length(y.position))), aes(x=xmin, y=y.position, yend=y.position-0.01)) +
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$dist_from_int_cent)+0.05, by=0.08, length.out=length(y.position))), aes(x=xmax, y=y.position, yend=y.position-0.01)) +
  geom_text(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(mixed_model_data$dist_from_int_cent)+0.05, by=0.08, length.out=length(y.position))), aes(label=as.character(p.adj.signif), y=y.position, x=(xmin+xmax)/2,), position=position_nudge(y=0.02)) +
  scale_y_continuous(limits=c(0,0.7), breaks=seq(0,0.6, by=0.2)) +
  scale_x_discrete(labels=c("Perennial", "Summer-dry", "More-seasons-dry")) +
  scale_fill_viridis(option="E", discrete=TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(margin=margin(0,0,2,0)),
        plot.subtitle = element_text(margin=margin(0,0,5,0)),
        axis.title.y = element_text(size=13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(size=13, color="black")) +
  xlab("") +
  ylab("Distance from interval centroid") +
  labs(title="", 
       subtitle=paste("Mixed model ANOVA, ", "F(",res.anova$NumDF[1], ",", res.anova$DenDF[1], ")","=",res.anova$F_value[1], ", ", "p=", res.anova$p[1], sep=""))
fig_6_plotlist <- append(fig_6_plotlist, list(p))
names(fig_6_plotlist)[length(fig_6_plotlist)] <- "6B"
remove(pwc, res.anova)

p <- ggplot(na.omit(SP_SP_dist_table %>% left_join(grouping_table[,c(1,4)], join_by(SITE==SITE)) %>% mutate(exp_group=factor(exp_group, levels=c("Perennial", "Summer dry", "More seasons dry")))), aes(x=factor(exp_group, levels=c("Perennial", "Summer dry", "More sesons dry")), y=dist_from_int_cent, fill=exp_group, alpha=plyr::revalue(factor(CAMPAIGN_Nr), c("11" = "2019-2021", "15" = "2021-2022", "19"="2022-2023")))) + 
  scale_fill_manual(values=c(viridis(3, option="E", direction=1))) +
  scale_x_discrete(labels=c("Perennial", "Summer-dry", "More-seasons-dry")) +
  scale_y_continuous(limits=c(0,0.7), breaks=seq(0,0.6, by=0.2)) +
  scale_alpha_discrete(range=c(0.2,1), labels=c("2022-2023", "2021-2022", "2019-2021"),
                       guide=guide_legend(override.aes = list(fill="black", alpha=c(0.9,0.4,0.2)))) +
  geom_boxplot(linewidth=0.8, median.linewidth=1.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        axis.title.y = element_text(size=13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(size=13, color="black")) +
  xlab("") +
  ylab("Distance from interval centroid") +
  labs(title="", subtitle="", fill="Seasonal occurrence \nof drying", alpha="Spring-to-spring interval")
fig_6_plotlist <- append(fig_6_plotlist, list(p))
names(fig_6_plotlist)[length(fig_6_plotlist)] <- "6A"

fig_6_plotlist <- fig_6_plotlist[order(names(fig_6_plotlist))]
ggarrange(plotlist=fig_6_plotlist, nrow=2, ncol=2, labels=c("A", "B", "C", "D"), font.label = list(size=25))
#ggsave("REV_FINAL_Fig_6.png", height = 8, width=16, bg="white")

# Plot each solo trajectories per group

sample_solo_plot_table <- data.frame()
for (i in c("DD365", "DD365+", "DD365++")){
  selected_site <- unique(SP_SP_dist_table %>% filter(dry_group==i) %>% pull(SITE))
  selected_campaign <- c(3:19)
  
  if(any(is.na(match(selected_site, unique(data$SITE))))){
    stop("One or more selected site number cannot be found in the dataset")}
  
  if(any(as.numeric(data %>% filter(SITE %in% selected_site, CAMPAIGN_Nr %in% selected_campaign) %>% group_by(SITE) %>% summarise(n_CAMP=length(CAMPAIGN_Nr)) %>% pull(n_CAMP))<2)){
    stop("One or more selected site have not enough observations in the dataset")}
  
  # Create filtered trajectory(s) of selected sites and campaigns to get eigenvalues and reduced number of dimensions
  # available only for sites with at least 2 samples
  
  FILTERED_DATASET_FILTERED_TRAJECTORIES <- list()
  sites <- data %>% filter(SITE %in% selected_site, CAMPAIGN_Nr %in% selected_campaign) %>% pull(SITE)
  surveys <- data %>% filter(SITE %in% selected_site, CAMPAIGN_Nr %in% selected_campaign) %>% pull(CAMPAIGN_Nr)
  dev.new()
  tr_def <- defineTrajectories(vegan::vegdist(as.data.frame(data %>% filter(SITE %in% selected_site, CAMPAIGN_Nr %in% selected_campaign) %>% select((samp_ID_ncol+1):ncol(data)))[,as.numeric(which(colSums(as.data.frame(data %>% filter(SITE %in% selected_site, CAMPAIGN_Nr %in% selected_campaign) %>% select((samp_ID_ncol+1):ncol(data))))!=0))], method=dist_metric), sites, surveys)
  tr <- trajectoryPCoA(tr_def)
  dev.off()
  FILTERED_DATASET_FILTERED_TRAJECTORIES <- append(FILTERED_DATASET_FILTERED_TRAJECTORIES, list(tr$eig))
  FILTERED_DATASET_FILTERED_TRAJECTORIES <- append(FILTERED_DATASET_FILTERED_TRAJECTORIES, list(data.frame("SITE"=sites, "CAMPAIGN_Nr"=surveys, tr$points)))
  names(FILTERED_DATASET_FILTERED_TRAJECTORIES) <- c("tr_eig", "tr_points")
  remove(tr)
  
  center_plot <- 0 # Create centered 2D plot? (0=NO, 1=YES)
  
  # Plot table
  plot_table <- rename(as.data.frame(FILTERED_DATASET_FILTERED_TRAJECTORIES$tr_points[,c(1:4)]), X=X1, Y=X2)
  
  if(center_plot==1){
    plot_table <- plot_table %>% left_join(as.data.frame(plot_table %>% group_by(SITE) %>% summarise(ini_x=head(X,1), ini_y=head(Y,1))), join_by(SITE==SITE)) %>% mutate(X=X-ini_x, Y=Y-ini_y) %>% select(SITE, CAMPAIGN_Nr, X,Y)}
  
  plot_table <- plot_table %>% left_join(grouping_table[,c(1,3)], join_by(SITE==SITE))
  
  if(any(unique(plot_table$SITE)==15)){
    sample_solo_plot_table <- rbind(sample_solo_plot_table, plot_table[which(plot_table$SITE==15),c(1:4)])
  }
  
  if(any(unique(plot_table$SITE)==45)){
    sample_solo_plot_table <- rbind(sample_solo_plot_table, plot_table[which(plot_table$SITE==45),c(1:4)])
  }
  
  if(any(unique(plot_table$SITE)==22)){
    sample_solo_plot_table <- rbind(sample_solo_plot_table, plot_table[which(plot_table$SITE==22),c(1:4)])
  }
  
  p <- ggplot(plot_table, aes(x=X, y=Y)) + 
    geom_path(arrow=arrow(type="closed", length = unit(0.8, "picas")), linewidth=0.8) +
    theme_bw() +
    geom_vline(xintercept = 0, lty=3) +
    geom_hline(yintercept = 0, lty=3) +
    scale_x_continuous(limits=c(-1,1)) +
    scale_y_continuous(limits=c(-1,1)) +
    xlab("PCoA 1") +
    ylab("PCoA 2") +
    labs(title=paste("2D trajectory plots of ", i, " sites", sep="")) +
    facet_wrap(~SITE,ncol=3, nrow=ceiling(length(unique(plot_table$SITE))/3))
  plot(p)
  #ggsave(paste("REV_FINAL_Fig_S", match(i, c("DD365", "DD365+", "DD365++"))+3, ".png", sep=""), height = ceiling(length(unique(plot_table$SITE))/3)*10/4, width=8, bg="white")
}

## TEST CONVERGENCE OF INTERPOLATED TRAJECTORIES  ##

interpolation_summary <- data.frame()
convergence_results <- data.frame()

for (i in c(unique(grouping_table$dry_group))){
  selected_site <- grouping_table %>% filter(dry_group==i) %>% pull(SITE)
  selected_campaign <- c(3:19)
  
  if(any(is.na(match(selected_site, unique(data$SITE))))){
    stop("One or more selected site number cannot be found in the dataset")}
  
  if(any(as.numeric(data %>% filter(SITE %in% selected_site, CAMPAIGN_Nr %in% selected_campaign) %>% group_by(SITE) %>% summarise(n_CAMP=length(CAMPAIGN_Nr)) %>% pull(n_CAMP))<2)){
    stop("One or more selected site have not enough observations in the dataset")}
  
  # Create filtered trajectory(s) of selected sites and campaigns to get eigenvalues and reduced number of dimensions
  # available only for sites with at least 2 samples
  
  FILTERED_DATASET_FILTERED_TRAJECTORIES <- list()
  sites <- data %>% filter(SITE %in% selected_site, CAMPAIGN_Nr %in% selected_campaign) %>% pull(SITE)
  surveys <- data %>% filter(SITE %in% selected_site, CAMPAIGN_Nr %in% selected_campaign) %>% pull(CAMPAIGN_Nr)
  dist_metric = "bray"
  tr_def <- defineTrajectories(vegan::vegdist(as.data.frame(data %>% filter(SITE %in% selected_site, CAMPAIGN_Nr %in% selected_campaign) %>% select((samp_ID_ncol+1):ncol(data)))[,as.numeric(which(colSums(as.data.frame(data %>% filter(SITE %in% selected_site, CAMPAIGN_Nr %in% selected_campaign) %>% select((samp_ID_ncol+1):ncol(data))))!=0))], method=dist_metric), sites, surveys)
  tr_int_def <- interpolateTrajectories(tr_def, times=as.numeric(names(which(table(data %>% filter(SITE %in% c(grouping_table$SITE[which(grouping_table$dry_group==i)]), between(CAMPAIGN_Nr, 3,19)) %>% pull(CAMPAIGN_Nr))>=floor(mean(table(data %>% filter(SITE %in% c(grouping_table$SITE[which(grouping_table$dry_group==i)]), between(CAMPAIGN_Nr, 3,19)) %>% pull(CAMPAIGN_Nr))))))))
  
  convergence_results <- rbind(convergence_results, data.frame("group"=i, "Mann_Kendall_tau"=as.numeric(round(trajectoryConvergence(tr_int_def, type = "multiple", add=TRUE)$tau,3)),"sign"=ifelse(trajectoryConvergence(tr_int_def, type = "multiple", add=TRUE)$p.value<0.05,"p<0.05",paste("p=", round(trajectoryConvergence(tr_int_def, type = "multiple", add=TRUE)$p.value,2), sep=""))))
  interpolation_summary <- rbind(interpolation_summary, as.data.frame(tr_def$metadata %>% group_by(sites) %>% summarise(n_obs=length(times), .groups="drop") %>% arrange(as.numeric(sites)) %>% mutate(Grouping="How long?", Group=i, SITE=sites, N_obs_orig=n_obs) %>% select(Grouping, Group, SITE, N_obs_orig)) %>% left_join(as.data.frame(tr_int_def$metadata %>% group_by(sites) %>% summarise(N_obs_int=length(times), .groups="drop")), join_by(SITE==sites)))
}

for (i in c(unique(grouping_table$exp_group))){
  selected_site <- grouping_table %>% filter(exp_group==i) %>% pull(SITE)
  selected_campaign <- c(3:19)
  
  if(any(is.na(match(selected_site, unique(data$SITE))))){
    stop("One or more selected site number cannot be found in the dataset")}
  
  if(any(as.numeric(data %>% filter(SITE %in% selected_site, CAMPAIGN_Nr %in% selected_campaign) %>% group_by(SITE) %>% summarise(n_CAMP=length(CAMPAIGN_Nr)) %>% pull(n_CAMP))<2)){
    stop("One or more selected site have not enough observations in the dataset")}
  
  # Create filtered trajectory(s) of selected sites and campaigns to get eigenvalues and reduced number of dimensions
  # available only for sites with at least 2 samples
  
  FILTERED_DATASET_FILTERED_TRAJECTORIES <- list()
  sites <- data %>% filter(SITE %in% selected_site, CAMPAIGN_Nr %in% selected_campaign) %>% pull(SITE)
  surveys <- data %>% filter(SITE %in% selected_site, CAMPAIGN_Nr %in% selected_campaign) %>% pull(CAMPAIGN_Nr)
  dist_metric = "bray"
  tr_def <- defineTrajectories(vegan::vegdist(as.data.frame(data %>% filter(SITE %in% selected_site, CAMPAIGN_Nr %in% selected_campaign) %>% select((samp_ID_ncol+1):ncol(data)))[,as.numeric(which(colSums(as.data.frame(data %>% filter(SITE %in% selected_site, CAMPAIGN_Nr %in% selected_campaign) %>% select((samp_ID_ncol+1):ncol(data))))!=0))], method=dist_metric), sites, surveys)
  tr_int_def <- interpolateTrajectories(tr_def, times=as.numeric(names(which(table(data %>% filter(SITE %in% c(grouping_table$SITE[which(grouping_table$exp_group==i)]), between(CAMPAIGN_Nr, 3,19)) %>% pull(CAMPAIGN_Nr))>=floor(mean(table(data %>% filter(SITE %in% c(grouping_table$SITE[which(grouping_table$exp_group==i)]), between(CAMPAIGN_Nr, 3,19)) %>% pull(CAMPAIGN_Nr))))))))
  
  convergence_results <- rbind(convergence_results, data.frame("group"=i, "Mann_Kendall_tau"=as.numeric(round(trajectoryConvergence(tr_int_def, type = "multiple", add=TRUE)$tau,3)),"sign"=ifelse(trajectoryConvergence(tr_int_def, type = "multiple", add=TRUE)$p.value<0.05,"p<0.05",paste("p=", round(trajectoryConvergence(tr_int_def, type = "multiple", add=TRUE)$p.value,2), sep=""))))
  interpolation_summary <- rbind(interpolation_summary, as.data.frame(tr_def$metadata %>% group_by(sites) %>% summarise(n_obs=length(times), .groups="drop") %>% arrange(as.numeric(sites)) %>% mutate(Grouping="When?", Group=i, SITE=sites, N_obs_orig=n_obs) %>% select(Grouping, Group, SITE, N_obs_orig)) %>% left_join(as.data.frame(tr_int_def$metadata %>% group_by(sites) %>% summarise(N_obs_int=length(times), .groups="drop")), join_by(SITE==sites)))
}

interpolation_summary <- interpolation_summary %>% mutate(diff=N_obs_int-N_obs_orig)

## TEST "REDUCED" TRAJECTORIES TO VALIDATE AVERAGE SEGMENT LENGTHS ##

reduced_sites <- c(unique(SP_SP_dist_table$SITE)[which(!unique(SP_SP_dist_table$SITE) %in% c(40,49,13))])
reference_site <- as.numeric(38)

REDUCED_DATASET_FULL_TRAJECTORIES <- list()
data_reduced <- data %>% filter(SITE %in% c(reduced_sites), CAMPAIGN_Nr %in% c(as.numeric(data %>% filter(SITE==reference_site) %>% pull(CAMPAIGN_Nr))))

data_reduced$SITE_ORDER <- rep(as.numeric(data_reduced %>% filter(SITE==reference_site) %>% pull(SITE_ORDER)), length(unique(data_reduced$SITE)))

data_reduced <- data_reduced[,-c(as.numeric(which(colSums(data_reduced[, (samp_ID_ncol+1):ncol(data_reduced)])==0))+samp_ID_ncol)]
sites <- data_reduced$SITE
surveys <- data_reduced$CAMPAIGN_Nr
dev.new()
tr_def <- defineTrajectories(vegan::vegdist(data_reduced[,(samp_ID_ncol+1):ncol(data_reduced)], method=dist_metric), sites, surveys)
tr <- trajectoryPCoA(tr_def)
dev.off()

REDUCED_DATASET_FULL_TRAJECTORIES <- append(REDUCED_DATASET_FULL_TRAJECTORIES, list(tr$eig))
REDUCED_DATASET_FULL_TRAJECTORIES <- append(REDUCED_DATASET_FULL_TRAJECTORIES, list(data.frame("SITE"=sites, "CAMPAIGN_Nr"=surveys, tr$points)))
REDUCED_DATASET_FULL_TRAJECTORIES <- append(REDUCED_DATASET_FULL_TRAJECTORIES, list(data.frame("SITE"=unique(sites),trajectoryLengths(tr_def)) %>% as.data.frame(row.names = 1:nrow(.))))
remove(tr)
segment_dist_IF <- data.frame(segmentDistances(tr_def)$Dinifin)
rownames(segment_dist_IF) <- str_replace_all(rownames(segment_dist_IF), "\\[", "_")
rownames(segment_dist_IF) <- str_replace_all(rownames(segment_dist_IF), "-", "_")
rownames(segment_dist_IF) <- paste("I", str_sub(rownames(segment_dist_IF),1,-2), sep="_")
colnames(segment_dist_IF) <- paste("F", str_sub(str_replace_all(colnames(segment_dist_IF), "\\.", "_"), 2, -2), sep="_")
REDUCED_DATASET_FULL_TRAJECTORIES <- append(REDUCED_DATASET_FULL_TRAJECTORIES, list(data.frame(segment_dist_IF)))
remove(segment_dist_IF)
segment_dist_FF <- data.frame(as.matrix(segmentDistances(tr_def)$Dfin))
rownames(segment_dist_FF) <- str_replace_all(rownames(segment_dist_FF), "\\[", "_")
rownames(segment_dist_FF) <- str_replace_all(rownames(segment_dist_FF), "-", "_")
rownames(segment_dist_FF) <- paste("F", str_sub(rownames(segment_dist_FF),1,-2), sep="_")
colnames(segment_dist_FF) <- paste("F", str_sub(str_replace_all(colnames(segment_dist_FF), "\\.", "_"), 2, -2), sep="_")
REDUCED_DATASET_FULL_TRAJECTORIES <- append(REDUCED_DATASET_FULL_TRAJECTORIES, list(data.frame(segment_dist_FF)))
remove(segment_dist_FF)
names(REDUCED_DATASET_FULL_TRAJECTORIES) <- c("tr_eig", "tr_points", "tr_lengths", "segment_dist_IF", "segment_dist_FF")

SP_SP_dist_table <- data.frame(data_reduced %>% filter(SITE %in% unique(grouping_table$SITE)) %>% select(SITE, YEAR_SEASON, CAMPAIGN_Nr, SEASON) %>% filter(SEASON=="SP"))
SP_SP_dist_table <- SP_SP_dist_table %>% 
  mutate(from=c(NA,as.numeric(data_reduced %>% filter(SITE %in% unique(grouping_table$SITE)) %>% select(SITE, CAMPAIGN_Nr, SEASON, SITE_ORDER) %>% filter(SEASON=="SP") %>% pull(CAMPAIGN_Nr))[1:(length(as.numeric(data_reduced %>% filter(SITE %in% unique(grouping_table$SITE)) %>% select(SITE, CAMPAIGN_Nr, SEASON, SITE_ORDER) %>% filter(SEASON=="SP") %>% pull(CAMPAIGN_Nr)))-1)]), to=CAMPAIGN_Nr) %>% filter(to!=3)
SP_SP_dist_table <- SP_SP_dist_table %>% left_join(data_reduced[,c(1,5,6)], join_by(SITE==SITE, from==CAMPAIGN_Nr))

SP_SP_dist_list <- NULL
for(i in c(1:nrow(SP_SP_dist_table))){
  SP_SP_dist_list <- append(SP_SP_dist_list, as.numeric(REDUCED_DATASET_FULL_TRAJECTORIES$segment_dist_IF[grep(paste("I_", SP_SP_dist_table[i,]$SITE, "_",SP_SP_dist_table[i,]$from,"_",  sep=""), rownames(REDUCED_DATASET_FULL_TRAJECTORIES$segment_dist_IF)),match(colnames(REDUCED_DATASET_FULL_TRAJECTORIES$segment_dist_IF)[grep(paste("F_", SP_SP_dist_table[i,]$SITE, "_", sep=""), colnames(REDUCED_DATASET_FULL_TRAJECTORIES$segment_dist_IF))][match(SP_SP_dist_table[i,]$to, str_after_last(colnames(REDUCED_DATASET_FULL_TRAJECTORIES$segment_dist_IF)[grep(paste("F_", SP_SP_dist_table[i,]$SITE, "_", sep=""), colnames(REDUCED_DATASET_FULL_TRAJECTORIES$segment_dist_IF))], "_"))],colnames(REDUCED_DATASET_FULL_TRAJECTORIES$segment_dist_IF))]))
}

SP_SP_dist_table$SP_SP_dist <- SP_SP_dist_list
SP_SP_dist_table <- SP_SP_dist_table %>% mutate(order_from=SITE_ORDER, order_to=c(as.numeric(SP_SP_dist_table %>% mutate(order_from=SITE_ORDER) %>% pull(order_from))[2:nrow(SP_SP_dist_table)]-1,0))

for (i in which(SP_SP_dist_table$order_to<SP_SP_dist_table$order_from)){
  SP_SP_dist_table$order_to[i] <- max(as.numeric(str_after_first(colnames(REDUCED_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))[which(!is.na(REDUCED_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE)))][-c(1, length(colnames(REDUCED_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))[which(!is.na(REDUCED_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE)))]))],"S")))
}

SP_SP_length_list <- NULL
SP_SP_rel_length_list <- NULL
for(i in c(1:nrow(SP_SP_dist_table))){
  SP_SP_length_list <- append(SP_SP_length_list, sum(as.numeric(as.data.frame(REDUCED_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))[,c(match(paste("S", SP_SP_dist_table[i,]$order_from, sep=""), colnames(REDUCED_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))):match(paste("S", SP_SP_dist_table[i,]$order_to, sep=""), colnames(REDUCED_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))))])))
  SP_SP_rel_length_list <- append(SP_SP_rel_length_list, sum(as.numeric(as.data.frame(REDUCED_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))[,c(match(paste("S", SP_SP_dist_table[i,]$order_from, sep=""), colnames(REDUCED_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))):match(paste("S", SP_SP_dist_table[i,]$order_to, sep=""), colnames(REDUCED_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))))]))/length(as.numeric(as.data.frame(REDUCED_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))[,c(match(paste("S", SP_SP_dist_table[i,]$order_from, sep=""), colnames(REDUCED_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))):match(paste("S", SP_SP_dist_table[i,]$order_to, sep=""), colnames(REDUCED_DATASET_FULL_TRAJECTORIES$tr_lengths %>% filter(SITE==SP_SP_dist_table[i,]$SITE))))]))) 
}
SP_SP_dist_table$SP_SP_length <- SP_SP_length_list
SP_SP_dist_table$SP_SP_rel_length <- SP_SP_rel_length_list

SP_SP_dist_table <- SP_SP_dist_table %>% left_join(grouping_table[,c(1,3)], join_by(SITE==SITE))

SP_SP_dist_table <- SP_SP_dist_table %>% reframe(SITE=SITE, 
                                                 YEAR_SEASON=YEAR_SEASON, 
                                                 SEASON=SEASON, 
                                                 CAMPAIGN_Nr=CAMPAIGN_Nr, 
                                                 SITE_ORDER=SITE_ORDER, 
                                                 dist_from_prev_SP=SP_SP_dist, 
                                                 length_from_prev_SP=SP_SP_length, 
                                                 rel_length_from_prev_SP=SP_SP_rel_length,
                                                 dry_group=dry_group)

SP_SP_dist_table <- SP_SP_dist_table %>% left_join(grouping_table[, c(1,4)], join_by(SITE==SITE))

avg_length_full_length_comp <- avg_length_full_length_comp %>% filter(!SITE %in% c(40,49,13))
avg_length_full_length_comp <- cbind(avg_length_full_length_comp, "full_length"=SP_SP_dist_table$length_from_prev_SP)
avg_length_full_length_comp <- as.data.frame(avg_length_full_length_comp %>% group_by(SITE) %>% summarise(avg_length=mean(avg_length), full_length=sum(full_length), .groups="drop") %>% left_join(grouping_table[,c(1,3,4)], join_by(SITE==SITE)))

avg_length_full_length_comp_cor <- data.frame()
for (i in unique(avg_length_full_length_comp$dry_group)){
  avg_length_full_length_comp_cor <- rbind(avg_length_full_length_comp_cor, data.frame("group"=i, "r"=round(as.data.frame(avg_length_full_length_comp %>% filter(dry_group==i) %>% cor_test(avg_length, full_length, method="spearman"))$cor, 2), "p"=ifelse(as.data.frame(avg_length_full_length_comp %>% filter(dry_group==i) %>% cor_test(avg_length, full_length, method="spearman"))$p<0.05, "p<0.05", paste("p=",round(as.data.frame(avg_length_full_length_comp %>% filter(dry_group==i) %>% cor_test(avg_length, full_length, method="spearman"))$p, 2), sep=""))))
}
for (i in unique(avg_length_full_length_comp$exp_group)){
  avg_length_full_length_comp_cor <- rbind(avg_length_full_length_comp_cor, data.frame("group"=i, "r"=round(as.data.frame(avg_length_full_length_comp %>% filter(exp_group==i) %>% cor_test(avg_length, full_length, method="spearman"))$cor, 2), "p"=ifelse(as.data.frame(avg_length_full_length_comp %>% filter(exp_group==i) %>% cor_test(avg_length, full_length, method="spearman"))$p<0.05, "p<0.05", paste("p=",round(as.data.frame(avg_length_full_length_comp %>% filter(exp_group==i) %>% cor_test(avg_length, full_length, method="spearman"))$p, 2), sep=""))))
}
avg_length_full_length_comp_cor <- rbind(avg_length_full_length_comp_cor, data.frame("group"="Overall", "r"=round(as.data.frame(avg_length_full_length_comp %>% cor_test(avg_length, full_length, method="spearman"))$cor, 2), "p"=ifelse(as.data.frame(avg_length_full_length_comp %>% cor_test(avg_length, full_length, method="spearman"))$p<0.05, "p<0.05", paste("p=",round(as.data.frame(avg_length_full_length_comp %>% cor_test(avg_length, full_length, method="spearman"))$p, 2), sep=""))))

ggplot(avg_length_full_length_comp, aes(x=full_length, y=avg_length, color=factor(dry_group, levels=c("DD365", "DD365+","DD365++", "Overall")))) + 
  geom_point(size=3, alpha=0.6) +
  geom_smooth(aes(color=dry_group), formula=y~x, method="lm", se=FALSE) +
  geom_smooth(method="lm", se=FALSE, aes(color="Overall"), formula=y~x) +
  scale_color_manual(values=c(viridis(3), "darkgrey"), guide="legend", breaks=c("DD365", "DD365+", "DD365++", "Overall"), 
                     labels=c(paste(avg_length_full_length_comp_cor[1,1], " (r=",avg_length_full_length_comp_cor[1,2], ", ",avg_length_full_length_comp_cor[1,3], ")",  sep=""),
                              paste(avg_length_full_length_comp_cor[2,1], " (r=",avg_length_full_length_comp_cor[2,2], ", ",avg_length_full_length_comp_cor[2,3], ")",  sep=""), 
                              paste(avg_length_full_length_comp_cor[3,1], " (r=",avg_length_full_length_comp_cor[3,2], ", ",avg_length_full_length_comp_cor[3,3], ")",  sep=""), 
                              paste(avg_length_full_length_comp_cor[7,1], " (r=",avg_length_full_length_comp_cor[7,2], ", ",avg_length_full_length_comp_cor[7,3], ")",  sep=""))) +
  theme_bw() +
  xlab("Trajectory length (arbitrarily reduced control dataset)") +
  ylab("Average trajectory segment length (original dataset)") +
  labs(title="Comparison of average segment length and full length", color="Duration of drying")
#ggsave("REV_FINAL_Fig_S2b.png", height = 5, width=10, bg="white")

ggplot(avg_length_full_length_comp, aes(x=full_length, y=avg_length, color=factor(exp_group, levels=c("Perennial", "Summer dry", "More seasons dry", "Overall")))) + 
  geom_point(size=3, alpha=0.6) +
  geom_smooth(aes(color=exp_group), formula=y~x, method="lm", se=FALSE) +
  geom_smooth(method="lm", se=FALSE, aes(color="Overall"), formula=y~x) +
  scale_color_manual(values=c(viridis(3, option="E"), "darkgrey"), guide="legend", breaks=c("Perennial", "Summer-dry", "More-seasons-dry", "Overall"), 
                     labels=c(paste(avg_length_full_length_comp_cor[4,1], " (r=",avg_length_full_length_comp_cor[4,2], ", ",avg_length_full_length_comp_cor[4,3], ")",  sep=""),
                              paste(avg_length_full_length_comp_cor[5,1], " (r=",avg_length_full_length_comp_cor[5,2], ", ",avg_length_full_length_comp_cor[5,3], ")",  sep=""), 
                              paste(avg_length_full_length_comp_cor[6,1], " (r=",avg_length_full_length_comp_cor[6,2], ", ",avg_length_full_length_comp_cor[6,3], ")",  sep=""), 
                              paste(avg_length_full_length_comp_cor[7,1], " (r=",avg_length_full_length_comp_cor[7,2], ", ",avg_length_full_length_comp_cor[7,3], ")",  sep=""))) +
  theme_bw() +
  xlab("Trajectory length (arbitrarily reduced control dataset)") +
  ylab("Average trajectory segment length (original dataset)") +
  labs(title="Comparison of average segment length and full length", color="Seasonal occurrence of drying")
#ggsave("REV_FINAL_Fig_S2a.png", height = 5, width=10, bg="white")

reduced_dataset_plotlist <- list()

# Compare reduced dataset full traj lengths in How long? groups
res.anova <- as.data.frame(SP_SP_dist_table %>% group_by(SITE) %>% summarise(tr_length=sum(length_from_prev_SP), dry_group=head(dry_group,1), .groups="drop")) %>% anova_test(tr_length ~ dry_group)
pwc <- tukey_hsd(as.data.frame(SP_SP_dist_table %>% group_by(SITE) %>% summarise(tr_length=sum(length_from_prev_SP), dry_group=head(dry_group,1), .groups="drop")), tr_length ~ dry_group)
pwc <- pwc %>% add_xy_position(x = "dry_group")

pwc$xmin <- as.numeric(sapply(pwc$group1, switch, "DD365"=1, "DD365+"=2, "DD365++"=3))
pwc$xmax <- as.numeric(sapply(pwc$group2, switch, "DD365"=1, "DD365+"=2, "DD365++"=3))

p <- ggplot(as.data.frame(SP_SP_dist_table %>% group_by(SITE) %>% summarise(tr_length=sum(length_from_prev_SP), dry_group=head(dry_group,1), .groups="drop")), aes(x=factor(dry_group, levels=c("DD365", "DD365+", "DD365++")), y=tr_length)) + 
  geom_boxplot(aes(fill=dry_group), show.legend = FALSE) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(as.data.frame(SP_SP_dist_table %>% group_by(SITE) %>% summarise(tr_length=sum(length_from_prev_SP), dry_group=head(dry_group,1), .groups="drop"))$tr_length)+0.3, by=0.4, length.out=length(y.position))), aes(y=y.position, x=xmin, xend=xmax)) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(as.data.frame(SP_SP_dist_table %>% group_by(SITE) %>% summarise(tr_length=sum(length_from_prev_SP), dry_group=head(dry_group,1), .groups="drop"))$tr_length)+0.3, by=0.4, length.out=length(y.position))), aes(x=xmin, y=y.position, yend=y.position-0.1)) +
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(as.data.frame(SP_SP_dist_table %>% group_by(SITE) %>% summarise(tr_length=sum(length_from_prev_SP), dry_group=head(dry_group,1), .groups="drop"))$tr_length)+0.3, by=0.4, length.out=length(y.position))), aes(x=xmax, y=y.position, yend=y.position-0.1)) +
  geom_text(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(as.data.frame(SP_SP_dist_table %>% group_by(SITE) %>% summarise(tr_length=sum(length_from_prev_SP), dry_group=head(dry_group,1), .groups="drop"))$tr_length)+0.3, by=0.4, length.out=length(y.position))), aes(label=p.adj.signif, y=y.position, x=(xmin+xmax)/2,), position=position_nudge(y=0.1)) +
  scale_y_continuous(limits=c(0,NA)) +
  scale_x_discrete(labels=c("DD365\nQuasi-perennial", "DD365+\nShort-term drying", "DD365++\nLong-term drying")) +
  scale_fill_viridis(option="D", discrete=TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(margin=margin(0,0,2,0)),
        plot.subtitle = element_text(margin=margin(0,0,0,0)),
        axis.title.y = element_text(size=13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(size=13, color="black")) +
  xlab("") +
  ylab("Trajectory length") +
  labs(title="", 
       subtitle=get_test_label(res.anova, detailed = TRUE))
reduced_dataset_plotlist <- append(reduced_dataset_plotlist, list(p))
names(reduced_dataset_plotlist)[length(reduced_dataset_plotlist)] <- "B"
remove(pwc, res.anova)

# Compare reduced dataset full traj lengths in When? groups
res.anova <- as.data.frame(SP_SP_dist_table %>% group_by(SITE) %>% summarise(tr_length=sum(length_from_prev_SP), exp_group=head(exp_group,1), .groups="drop")) %>% anova_test(tr_length ~ exp_group)
pwc <- tukey_hsd(as.data.frame(SP_SP_dist_table %>% group_by(SITE) %>% summarise(tr_length=sum(length_from_prev_SP), exp_group=head(exp_group,1), .groups="drop")), tr_length ~ exp_group)
pwc <- pwc %>% add_xy_position(x = "exp_group")

pwc$xmin <- as.numeric(sapply(pwc$group1, switch, "Perennial"=1, "Summer dry"=2, "More seasons dry"=3))
pwc$xmax <- as.numeric(sapply(pwc$group2, switch, "Perennial"=1, "Summer dry"=2, "More seasons dry"=3))

p <- ggplot(as.data.frame(SP_SP_dist_table %>% group_by(SITE) %>% summarise(tr_length=sum(length_from_prev_SP), exp_group=head(exp_group,1), .groups="drop")), aes(x=factor(exp_group, levels=c("Perennial", "Summer dry", "More seasons dry")), y=tr_length)) + 
  geom_boxplot(aes(fill=exp_group), show.legend = FALSE) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(as.data.frame(SP_SP_dist_table %>% group_by(SITE) %>% summarise(tr_length=sum(length_from_prev_SP), dry_group=head(dry_group,1), .groups="drop"))$tr_length)+0.3, by=0.4, length.out=length(y.position))), aes(y=y.position, x=xmin, xend=xmax)) + 
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(as.data.frame(SP_SP_dist_table %>% group_by(SITE) %>% summarise(tr_length=sum(length_from_prev_SP), dry_group=head(dry_group,1), .groups="drop"))$tr_length)+0.3, by=0.4, length.out=length(y.position))), aes(x=xmin, y=y.position, yend=y.position-0.1)) +
  geom_segment(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(as.data.frame(SP_SP_dist_table %>% group_by(SITE) %>% summarise(tr_length=sum(length_from_prev_SP), dry_group=head(dry_group,1), .groups="drop"))$tr_length)+0.3, by=0.4, length.out=length(y.position))), aes(x=xmax, y=y.position, yend=y.position-0.1)) +
  geom_text(data=pwc[which(!pwc$p.adj.signif=="ns"),] %>% mutate(y.position=seq(from=max(as.data.frame(SP_SP_dist_table %>% group_by(SITE) %>% summarise(tr_length=sum(length_from_prev_SP), dry_group=head(dry_group,1), .groups="drop"))$tr_length)+0.3, by=0.4, length.out=length(y.position))), aes(label=p.adj.signif, y=y.position, x=(xmin+xmax)/2,), position=position_nudge(y=0.1)) +
  scale_y_continuous(limits=c(0,NA)) +
  scale_x_discrete(labels=c("Perennial", "Summer-dry", "More-seasons-dry")) +
  scale_fill_viridis(option="E", discrete=TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(margin=margin(0,0,2,0)),
        plot.subtitle = element_text(margin=margin(0,0,0,0)),
        axis.title.y = element_text(size=13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(size=13, color="black")) +
  xlab("") +
  ylab("Trajectory length") +
  labs(title="", 
       subtitle=get_test_label(res.anova, detailed = TRUE))
reduced_dataset_plotlist <- append(reduced_dataset_plotlist, list(p))
names(reduced_dataset_plotlist)[length(reduced_dataset_plotlist)] <- "A"
remove(pwc, res.anova)

reduced_dataset_plotlist <- reduced_dataset_plotlist[order(names(reduced_dataset_plotlist))]
ggarrange(plotlist=reduced_dataset_plotlist, nrow=2, ncol=1, labels=c("A", "B"), font.label = list(size=25))
#ggsave("REV_FINAL_Fig_S3.png", height = 8, width=8, bg="white")

ggplot(na.omit(DD_365_FINAL %>% left_join(grouping_table[, c(1,3)], join_by(SITE==SITE))), aes(x=factor(dry_group), y=DD_365, fill=dry_group, group=reorder(factor(SITE), DD_365))) +
  geom_boxplot(width=1, position = position_dodge2(preserve = "single")) +
  geom_text(data=data.frame(na.omit(DD_365_FINAL %>% left_join(grouping_table[, c(1,3)], join_by(SITE==SITE))) %>% group_by(dry_group, SITE) %>% summarise(mean_DD=mean(DD_365), .groups="drop") %>% arrange(mean_DD, .by_group=TRUE)), 
            aes(label=SITE, x=factor(dry_group), y=-6, group=reorder(factor(SITE), mean_DD)), 
            position = position_dodge2(width=1, preserve = "single"), size=3) +
  annotate("text", label="SITE:", x=0.51, y=-6, size=3) +
  scale_fill_manual(values=viridis(3)) + 
  scale_y_continuous(limits=c(-6,NA)) +
  theme_bw() +
  theme(axis.title.y = element_text(size=13),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        panel.grid.major.x = element_blank()) +
  ylab("Number of dry days\nwithin the 365 days before sampling") +
  xlab("") +
  labs(fill="Duration of drying")
#ggsave("REV_FINAL_Fig_S1.png", height = 6, width=12, bg="white")