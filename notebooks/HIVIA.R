library(dplyr)
library(tidyr)
library(ggplot2)
library(PRROC)
library(stringr)
library(xtable)
library(psych)
library(optparse)

set.seed(100)
option_list=list(
  make_option(c('--adaptscore_training'), type="character", default=NULL, help="CSV file with predicted adapted score for the training samples via CV and HLAboth model", metavar="FILE"), 
  make_option(c('--adaptscore_training_seed2'), type="character", default=NULL, help="CSV file with predicted adapted score for the training samples via CV with different seed", metavar="STRING"), 
  make_option(c('--adaptscore_training_hla1_without_clin'), type="character", default=NULL, help="CSV file with predicted adaptation score for training samples via CV and HLA1 model without any clinical factors", metavar="FILE"),
  make_option(c('--adaptscore_leftout'), type="character", default=NULL, help="CSV file with predicted adapted score for the leftout samples", metavar="STRING"), 
  make_option(c('--adaptscore_acute'), type="character", default=NULL, help="CSV file with predicted adaptation score for acute data set", metavar="FILE"),
  make_option(c('--saav_preds_training'), type="character", default=NULL, help="RDS file with conditional probabilities for every site in the training set via CV and HLAboth model", metavar="FILE"),
  make_option(c("--coreceptor_raw"), type="character", default=NULL, help="CSV file with output from geno2pheno_coreceptor_tool", metavar="STRING"), 
  make_option(c("--coreceptor_processed"), type="character", default=NULL, help="CSV file with processed output from geno2pheno_coreceptor_tool", metavar="STRING"), 
  make_option(c("--hla"), type="character", default=NULL, help="CSV file with HLA alleles for every sample.", metavar="STRING"), 
  make_option(c("--clinical"), type="character", default=NULL, help="CSV file with clinical information", metavar="STRING"), 
  make_option(c("--sign_coeffs_dir"), type="character", action= "NULL", default=NULL, help="Directory where the coefficients of each per site model is stored"), 
  make_option(c("--sign_coeffs_suffix"), type="character", action= "NULL", default=NULL, help="(Suffix) Path where the coefficients of each per site model is stored"),
  make_option(c("--output_dir"), type="character", action= "NULL", default=NULL, help="Path where results for the manuscript are stored")

)
#functions
convert_y <-function(y, pred_saav){
  # same length?
  # all y the same
  
  if (unique(y) %in% pred_saav  || !("OTHER" %in% pred_saav)){
    return(y)
  } else {
    return(rep("OTHER", length(y)))
  }
}

# Combines the leftout predictions with the CV predictions.
# @param df1 first dataset
# @param df2 second dataset; cv dataset with column rand_id
# @param name1 name of the first dataset
# @param name2 name of the second dataset
# @return A dataframe with the combined data and a new column dataset

combine_leftout_cv <- function(df1, df2, name1, name2){

	df <- bind_rows(list(leftout =df1, cv = df2 %>% filter(is.na(rand_id))%>% select(-rand_id)) , .id = 'dataset')
	return(df)
}

get_mean_adaptscore <-function(df){
	stopifnot(c("dataset", "adapt_trans") %in% colnames(df)) 
	mu <-df %>% 
		group_by(dataset) %>% 
		# adapt trans - transformed adapt score
		summarize(grp.mean = mean(adapt_trans))
	return(mu)
}
# set option parser
opt_parser<-OptionParser(option_list=option_list)
# read arguments
opt<-parse_args(opt_parser)

now <- format(Sys.time(), "%Y%m%d_%H%M%S")
#create output directory if not existing yet
output_dir<-file.path(opt$output_dir, now)
ifelse(!dir.exists(output_dir), dir.create(output_dir, recursive = TRUE), FALSE)


# open log file
sink(file.path(output_dir, 'log.txt'))
# LOAD INPUT DATA

adapt_leftout <- read.csv(opt$adaptscore_leftout)
adapt_cv_42 <- read.csv(opt$adaptscore_training_seed2)
adapt_cv<-read.csv(opt$adaptscore_training)
adapt_cv_hla_without_clin <- read.csv(opt$adaptscore_training_hla1_without_clin) %>% filter(is.na(rand_id))
pred_hlaboth<-readRDS(opt$saav_preds_training)
acute<-read.csv(opt$adaptscore_acute)  %>% mutate(patient_id = as.character(patient_id))
coreceptor_raw<-read.csv(opt$coreceptor_raw)
coreceptor_raw <- coreceptor_raw %>% rename(patient_id = header)
coreceptor_df<-read.csv(opt$coreceptor_processed)
hla_df<-read.csv(opt$hla)
clinical_df<-read.csv(opt$clinical)


#'
#' # Characteristics of adaptation score
#' Comparison of adaptation scores between chronic and less-chronic data set (cutoff 500 cd4count). 
#'
## cd4low vs cd4high
# bind cd4high and cd4low dataset together, remove the random predictions from the cd4low data set
adapt_df<-combine_leftout_cv(df1=adapt_leftout, df2=adapt_cv)
# compute mean adaptation for each dataset


mu<-get_mean_adaptscore(adapt_df)
print("Mean adaptation of training:")
mu %>% filter(dataset =="cv")
print("Mean adaptation of leftout:")
mu %>% filter(dataset == "leftout")
# plot
p<- adapt_df %>% 
	# adaptation score stratified on data_fit
	ggplot(aes(x = adapt_trans, fill = dataset)) + 
	# as histogram
		geom_histogram(alpha = 0.5, position = 'dodge') + 
	# add a line 
		geom_vline(mu, mapping = aes(xintercept=grp.mean, color = dataset), linetype = 'dashed') + 
	# flip coordinates
		coord_flip() +  
		theme_minimal() +
	# put legends on the top
#		theme(legend.position="top") +
	# label x_axis
		labs(x = 'adaptation score') + 
		theme(legend.position = "none",text=element_text(size=20),
         axis.title.y=element_text(size=20),
         axis.title.x=element_text(size=20))
	# relabel 
ggsave(filename = file.path(output_dir, '/chronic_lowCD4_vs_chronic_highCD4.pdf'), plot = p, width = 5.33, height=4, unit= 'in')
p_val = wilcox.test(adapt_df %>% filter(dataset == 'cv') %>% select(adapt_trans) %>% pull, adapt_df %>% filter(dataset == 'leftout') %>% select(adapt_trans) %>% pull, alternative = 'greater',paired =FALSE)$p.value

print('Do patients with lower CD4 count have more adapted viruses than patients with higher CD4 count?')
print(p_val)
print('(unpaired Wilcoxon test)')

#'
#' Chronic variants (low CD4 count) are not more adapted than less chronic variants (cd4 >500). (p-value = `r format.pval(p_val)`)
#' 
#' Comparison between adaptation score for random and true data set
rand_mean_df<-adapt_cv%>% filter(!is.na(rand_id))  %>% group_by(patient_id) %>% summarize(adapt_trans = mean(adapt_trans))
adapt_df<-bind_rows(list(rand = rand_mean_df, cv= adapt_cv %>% filter(is.na(rand_id))),.id = "dataset")
mu <-adapt_df %>% group_by(dataset) %>% summarize(grp.mean = mean(adapt_trans))
print("Mean adaptation of training:")
mu %>% filter(dataset =="cv")
print("Mean adaptation of random data set:")
mu %>% filter(dataset == "rand")

p<- adapt_df %>% 
	# adaptation score stratified on data_fit
	ggplot(aes(x = adapt_trans, fill = dataset)) + 
	# as histogram
		geom_histogram(alpha = 0.5, position = 'dodge') + 
	# add a line 
		geom_vline(mu, mapping = aes(xintercept=grp.mean, color = dataset), linetype = 'dashed') + 
	# flip coordinates
		coord_flip() +  
		theme_minimal() +
	# put legends on the top
#		theme(legend.position="top") +
	# label x_axis
		labs(x = 'adaptation score') +
		theme(legend.position = "none",text=element_text(size=20),
         axis.title.y=element_text(size=20),
         axis.title.x=element_text(size=20))
	# relabel 
#		scale_fill_manual(name = 'dataset', 
#			breaks= c('chronic_lowCD4', 'random'), 
#			values = c('blue', 'red'))
ggsave(filename = file.path(output_dir,'chronic_lowCD4_vs_random.pdf'), plot = p, width = 5.33, height=4, unit= 'in')
p_val = wilcox.test(adapt_df %>% filter(dataset == 'cv') %>% arrange(patient_id) %>% select(adapt_trans) %>% pull, adapt_df %>% filter(dataset == 'rand') %>% arrange(patient_id) %>% select(adapt_trans) %>% pull, alternative = 'greater',paired =TRUE)$p.value

## Comparison of adaptation scores between true and random HLA
# same computation -different code
df<-left_join(rand_mean_df,adapt_cv %>% filter(is.na(rand_id)) %>% select(patient_id, adapt_trans), by = 'patient_id')
p_val <-df %>% summarize(p_val = wilcox.test(adapt_trans.x, adapt_trans.y, paired = TRUE, alternative = 'less')$p.value)

#'  Adaptation in chronic_lowCD4 is higher than of random dataset (p-value `r format.pval(p_val)`)
#'
print('Are viruses more adapted to the real HLA host profile compared to a random created HLA profile? ') 
print(p_val)
print('(paired Wilcoxon test)')

#boxplot all together
rand_mean_df<-adapt_cv%>% filter(!is.na(rand_id))  %>% group_by(patient_id) %>% summarize(adapt_trans = mean(adapt_trans))
adapt_df<-bind_rows(list( 'chronic_highCD4'= adapt_leftout, 'chronic_lowCD4'= adapt_cv %>% filter(is.na(rand_id)), random = rand_mean_df),.id = "dataset")
p<-adapt_df %>% 
	ggplot(aes(x=dataset, y  = adapt_trans)) + 
	geom_boxplot() +
	labs(y = 'adaptation score', x = 'dataset') +
	theme_minimal() +
	theme(text=element_text(size=20),
         axis.title.y=element_text(size=20),
         axis.title.x=element_text(size=20)) +
	coord_flip()

ggsave(filename = file.path(output_dir, 'boxplot_datasets.pdf'), plot = p, width = 5.33, height=4, unit= 'in')

# Percentage of predicted adapted samples
adapt_df %>% group_by(dataset) %>% mutate(tot_size = n()) %>% filter(adapt_trans>0.1) %>% summarize(n_adapt=n(), tot_size = mean(tot_size)) %>% group_by(dataset) %>% summarize(perc = n_adapt/tot_size)

#' # Datasets - Summary statistics
format_stats <- function(x){
	mean <-round(mean(x, na.rm=TRUE),2)
	sd <-round(sd(x, na.rm = TRUE),2)
	res<-paste0(mean, '±', sd)
	return(res)

}
adapt_df<-bind_rows(list('chronic_highCD4'=adapt_leftout, 'chronic_lowCD4'= adapt_cv %>% filter(is.na(rand_id))), .id = 'dataset')
comb_df<-adapt_df %>% 
  left_join(coreceptor_df, by = "patient_id") %>%
  left_join(coreceptor_raw, by = "patient_id") %>%
  left_join(clinical_df, by = "patient_id") %>%
  left_join(hla_df, by= 'patient_id')
numeric_stats<-comb_df %>% group_by(dataset) %>% select(age, sex, coarse_race, cd4_count, log_vl, adapt_trans) %>% summarise_if(is.numeric, ~format_stats(.))
sex_stats<-comb_df %>% 
	group_by(dataset,sex) %>%
	summarize(sex_count = n()) %>%
	group_by(dataset) %>%
	mutate(sex_count = sex_count/sum(sex_count)) %>%
	spread(sex, sex_count)
	
ethnic_stats<- comb_df %>% 
	group_by(dataset,coarse_race) %>%
	summarize(ethnicity_count = n()) %>%
	group_by(dataset) %>%
	mutate(ethnicity_count = ethnicity_count/sum(ethnicity_count)) %>%
	spread(coarse_race, ethnicity_count)

coreceptor_stats<-comb_df %>% 
	group_by(dataset,coreceptor_EU) %>%
	summarize(count = n()) %>%
	group_by(dataset) %>%
	mutate(count = count/sum(count)) %>%
	spread(coreceptor_EU, count)
	
stats_df<-left_join(numeric_stats, sex_stats, by='dataset') %>%
	left_join(ethnic_stats, by = 'dataset') %>%
	left_join(coreceptor_stats)

coreceptor_adapt<-comb_df %>% group_by(dataset, coreceptor_EU) %>% summarize(adapt = paste0(round(mean(adapt_trans),2), '±', round(sd(adapt_trans),2)), max = max(adapt_trans))
#' Summary statistics 
print('summary statistics')
stats_df %>% xtable
#t(stats_df) %>% xtable
print('average adaptation stratified according to coreceptor usage')
#' Average adaptation stratified according to coreceptor usage and dataset
coreceptor_adapt %>% xtable



adapt_robust<-left_join(adapt_cv %>% filter(is.na(rand_id)) %>% select(patient_id, adapt_trans), adapt_cv_42 %>% filter(is.na(rand_id)) %>% select(patient_id, adapt_trans), by = 'patient_id')

print('correlation betweeen two different CV runs (seed change)')
adapt_robust %>% summarize(cor = cor(adapt_trans.x, adapt_trans.y))




avg_sd <- adapt_robust %>% group_by(patient_id) %>% mutate(sd_= sd(c(adapt_trans.x, adapt_trans.y)))  %>% ungroup %>% summarize(avg_sd = mean(sd_))
print('average standard deviation for adaptation score is')
print(avg_sd)


#' # SAAV performance


pred_df_hlaboth <- pred_hlaboth %>% 
				filter(is.na(rand_id)) %>%
        ungroup() %>% 
        gather(key="pred_saav", value = "pred", starts_with("P.Y")) %>%
        mutate(pred_saav=gsub("\\.$", "", gsub("P\\.Y\\.\\.\\.", "", pred_saav))) %>%
        mutate(pred_saav = gsub("\\.", "-", pred_saav)) %>% 
        group_by(patient_id, saav_id,model) %>%
        filter(!is.na(y), !is.na(pred)) %>%
        mutate(y = convert_y(y, pred_saav)) %>% 
        ungroup

perf_df_hlaboth <- pred_df_hlaboth %>% 
  group_by(protein, model, saav_id, pred_saav) %>% 
  mutate(y_bin = ifelse(y ==pred_saav, 1, 0)) %>% 
  filter(!is.na(y_bin), !is.na(pred))%>%
  summarize(PRROC = pr.curve(scores.class0=pred[y_bin==1], scores.class1=pred[y_bin==0])$auc.davis.goadrich, pr_bline = sum(y_bin)/n())

saveRDS(perf_df_hlaboth, file.path(output_dir, 'per_site_models_preds_performance.RDS'))
#averaged over all SAPs per site
#' Top ten sites with well performing models
perf_sites <- perf_df_hlaboth %>% group_by(saav_id) %>% filter(model == 'full') %>% summarize(mean_PRROC = mean(PRROC), mean_pr_bline = mean(pr_bline)) %>% filter(mean_PRROC > mean_pr_bline) %>% mutate(diff = mean_PRROC - mean_pr_bline) %>% arrange(desc(diff)) %>% head(10)

total_site_models<-perf_df_hlaboth %>% group_by(saav_id) %>% filter(model == 'full') %>% summarize(mean_PRROC = mean(PRROC), mean_pr_bline = mean(pr_bline)) %>% dplyr::count() %>% pull

perf_site_models<-perf_df_hlaboth %>% group_by(saav_id) %>% filter(model == 'full') %>% summarize(mean_PRROC = mean(PRROC), mean_pr_bline = mean(pr_bline)) %>% filter(mean_PRROC > mean_pr_bline) %>% dplyr::count() %>% pull

print(paste0('Out of ', total_site_models, ' site models ', perf_site_models, ' have higher performance than the Precision recall baseline model'))


print('Top ten sites with best performing models (averaged over all SAPs per site')
perf_sites

# over all variants
#' Top ten saav models (with higher performance than baseline)
total_saav_models<-perf_df_hlaboth %>% ungroup %>% filter(model=="full") %>% dplyr::count() %>% pull
perf_saav_models<-perf_df_hlaboth %>% ungroup %>% filter(model =='full', PRROC > pr_bline) %>% dplyr::count() %>% pull
print(paste0('Out of ', total_saav_models, ' saav models ', perf_saav_models, ' have higher performance than the Precision recall baseline model'))
perf_preds <- perf_df_hlaboth %>% ungroup %>% filter(model =='full', PRROC > pr_bline) %>% mutate(diff = PRROC - pr_bline) %>% arrange(desc(diff)) %>% head(10) %>% select(-protein, -model, -diff)

print('Top ten saav models (wrt precision_recall_baseline)')
perf_preds


compute_odds_df<-function(df){
	odd_df<- df %>%
     ungroup() %>% 
     gather(key="pred_saav", value = "pred", starts_with("P.Y")) %>%
     mutate(pred_saav=gsub("\\.$", "", gsub("P\\.Y\\.\\.\\.", "", pred_saav))) %>%
     mutate(pred_saav = gsub("\\.", "-", pred_saav)) %>%
     group_by(patient_id, saav_id, rand_id,model) %>%
     filter(!is.na(y), !is.na(pred)) %>%
     mutate(y = convert_y(y, pred_saav)) %>% 
     ungroup %>% 
     filter(y == pred_saav) %>% 
     tidyr::spread(model, pred) %>% 
	    group_by(patient_id, rand_id, saav_id) %>% 
	    #summarize(odd = ifelse(abs(full -null_clin_rand)>0.1, full/null_clin_rand,1)) %>%
	    summarize(odd = max(full, 0.0001)/max(null_clin_rand, 0.0001), full=full,null = null_clin_rand, y = y)
	return(odd_df)
}
# Retrieve patient_id with highest adaptation
pid_highest_adapt<-adapt_cv %>% filter(is.na(rand_id)) %>% arrange(desc(adapt_trans)) %>% head(1) %>% select(patient_id) %>% pull
odds_df <-compute_odds_df(pred_hlaboth %>% filter(patient_id %in% pid_highest_adapt, is.na(rand_id)))
saveRDS(odds_df, file=file.path(output_dir, "odds.rds"))

## MAKE LOGO
SIGN_CUTOFF = 0.1
ADAPT_CUTOFF = 0
SEQ_NUMS = 100

logo_df<- odds_df %>% 
# select patient with highest adaptation
filter(patient_id %in% pid_highest_adapt, is.na(rand_id)) %>% 
# shift odds using 1 as cutoff
mutate(odds_shift= odd -1) %>% 
# se;ect only those sites with effect > as variation of 0.1
filter(abs(odds_shift)>SIGN_CUTOFF) %>%
# scale odds by maximum contribution
 mutate(odds_norm = round(odds_shift/max(abs(odds_shift)),2)) %>% 
 # select only those contributing at least 1 %
 filter(abs(odds_norm) >=0.01) %>% 	
 mutate(y_logo = ifelse(odds_norm>ADAPT_CUTOFF, toupper(y), tolower(y))) %>% 
# created signal
mutate(n_gap= SEQ_NUMS-(abs(odds_norm)*SEQ_NUMS), n_y=  abs(odds_norm)*SEQ_NUMS) %>% 
group_by(saav_id) %>% 
mutate(sig=paste0(c(rep("-", (SEQ_NUMS-(abs(odds_norm)*SEQ_NUMS))), rep(y_logo, abs(odds_norm)*SEQ_NUMS)), collapse="")) %>%
ungroup

print("after_logo_df")
logo_ticks <- logo_df %>% select(saav_id) %>% pull
logo_seqs <- logo_df %>% select(sig) %>% pull

logo_align <-lapply(logo_seqs, function(x){unlist(strsplit(x, split=""))})
logo_align <-do.call(cbind, logo_align)

library(seqinr)


fasta_align <-apply(logo_align,1, function(x){paste(x, collapse="")})

write.fasta(as.list(fasta_align), file.out=file.path(output_dir, "logo.fasta"), names = seq(1, SEQ_NUMS))
write(logo_ticks, file=file.path(output_dir, "logo_ticks.txt"))

#TODO call to weblogo
command_prefix <-"weblogo --format PNG --size large < "
input <- file.path(output_dir, "logo.fasta")
output <- file.path(output_dir, "logo.png")
x_ticks<-paste0(logo_ticks, collapse=", ")
command<-paste0(command_prefix, input, " -o ", output, " -A protein --composition 'equiprobable' -X YES --ticmarks 3.0 --errorbars NO --scale-width NO --resolution 800 -P '' --annotate '", x_ticks, "' -Y NO --alphabet 'AaRrNnDdCcEeQqGgHhIiLlKkMmFfPpSsTtWwYyVv*#' --composition equiprobable --color blue 'ARNDCEQGHILKMFPSTWYV*' positive --color orange 'arndceqghilkmfpstwyv#' negative", sep ="")
system(command)


#' # Adaptation score and clinical variables

#' Correlation of adaptation score and VL, CD4 count and coreceptor usage 
comb_df %>%  
select(adapt_trans, log_vl, cd4_count, FPR) %>% 
psych::corr.test(method='spearman')

comb_df %>%  
select(adapt_trans, log_vl, cd4_count, FPR) %>% 
psych::corr.test(method='pearson')

comb_df_new <-comb_df %>% 	mutate (adapt_level = case_when(
		adapt_trans > 0.1 ~ 'adapted',
		adapt_trans < -0.1 ~ 'non-adapted', 
		TRUE ~ 'INTERMEDIATE' ))


var_adapt_plot<-comb_df_new  %>% mutate(
	CD4 = ifelse(cd4_count>200, 'cd4>200', 'cd4<=200')) %>%
	rename(coreceptor = coreceptor_EU) %>%
	select(adapt_trans, CD4, coreceptor) %>%
	gather(condition, level, CD4, coreceptor, -adapt_trans) %>% ggplot(aes(x=condition, y = adapt_trans, col = level)) + 
	geom_boxplot(width= 0.5, position= position_dodge(0.7)) + 
	theme_minimal() + 
	labs( y  ='adaptation score') + 
	theme(axis.title.x = element_blank(), text=element_text(size=18),
         axis.title.y=element_text(size=18))

ggsave(filename = file.path(output_dir, 'adaptation_stratified.pdf' ), plot = var_adapt_plot,width = 5.33, height=4, unit= 'in' )




p_val<-comb_df_new %>% mutate(CD4 = ifelse(cd4_count>200, 'cd4>200', 'cd4<=200')) %>% 
summarize(cd4_pval = wilcox.test(adapt_trans[CD4=='cd4<=200'], adapt_trans[CD4=='cd4>200'], alternative='greater')$p.value)
saveRDS(comb_df_new, file.path(output_dir, "comb_df_new.rds"))
n_AIDS<- comb_df_new %>% filter(cd4_count <= 200) %>% dplyr::count() %>% pull
n_nonAIDS <- comb_df_new %>% filter(cd4_count >200) %>% dplyr::count() %>% pull

#' Do AIDS patients (n = `r n_AIDS`) harbor higher adapted viruses compared to non-AIDS (n = `r n_nonAIDS`)? (p-value =  `r format.pval(p_val)`)

print('Do AIDS patients have higher adapted viruses compared to non-AIDS patients?')
print(p_val)

p_val<-comb_df_new %>% 
summarize(vl_pval = wilcox.test(adapt_trans[log_vl<log10(2000)], adapt_trans[log_vl>log10(10000)], alternative='less')$p.value)

n_controller<- comb_df_new %>% filter(log_vl <= log10(2000)) %>% dplyr::count() %>% pull
n_noncontroller <- comb_df_new %>% filter(log_vl >= log10(10000)) %>% dplyr::count() %>% pull


#' Have patients with low viral load less adapted viruses?  (p-value `r format.pval(p_val)`).
print('Have patients with low viral load less adapted viruses?')
print(p_val)

#' However, sample size is pretty low (controller : `r n_controller`, non-controller = `r n_noncontroller`)

p_val<-comb_df_new %>% 
summarize(coreceptor_pval = wilcox.test(adapt_trans[coreceptor_EU == 'X4'], adapt_trans[coreceptor_EU =='R5'], alternative='greater')$p.value)

n_X4<- comb_df_new %>% filter(coreceptor_EU == 'X4') %>% dplyr::count() %>% pull
n_R5 <- comb_df_new %>% filter(coreceptor_GER == 'R5') %>% dplyr::count() %>% pull

#' Are X4 variants (n=`r n_X4`) more adapted than R5 variants (n=`r n_R5`) (p-value = `r format.pval(p_val)`)
print('Are X4 variants more adapted than R5 variants?')
print(p_val)

# update adapt level
perform_R5_matchedPairAnalysis<-function(seed, df){
	set.seed(seed)
	# extract samples with R5 coreceptor, non-adapted, in decreasing cd4 count order for given model
	R5_low_adapt<-df %>% filter(coreceptor_EU =='R5', adapt_level == 'non-adapted') %>% arrange(desc(cd4_count))
	# extract samples with R5 coreceptor, non-adapted, in decreasing cd4 count order for given model
	R5_high_adapt <- df%>% filter(coreceptor_EU=='R5', adapt_level == 'adapted') %>% arrange(desc(cd4_count))
	#matched_ys<-NULL
	matched_y<-NULL
	matched_fpr<-NULL
	matched_seeds <- NULL
	matched_x<-NULL 
	# iterate over all low R5 samples
	for (i in 1:nrow(R5_low_adapt)){
		# store current cd4 count
		tmp_cd4_count<-R5_low_adapt[i, "cd4_count"]
		# filter all potential matches from the R5_high dataset that have not been selected yet
		pot_matches<-R5_high_adapt %>% filter(cd4_count <=tmp_cd4_count + 50, cd4_count >= tmp_cd4_count - 50, !(patient_id %in% matched_y))
		# if match was found
		if(nrow(pot_matches)!= 0){
			# store the patient id of first match
			matched_y<-c(matched_y, pot_matches[1, "patient_id"])
			# store the fpr of first match
			matched_fpr<-c(matched_fpr, pot_matches[1, "FPR"])
		} else {
		# add NA
		matched_y<-c(matched_y, NA)
		matched_fpr<-c(matched_fpr, NA)
		}
		# store the R5 low adapt patient id
		matched_x<-c(matched_x, R5_low_adapt[i,'patient_id'])
	}
	matched_df<- data.frame( matched_id = matched_y, matched_fpr = matched_fpr, patient_id = matched_x)
	R5FPR_df<-R5_low_adapt %>% left_join(matched_df, by = 'patient_id') %>% filter(!is.na(matched_fpr)) %>% select(FPR, matched_fpr) %>% gather(variable, value) %>% mutate(adapt_level= ifelse(grepl('matched', variable), 'adapted', 'non-adapted')) %>% mutate(variable = 'R5-FPR') 
	matched_df <- R5_low_adapt %>% left_join(matched_df, by = 'patient_id')
	return(list("matched_df" = matched_df, "R5FPR_df"=R5FPR_df))
}

R5FPR_res<-perform_R5_matchedPairAnalysis(seed=100, comb_df_new)
R5FPR_df <-R5FPR_res[["R5FPR_df"]]

#' # Comparison of clinical variables wrt adaptation level
n_adapt<- comb_df_new %>% group_by(adapt_level) %>% dplyr::count()
#' There are `r n_adapt %>% filter(adapt_level == 'adapted') %>% select(n) %>% pull` adapted viruses and `r n_adapt %>% filter(adapt_level == 'non-adapted') %>% select(n) %>% pull` non-adapted viruses

p<-comb_df_new %>%
		filter(adapt_level !='INTERMEDIATE') %>%
		rename(CD4 = cd4_count, VL = log_vl) %>% 
		select(adapt_level, CD4, VL, FPR) %>% 
		gather(variable, value, -adapt_level ) %>%
		bind_rows(R5FPR_df) %>% 
		ggplot(aes(x=variable, y = value, col = adapt_level)) + geom_boxplot() + 
		facet_wrap(~variable, scale='free') + 
		theme_minimal()+
		theme(legend.position="top", text=element_text(size=18), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x  = element_blank(), axis.title.y =element_blank())
		
ggsave(file.path(output_dir, 'boxplot_adapt_level_new.pdf'), p, width = 5.33, height=4, unit= 'in' )



matched_df<- R5FPR_res[["matched_df"]]

## Lower FPR in higher adapted R5 
 	p_val <-matched_df %>% summarize(p_val=wilcox.test(matched_fpr,FPR, alternative = "less", paired=TRUE)$p.value)



#' Do adapted R5-variants have lower FPR compared to non-adapted R5-variants (n = `r n`) (p-value = `r format.pval(p_val)`) 
print('Do adapted R5-variants have lower FPR compared to non-adapted R5-variants')
print(p_val)

p_val<-comb_df_new %>% summarize(FPR_pval = wilcox.test(FPR[adapt_level == 'adapted'], FPR[adapt_level == 'non-adapted'], alternative = 'less')$p.value)

#' Do adapted viruses have lower FPR than non-adapted? (p-value `r format.pval(p_val)`)
print('Do adapted viruses have lower FPR compared to non-adapted viruses')
print(p_val)

p_val<-comb_df_new %>% summarize(cd_count = wilcox.test(cd4_count[adapt_level == 'adapted'], cd4_count[adapt_level == 'non-adapted'], alternative = 'less')$p.value)

#' Do adapted viruses have lower cd4_count than non-adapted? (p-value `r format.pval(p_val)`)
print('Do adapted viruses have lower CD4 count than non-adapted viruses?')
print(p_val)

p_val<-comb_df_new %>% summarize(log_vl_pval = wilcox.test(log_vl[adapt_level == 'adapted'], log_vl[adapt_level == 'non-adapted'], alternative = 'greater')$p.value)

#' Do adapted viruses have higher viral loads than non-adapted? (p-value `r format.pval(p_val)`)
print('Do adapted viruses have higher viral loads than non-adapted?')
print(p_val)

#' # Adaptation scores (HLA1 based) of different datasets

p_val<-wilcox.test(adapt_cv_hla_without_clin %>% select(adapt_trans) %>% pull, acute %>% select(adapt_trans) %>% pull, alternative = 'greater')$p.value
	
# Acute infected patients have less adapted viruses as chronic patients (p-value `r format.pval(p_val)`)
print('chronic infected patients have more adapted viruses compared to acute-infected patients?')
print(p_val)

comb_df <- bind_rows(list('acute'= acute, 'chronic' = adapt_cv_hla_without_clin), .id = 'dataset')
mu <-comb_df %>% group_by(dataset) %>% summarize(grp.mean = mean(adapt_trans))
print("Mean adaptation of hla only model:")
mu %>% filter(dataset =="chronic")
print("Mean adaptation of acute:")
mu %>% filter(dataset == "acute")
p<-comb_df %>% ggplot(aes(x = adapt_trans, fill = dataset)) + 
	# as histogram
		geom_histogram(alpha = 0.5, position = 'dodge') + 
	# add a line 
		geom_vline(mu, mapping = aes(xintercept=grp.mean, color = dataset), linetype = 'dashed') + 
	# flip coordinates
		coord_flip() +  
		theme_minimal() +
	# put legends on the top
#		theme(legend.position="top") +
	# label x_axis
		labs(x = 'adaptation score') + 
		theme(legend.position = "none",text=element_text(size=20),
         axis.title.y=element_text(size=20),
         axis.title.x=element_text(size=20))

ggsave(filename = file.path(output_dir, 'acute_vs_chronic.pdf'), plot = p, width = 5.33, height=4, unit= 'in' )

saav_ids<- list.files(path = opt$sign_coeffs_dir, pattern="^[0-9]*{3}$")
coeffs<-lapply(saav_ids, function(s){
	fname <- file.path(opt$sign_coeffs_dir, s, opt$sign_coeffs_suffix)
	df <-read.csv(fname)
	top10<-df %>% filter(!grepl('Intercept', X)) %>%  arrange(desc(abs(Estimate))) %>% head(10)
	return(top10)
})
names(coeffs)<-saav_ids
coeffs_df<-bind_rows(coeffs, .id = 'site')
sign.coeffs<-coeffs_df %>% filter(l.95..CI * u.95..CI > 0)

print('Top  significant coefficients per site model')
print(sign.coeffs)

## compare heterologous and autologous

# for each viral sequence store the amino acid (y) at each site of interest (saav). 
# select all true samples
viral_ys<-pred_hlaboth %>% 
	filter(is.na(rand_id), model =="full") %>%
	# for each virus/patient
	group_by(patient_id) %>% 
	# select only the id, the site, and the amino acid (drop the predictions)
	select(patient_id, saav_id, y) %>% 
	# order wrt to patient and site
	arrange(patient_id, saav_id)%>% 
	# rename patient id to viral id
	rename(viral_id=patient_id)

# create new data frame where each patient hla profile (and predictions) is combined with all viruses 
new_df<-pred_hlaboth %>% 
	#take only true samples	
	filter(is.na(rand_id)) %>%
	# drop the y column containing the amino acid of the autologous virus
	select(-one_of("y")) %>% 
	# add all viral information for each site
	left_join(viral_ys, by="saav_id")

# compute adaptation for autologous and heterologous viruses
adapt_score<-new_df %>%
      # remove any existing groupings as a precaution 
      ungroup() %>% 
      # convert prediction outcome from wide format to long format. 
      # column pred contains the predicted likelihoods
      # column pred_saav contains the corresponding amino acid at this position
      # There are many possible amino acids and corresponding likelihood predictions for the same site
      gather(key="pred_saav", value = "pred", starts_with("P.Y")) %>%
      # remove dots
      mutate(pred_saav=gsub("\\.$", "", gsub("P\\.Y\\.\\.\\.", "", pred_saav))) %>%
      # substitute remaining dot for dash
      mutate(pred_saav = gsub("\\.", "-", pred_saav)) %>%
      # for each patient, variant site, virus and model (full and null)
      group_by(patient_id, viral_id, saav_id, model) %>%
      # keep only those where there is true information about the amino acid at the variant site and there is a prediction
      filter(!is.na(y), !is.na(pred)) %>%
      # convert the true outcome if necessary to OTHER
      mutate(y = convert_y(y, pred_saav)) %>% 
      # ungroup
      ungroup %>% 
      # keep only those predictions where we have predictions for the true amino acid
      filter(y == pred_saav) %>% 
      # convert long format to wide s.t. we have new columns full and null_clin_rand
      tidyr::spread(model, pred) %>% 
      # for each patient, for each virus, for each site
	   group_by(patient_id, viral_id, saav_id) %>% 
#	    summarize(odd = ifelse(abs(full -null_clin_rand)>0.1, full/null_clin_rand,1)) %>%
      # compute the odds between the likelihoods of the full and the null model
      summarize(odd = max(full, 0.0001)/max(null_clin_rand, 0.0001)) %>%  
      # ungroup
	    ungroup %>%
	    # group by sample -> iterate over all variant sites
	    group_by(patient_id, viral_id) %>%
	    # keep only odds which are not NA
	    # TODO: check and warn if there are NA odds
	    filter(!is.na(odd)) %>%
	    # multipliy the odds over all variant sites for a sample
	    summarize(adapt = prod(odd)) %>%
	    # transform the adaptation score
	    mutate(adapt_trans = atan(log(adapt))*2/pi) %>%
	    # sort in decreasing order
	    arrange(desc(adapt))
		
# data frame to compare autologous and heterologous adaptation	
auto_vs_hetero<- adapt_score %>% 
	# order according to patient
	arrange(patient_id) %>% 
	# for each patient
	group_by(patient_id) %>% 
	# summarize adaptation score for autologous (patient_id equals viral_id) and heterologous viruses (mean)
	summarize(heterologous=mean(adapt_trans[-(patient_id==viral_id)]), autologous =adapt_trans[patient_id==viral_id]) %>%
	# transfer to long format
	gather(key="virus", value="adapt_trans", ends_with("logous"))
	
saveRDS(auto_vs_hetero, file=file.path(output_dir, "auto_vs_hetero.RDS"))
# compute paired wilcoxon test comparing the adaptation per patient from autologous and heterologous 
auto_vs_hetero %>% summarize(p_val = wilcox.test(adapt_trans[virus=="heterologous"], adapt_trans[virus=="autologous"], paired=TRUE, alternative="less")$p.value)

mu <-auto_vs_hetero %>% group_by(virus) %>% summarize(grp.mean = mean(adapt_trans)) 
print("Mean adaptation of autologous:")
mu %>% filter(virus =="autologous")
print("Mean adaptation of heterologous:")
mu %>% filter(virus == "heterologous")
p<- auto_vs_hetero %>% 
	arrange(desc(adapt_trans)) %>%
	# adaptation score stratified on data_fit
	ggplot(aes(x = adapt_trans, fill = virus)) + 
	# as histogram
		geom_histogram(alpha=0.5, position = "dodge") +  
			# add a line 
		geom_vline(mu, mapping = aes(xintercept=grp.mean, color = virus), linetype = 'dashed') + 
	# flip coordinates  
		theme_minimal() +
	# put legends on the top
#		theme(legend.position="top") +
	# flip coordinates
		coord_flip() +  
	# label x_axis
		labs(x = 'adaptation') + 
		theme(legend.position = "none",text=element_text(size=20),
         axis.title.y=element_text(size=20),
         axis.title.x=element_text(size=20))
		 
ggsave(filename = file.path(output_dir, 'auto_vs_hetero.pdf'), plot = p, width = 5.33, height=4, unit= 'in' )		 
			
sink()

