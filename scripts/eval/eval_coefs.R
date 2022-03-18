# for all model_versions
# save significant as well
library(dplyr)
dir<- 'output/[name]/models/gag/categorical/final_models/hlaboth_clin/'
dir<-'output/[name]/models/gag/categorical/final_models/hla1_clin/'
suffix <- 'saav_freq_0.01/hla_freq_0.01/cd4_cutoff500/seed_100/full/fit_brms_phylo-raxml_seqtype-nt_glmm_coeffs.csv'
saav_ids<- list.files(path = dir)
coeffs<-lapply(saav_ids, function(s){
	fname <- file.path(dir, s, suffix)
	df <-read.csv(fname)
	top10<-df %>% filter(!grepl('Intercept', X)) %>%  arrange(desc(abs(Estimate))) %>% head(10)
	return(top10)
})
names(coeffs)<-saav_ids
coeffs_df<-bind_rows(coeffs, .id = 'site')
write.csv(coeffs_df, row.names=FALSE, file = 'output/[name]/results/sign_coeffs.csv')
sign.coeffs<-coeffs_df %>% filter(l.95..CI * u.95..CI > 0)
