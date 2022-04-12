# HIV immunoadaptation (HIVIA) tool
This is the code basis for the paper 'Insights to HIV-1 coreceptor usage by estimating HLA adaptation using Bayesian generalized linear mixed models'. 
Since HLA information contains private sensitive data that we are not allowed to publish, we can't publish neihter the training HLA data nor the final models. 

## Content

### Snakemake pipeline
The repository contains
* a pipeline for learning HLA adaptation based on the p24 protein and the HLA I and II alleles, as well as clinical factors like age, sex, and ethnicity 
```
 ./Snakefile
```
 
* a prediction pipeline for predicting the HLA adaptation based on the HLA profile and the p24 sequence, as well as clinical factors. 
```
./Snakefile_predict
```

The folder `./rules` contains the different parts of the pipeline. 
The folder `./scripts` contains the source code for the rules. 
### Config files
For each pipeline there is a corresponding config file: 
```
./config.yml
./config_predict.yml
```
### Environments
All environments are stored in the directory `./envs`. 
The conda environment is defined in the file 
```
./envs/conda/HIVIA.yml
```
The Snakemake environment is defined in the directory
```
./envs/snakemake/
```
There are different environments depending on the computational environment (GRID Engine, SLURM, cluster, etc. )

### Run pipepline

The command to run the pipeline (on a cluster) is as follows
```
snakemake -n --configfile config.yml --profile envs/snakemake/[cluster] --nolock

```
### Code for manuscript figures and tables
The file `./notebooks/HIVIA.R` contains the code to produce all figures and tables and p-values in the paper. 


## Key idea
The key elements of this project is to compute the adaptation score of a viral sequence based on its HLA profile, age, sex, and ethnicity. 
The adaptation score is computed as follows: 


For each frequent variant site s<sub>i</sub>, two Bayesian generalized linear mixed models are computed: 
* HLA model
  + age
  + sex
  + ethnicity
  + phylogeny
  + HLA information (binarized)
* baseline model
  + age
  + sex
  + ethnicity
  + phylogeny

The modelling is performed with the script `./scripts/modeling/fit_glmm.R`.
The computation of the adaptation score is performed in the script `./scripts/eval/compute_adaptation.R.`
The call for fitting the HLA model using Bayesian GLMM is the following: 
```
fit<- brm(formula = 'y ~ HLA_allele_1 + ... + HLA_allele_i + ... HLA_allele_n + age + sex + ethnicity + (1|patient_id)', cov_ranef = list(patient_id = ape::vcv.phylo(tree)) , data = df, family = 'categorical', prior = c(set_prior('horseshoe(df = 1, par_ratio = 0.01)', class="b")), seed = 100, control = list(adapt_delta = 0.99))
```

