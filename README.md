# HIV immunoadaptation (HIVIA) tool
This is the code basis for the paper 'Insights to HIV-1 coreceptor usage by estimating HLA adaptation using Bayesian generalized linear mixed models'. 
Since HLA information is private sensitive data, which we are not allowed to publish, we can't publish the training HLA data nor the final models. 

The code contains
* a pipeline for learning HLA adaptation based on the p24 protein and the HLA I and II alleles, as well as clinical factors like age, sex, and ethnicity (
* ```
* ./Snakefile
* ```)
* a prediction pipeline for predicting the HLA adaptation based on the HLA profile and the p24 sequence, as well as clinical factors. 
(```
./Snakefile_predict
```

For each pipeline there is a corresponding config file: 
```
./config.yml
./config_predict.yml
```


The conda environment is defined in the file 
```
```
The Snakemake environment is defined in the file
```
```

The command to run the pipeline is as follows
```
snakemake --config_files --profile 

```

