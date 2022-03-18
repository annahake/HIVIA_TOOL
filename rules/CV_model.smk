# Copyright 2019 Anna Hake.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms. 

  
rule stratified_kCV_folds:
  input: 
    df = rules.select_variables.output
  output:
    temp(expand(['{{run_id}}/data/processed/datasets/{{protein}}/{{outcome_type}}/{{model_version}}/{{saav}}/saav_freq_{{saav_cutoff}}/cd4_cutoff{{cd4}}/seed_{{seed}}/fold_{fold}/train.rds', '{{run_id}}/data/processed/datasets/{{protein}}/{{outcome_type}}/{{model_version}}/{{saav}}/saav_freq_{{saav_cutoff}}/cd4_cutoff{{cd4}}/seed_{{seed}}/fold_{fold}/test.rds'],fold=FOLDS))
  params:
    k = config["modeling"]["kfold"],
    dir = '{run_id}/data/processed/datasets/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/cd4_cutoff{cd4}/seed_{seed}/',
    group = config['modeling']['outcome'],
    script_dir=config['script_dir']
  log: 
    'logs/datasets/{run_id}/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/cd4_cutoff{cd4}/seed_{seed}/stratified_kCV_folds.log'
  shell:
    'Rscript {params.script_dir}/modeling/stratified_kCV_folds.R -i {input.df} -s {wildcards.seed} -d {params.dir} -k {params.k} --group {params.group} --log {log}'

rule CV_get_rand_test: 
  input: 
    rand_df = rules.create_random_dataset.output, 
    test_df = '{run_id}/data/processed/datasets/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/cd4_cutoff{cd4}/seed_{seed}/fold_{fold}/test.rds'
  params:
    script_dir = config['script_dir']
  output: 
    temp('{run_id}/data/processed/datasets/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/cd4_cutoff{cd4}/seed_{seed}/fold_{fold}/test_rand.rds')
  shell: 
    "Rscript {params.script_dir}/modeling/get_rand_test.R -i {input.rand_df} -t {input.test_df} -o {output} --cv cv" 


rule convert_train_hla2binary:
  input:
    '{run_id}/data/processed/datasets/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/cd4_cutoff{cd4}/seed_{seed}/fold_{fold}/train.rds'
  params:
    script_dir=config['script_dir'],
    hla = lambda wildcards: config['hla'][wildcards.model_version]
  output: 
    temp('{run_id}/data/processed/datasets/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/fold_{fold}/train_bin.rds')
  log:
    'logs/make_binary_matrix/{run_id}/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/fold_{fold}/log.log'
  run:
    hla = ",".join(params.hla)
    exec = 'Rscript {params.script_dir}/modeling/convert_train_hla2binary.R'
    exec += ' -i {input} -c {wildcards.hla_cutoff} -o {output}'
    exec += ' --hla ' + hla
    exec += ' -l {log}'
    shell(exec)


rule CV_convert_test_hla2binary: 
  input: 
    df  = '{run_id}/data/processed/datasets/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/cd4_cutoff{cd4}/seed_{seed}/fold_{fold}/{test}.rds'
  params:
    script_dir = config['script_dir'],
    hla = lambda wildcards: config['hla'][wildcards.model_version]
  output:
    temp('{run_id}/data/processed/datasets/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/cd4_cutoff{cd4}/seed_{seed}/fold_{fold}/{test}_bin.rds')
  wildcard_constraints:
    test="(test|test_rand)"
  run:
    hla = ",".join(params.hla)
    exec = 'Rscript {params.script_dir}/eval/convert_test_hla2binary.R -i {input.df} --hla ' + hla + ' -o {output}'
    shell(exec)
    
    
rule CV_fit_SAAV_glmm:  
  input: 
    df = '{run_id}/data/processed/datasets/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/fold_{fold}/train_bin.rds',
    tree = ancient(rules.run_raxml.output)
  params:
    script_dir = config['script_dir'],
    # TODO convert method and prior to wildcards (even if single) set in the config 
    method = 'brms', 
    prior = lambda wildcards: config['modeling']['prior'][wildcards.model]['name'],
    par_ratio = lambda wildcards: config['modeling']['prior'][wildcards.model]['par_ratio'],
    adapt_delta=config["modeling"]["adapt_delta"],
    fixed = lambda wildcards: config['modeling']['fixed'][wildcards.model_version], 
    random = lambda wildcards: config['modeling']['random'],
    outcome = lambda wildcards: config['modeling']['outcome'],
    coeffs = '{run_id}/models/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/fold_{fold}/{model}/fit_brms_phylo-raxml_seqtype-nt_glmm_coeffs.csv', 
    ref_map = rules.get_ref_map.output.ref_map
  benchmark: 
    "benchmarks/{run_id}_CV_fit_SAAV_glmm_{protein}_{outcome_type}_{model_version}_{model}_{saav}_{cd4}_{seed}_{fold}_{saav_cutoff}_{hla_cutoff}_benchmark.txt"
  log:
    'logs/fit_brms/{run_id}/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/fold_{fold}/{model}/log.txt'
  output: 
  	temp('{run_id}/models/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/fold_{fold}/{model}/fit_brms_phylo-raxml_seqtype-nt_glmm.rds')
#  shell: 
#    'Rscript {params.script_dir}/modeling/fit_glmm.R -i {input.df} -t {input.tree} --method {params.method} --model {wildcards.model} -p {params.prior} -s {wildcards.seed} --adapt_delta {params.adapt_delta} -o {output} -l {log}'
  run:
    import os
    my_env=dict(os.environ)
    ENV_VARS=['HOME', 'TMPDIR','TMP','HOSTNAME','PATH','SGE_STDERR_PATH', 'SGE_O_WORKDIR' ,'SGE_O_HOST','SGE_O_HOME', 'PYTHONPATH']
    #with open(log[0], 'w') as l:
    #  [l.write(env_var+":"+ my_env.get(env_var, "")+"\n") for env_var in ENV_VARS]
    for key, value in my_env.items():
      print(key, value)
      
    fixed = ",".join(params.fixed)
#    random = ",".join(params.random)
#    outcome = ",".join(params.outcome)
#    
    exec = "Rscript {params.script_dir}/modeling/fit_glmm.R"
    exec += " -i {input.df} -t {input.tree} --method {params.method}"
    exec += " --model {wildcards.model} -p {params.prior} -s {wildcards.seed}"
    exec += " --adapt_delta {params.adapt_delta} -o {output}"
    exec += " -f " + fixed
    exec += " -r {params.random}"
    exec += " -y {params.outcome}"
    exec += " --coeffs {params.coeffs}"
    exec += " --outcome_type {wildcards.outcome_type}"
    exec += " --ref_map {params.ref_map}"
    exec += " --saav {wildcards.saav}"
    exec += " --par_ratio {params.par_ratio}"
      
    shell(exec)

rule CV_classify_SAAV:
  input:
    fit = rules.CV_fit_SAAV_glmm.output,
    test = '{run_id}/data/processed/datasets/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/cd4_cutoff{cd4}/seed_{seed}/fold_{fold}/test_bin.rds',
    rand = '{run_id}/data/processed/datasets/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/cd4_cutoff{cd4}/seed_{seed}/fold_{fold}/test_rand_bin.rds'
  params: 
    script_dir = config['script_dir']
  output: 
    '{run_id}/results/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/fold_{fold}/{model}/eval.csv'
  shell: 
    "Rscript {params.script_dir}/eval/CV_classify_SAAV.R -f {input.fit} -t {input.test} -r {input.rand} -o {output} --model {wildcards.model} --model_version {wildcards.model_version} -p {wildcards.protein} -s {wildcards.saav} --seed {wildcards.seed}"



def CV_checkpoint_aggregate(wildcards):
  checkpoint_output = checkpoints.create_datasets.get(**wildcards).output[0]
  split_files = expand('{run_id}/results/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/fold_{fold}/{model}/eval.csv', 
    protein=wildcards.protein,
    saav=glob_wildcards(os.path.join(checkpoint_output, '{saav}/full_data/df.rds')).saav,
    cd4 = wildcards.cd4, 
    seed = wildcards.seed, 
    fold = FOLDS,
    model = MODELS, 
    model_version = wildcards.model_version, 
    outcome_type = wildcards.outcome_type, 
    saav_cutoff = wildcards.saav_cutoff, 
    hla_cutoff = wildcards.hla_cutoff, 
    run_id = RUN_ID
  )
  #print(split_files)
  return(split_files)

rule write_fnames2file:
	input: CV_checkpoint_aggregate
	output: temp('{run_id}/results/{protein}/{outcome_type}/{model_version}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/filenames.txt')
	run:
		with open(output[0], 'w') as f: 
			f.write('\n'.join(input))
  

rule CV_aggregate: 
	input: rules.write_fnames2file.output
	params: 
	  script_dir = config['script_dir'],
	  output_pred = '{run_id}/results/{protein}/{outcome_type}/{model_version}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/CV/pred.rds'
	output: '{run_id}/results/{protein}/{outcome_type}/{model_version}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/CV/adapt.csv'
	run:
	  exec = "Rscript {params.script_dir}/eval/compute_adaptation_score.R -i {input} --output_adapt {output} --output_pred {params.output_pred} --outcome_type {wildcards.outcome_type} --rand"
	  shell(exec)
