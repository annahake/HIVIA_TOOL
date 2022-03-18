# Copyright 2019 Anna Hake.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms. 

  
rule create_saav_list:
  input: '{run_id}/results/{protein}/{outcome_type}/{model_version}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/CV/adapt.csv'
  output: '{run_id}/models/{protein}/{outcome_type}/final_models/{model_version}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/sites.csv'
  run:
    path = wildcards.run_id + '/'+'models/' +  wildcards.protein + '/' + wildcards.outcome_type +'/final_models/' + wildcards.model_version + '/*'
    SAAVs_paths=glob.glob(path)
    SAAVs = [output.split('/')[-1] for output in SAAVs_paths]
    pd.DataFrame(SAAVs,columns=['sites']).to_csv(output[0], index = None)
  
  
rule final_model_convert_train_hla2binary:
  input:
    rules.select_variables.output
  params:
    script_dir=config['script_dir'],
    hla = lambda wildcards: config['hla'][wildcards.model_version]
  output: 
    temp('{run_id}/models/{protein}/{outcome_type}/final_models/{model_version}/{saav}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/df_full_bin.rds')
  log:
    'logs/make_binary_matrix/{run_id}/{protein}/{outcome_type}/{model_version}/{saav}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/log.log'
  run:
    hla = ",".join(params.hla)
    exec = 'Rscript {params.script_dir}/modeling/convert_train_hla2binary.R'
    exec += ' -i {input} -c {wildcards.hla_cutoff} -o {output}'
    exec += ' --hla ' + hla
    exec += ' -l {log}'
    shell(exec)

    
rule final_model_fit:
  input: 
    df = rules.final_model_convert_train_hla2binary.output,
    tree = ancient(rules.run_raxml.output), 
  params: 
    # script directory
    script_dir = config['script_dir'],
    # fit params
    method = 'brms', 
    prior = lambda wildcards: config['modeling']['prior'][wildcards.model]['name'],
    par_ratio = lambda wildcards: config['modeling']['prior'][wildcards.model]['par_ratio'],
    adapt_delta=config["modeling"]["adapt_delta"],
    fixed = lambda wildcards: config['modeling']['fixed'][wildcards.model_version], 
    random = lambda wildcards: config['modeling']['random'], 
    outcome = lambda wildcards: config['modeling']['outcome'],
    coeffs = '{run_id}/models/{protein}/{outcome_type}/final_models/{model_version}/{saav}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/{model}/fit_brms_phylo-raxml_seqtype-nt_glmm_coeffs.csv',
    ref_map = rules.get_ref_map.output.ref_map
  output:
    '{run_id}/models/{protein}/{outcome_type}/final_models/{model_version}/{saav}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/{model}/fit_brms_phylo-raxml_seqtype-nt_glmm.rds'
  log: 
    'logs/fit_brms/{run_id}/{protein}/{outcome_type}/final_models/{model_version}/{saav}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/{model}/log.txt'
  run: 
    import os
    my_env=dict(os.environ)
    ENV_VARS=['HOME', 'TMPDIR','TMP','HOSTNAME','PATH','SGE_STDERR_PATH', 'SGE_O_WORKDIR' ,'SGE_O_HOST','SGE_O_HOME', 'PYTHONPATH']
    for key, value in my_env.items():
      print(key, value)
    
    fixed = ",".join(params.fixed)
    #random = ",".join(params.random)
    #outcome = ",".join(params.outcome)
    
    # build fit
    exec = "Rscript {params.script_dir}/modeling/fit_glmm.R"
    exec += " -i {input.df} -t {input.tree} --method {params.method}"
    exec += " --model {wildcards.model} -p {params.prior} -s {wildcards.seed}"
    exec += " --adapt_delta {params.adapt_delta} -o {output}"
    exec += " -f " + fixed
    exec += " -r {params.random}"
    exec += " -y {params.outcome}"
    exec += " --coeffs {params.coeffs}"
    exec += " --outcome_type {wildcards.outcome_type}"
    exec += " --model_version {wildcards.model_version}"
    exec += " --protein {wildcards.protein}"
    exec += " --saav {wildcards.saav}"
    exec += " --ref_map {params.ref_map}"
    exec += " --par_ratio {params.par_ratio}"
    shell(exec)

def checkpoint_aggregate(wildcards):
  checkpoint_output = checkpoints.create_datasets.get(**wildcards).output[0]
  split_files = expand('{run_id}/models/{protein}/{outcome_type}/final_models/{model_version}/{saav}/saav_freq_{saav_cutoff}/hla_freq_{hla_cutoff}/cd4_cutoff{cd4}/seed_{seed}/{model}/fit_brms_phylo-raxml_seqtype-nt_glmm.rds', 
    protein=wildcards.protein, 
    model=MODELS,
    saav=glob_wildcards(os.path.join(checkpoint_output, '{saav}/full_data/df.rds')).saav,
    cd4 = wildcards.cd4, 
    model_version = wildcards.model_version, 
    outcome_type = wildcards.outcome_type, 
    saav_cutoff = wildcards.saav_cutoff, 
    hla_cutoff = wildcards.hla_cutoff, 
    seed = wildcards.seed, 
    run_id = wildcards.run_id
  )
  print(split_files)
  return split_files



rule aggregate: 
	input: checkpoint_aggregate
	output: temp(touch("{run_id}/results/{protein}/cd4_cutoff{cd4}_{seed}_{model}_{outcome_type}_{model_version}_{saav_cutoff}_{hla_cutoff}_checkpoint_done_par_ratio.txt"))
