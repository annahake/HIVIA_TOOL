# Copyright 2019 Anna Hake.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms. 

# define local rules which are not run on cluster
localrules: clean_coreceptor_data, clean_clinical_data, collect_hla_alleles, clean_hla_table, fetch_sequences, clean_subtypes, select_seq, clean_alignment, remove_refseq, convert_fasta_to_phylip,get_aa_ref_seq, codon_align

# PREPROCESSING

rule clean_coreceptor_data:
  input: rules.g2p_coreceptor.output
  output: "data/processed/coreceptor.csv"
  log: "logs/coreceptor/clean_coreceptor_data.log"
  params: 
    script_dir = config['script_dir']
  shell:
    "Rscript {params.script_dir}/prep/coreceptor_data/clean_coreceptor_data.R -i {input} -o {output} -l {log}"

## subtype sequences with the comet tool
rule subtype_sequences:
	input:
		sequences = 'data/interim/env'+str(config['fetch_seq']['consensus_cutoff'])+config['fetch_seq']['output_fname_suffix'] +config['fetch_seq']['output_fname_extension']
	output: 
		subtypes = 'data/interim/env'+str(config['fetch_seq']['consensus_cutoff'])+config['fetch_seq']['output_fname_suffix']+config['subtype_seq']['output_fname_suffix'] +config['subtype_seq']['output_fname_extension']
	params:
		api=config['api']['comet'],
		script_dir = config['script_dir']
	shell:
		"Rscript {params.script_dir}/prep/sequence_data/subtype_sequences.R -i {input.sequences} -o {output.subtypes} -a {params.api}"

## clean the subtype information from the comet tool
rule clean_subtypes:
  input:
    subtypes = rules.subtype_sequences.output.subtypes
  output:
    cleaned_subtypes = 'data/interim/env'+str(config['fetch_seq']['consensus_cutoff'])+config['fetch_seq']['output_fname_suffix']+config['subtype_seq']['output_fname_suffix'] +'_cleaned.csv'
  params:
    script_dir = config['script_dir']
  shell:
    'Rscript {params.script_dir}/prep/sequence_data/clean_subtype_table.R -s {input.subtypes}'

#TODO check if nmer tbl still needed
rule select_seq:
  input:
    sequences = rules.select_seq_clin.output.sequences,
    cleaned_subtypes = rules.clean_subtypes.output.cleaned_subtypes
  output:
    selected_sequences = 'data/interim/{protein}'+str(config['fetch_seq']['consensus_cutoff'])+config['fetch_seq']['output_fname_suffix'] + config['select_seq']['output_fname_suffix']+config['select_seq']['output_fname_extension'],
    nmer_tbl='tables/{protein}_nmer_distribution.tex'
  params:
    script_dir = config['script_dir']
  log: 'logs/excluded/{protein}/excluded.csv'
  shell:
    "Rscript {params.script_dir}/prep/sequence_data/select_sequences.R -f {input.sequences} -s {input.cleaned_subtypes} -o {output.selected_sequences} -l {log} -t {output.nmer_tbl}"
#invoked?
rule select_seq_leftout:
  input:
    sequences = rules.select_seq_clin.output.leftout,
    cleaned_subtypes = rules.clean_subtypes.output.cleaned_subtypes
  output:
    leftout_selected_seqs = 'data/interim/{protein}'+str(config['fetch_seq']['consensus_cutoff'])+config['fetch_seq']['output_fname_suffix'] + config['select_seq']['output_fname_suffix']+'_leftout_clinical'+config['select_seq']['output_fname_extension'],
    #nmer_tbl='tables/{protein}_nmer_distribution.tex'
  params:
    script_dir = config['script_dir']
  #log: 'logs/excluded/{protein}/excluded.csv'
  shell:
    "Rscript {params.script_dir}/prep/sequence_data/select_sequences.R -f {input.sequences} -s {input.cleaned_subtypes} -o {output.selected_sequences}"
### SEPARATE HERE IN SPECIFIC AND GENERAL

# get the subtype reference sequence from LANL NewAlign Tool

rule get_reference_sequence:
  # define the output files
  output:
    # fasta file with the reference sequence of interest
    refseq_file = temp("{run_id}/data/interim/{protein}_"+config['ref_seq']['query']+"_"+config['ref_seq']['ref_seq_fname']),
    # fasta file with the corresponding full alignment
    alignment_file = "{run_id}/data/interim/{protein}_"+config['ref_seq']['align_fname']
  # set parameters that are fixed and not need to be generated beforehand
  params:
    alignment_type = config['ref_seq']['alignment_type'],
    suborganism = config['ref_seq']['suborganism'],
    subtype = config['ref_seq']['subtype'],
    region_definition = config['ref_seq']['region_definition'],
    region = lambda wildcards: config['ref_seq']['region'][wildcards.protein],
    basetype = config['ref_seq']['basetype'],
    year = config['ref_seq']['year'],
    format = config['ref_seq']['format'],
    query = config['ref_seq']['query'],
    api_source = config["api"]['lanl_align'],
    script_dir = config['script_dir']
  shell:
    "Rscript {params.script_dir}/prep/sequence_data/get_reference_sequence.R -a {output.alignment_file} -o {output.refseq_file} -t {params.alignment_type} -s {params.suborganism} --subtype {params.subtype} -d {params.region_definition} -r {params.region} -b {params.basetype} -y {params.year} -f {params.format} -q {params.query} --api_source {params.api_source}"
    

## align all sequences
rule align_seq:
  input:
    selected_sequences = rules.select_seq.output.selected_sequences,
    reference_sequence = rules.get_reference_sequence.output.refseq_file
  output:
    alignment = temp('{run_id}/data/interim/{protein}'+str(config['fetch_seq']['consensus_cutoff'])+config['fetch_seq']['output_fname_suffix']+config['subtype_seq']['output_fname_suffix'] + config['select_seq']['output_fname_suffix']+config['align_seq']['output_fname_suffix'] + config['align_seq']['output_fname_extension'])
  params:
    script_dir = config['script_dir']
  shell:
    "mafft --add {input.selected_sequences} --quiet {input.reference_sequence}  > {output.alignment}"

#rule align_leftout_seq:
#  input:
#    selected_sequences = rules.select_seq_leftout.output.leftout_selected_seqs,
#    reference_sequence = rules.get_reference_sequence.output.refseq_file
#  output:
#    alignment = temp('data/interim/{protein}'+str(config['fetch_seq']['consensus_cutoff'])+config['fetch_seq']['output_fname_suffix']+config['subtype_seq']['output_fname_suffix'] + config['select_seq']['output_fname_suffix']+config['align_seq']['output_fname_suffix'] + '_leftout_clinical'+ config['align_seq']['output_fname_extension'])
#  params:
#    script_dir = config['script_dir']
#  shell:
#    "mafft --add {input.selected_sequences} --quiet {input.reference_sequence}  > {output.alignment}"
    
## get HXB2 reference protein sequence
rule get_aa_ref_seq:
  params:
    id = lambda wildcards: config['aa_ref_seq']['id']['protein'][wildcards.protein],
    db = config['aa_ref_seq']['database'],
    protein = lambda wildcards: wildcards.protein,
    script_dir = config['script_dir']
  output:
    ref_seq_aa = temp('{run_id}/data/interim/{protein}_' + config['aa_ref_seq']['database'] + '.fasta')
  shell:
    "Rscript {params.script_dir}/prep/sequence_data/get_reference_sequence2.R -i {params.id} -n {params.protein} -d {params.db} -o {output}"
    
## codon align all sequences
rule codon_align:
  input:
    dna_alignment = rules.align_seq.output.alignment
  params:
    frame = config["codon_align"]["frame"],
    compensating_num = config["codon_align"]["compensating_num"],
    input_type = config["codon_align"]["input_type"],
    output_format = config["codon_align"]["output_format"],
    output_seqs_type = config["codon_align"]["output_seqs_type"],
    protein = lambda wildcards: wildcards.protein,
    output_dir = '{run_id}/data/interim/',
    api_source = config['api']['codon_align'],
    script_dir = config['script_dir']
  output:
    codon_align_nt = "{run_id}/data/interim/{protein}_nt_codonalign"+"."+ config["codon_align"]["output_format"], 
    codon_align_aa = "{run_id}/data/interim/{protein}_aa_codonalign"+"."+ config["codon_align"]["output_format"]
  shell:
    "Rscript {params.script_dir}/prep/sequence_data/get_codon_alignment.R -i {input.dna_alignment} -f {params.frame} -c {params.compensating_num} -t {params.input_type} --output_format {params.output_format} --output_seqs_type {params.output_seqs_type} -n {params.protein} --output_dir {params.output_dir} --api_source {params.api_source}"

## clean codon alignment
rule clean_alignment:
  input:
    codon_align_aa = rules.codon_align.output.codon_align_aa
  output:
    ref_seq = "{run_id}/data/interim/{protein}_aa_codonalign_"+config['ref_seq']['query']+"_"+config['ref_seq']['ref_seq_fname']	,
    alignment = temp("{run_id}/data/interim/{protein}_aa_codonalign_cleaned.csv"),
    seq_freq_mat = "{run_id}/data/interim/{protein}_aa_codonalign_seqfreqmat.csv",
    pos_freq_mat = "{run_id}/data/interim/{protein}_aa_codonalign_posfreqmat.csv"
  log: 'logs/codonalign/{run_id}/{protein}/clean_alignment.log'
  params:
    script_dir = config['script_dir']
  run: 
    exec = "Rscript {params.script_dir}/prep/sequence_data/clean_codon_align.R "
    exec+= "-i {input.codon_align_aa} "
    exec+= "-r {output.ref_seq} "
    exec+= "-a {output.alignment} "
    exec+= "-l {log} "
    exec+= "-s {output.seq_freq_mat} "
    exec+= "-p {output.pos_freq_mat} "
    shell(exec)
    

## create reference map between subtype numbering and HXB2
rule get_ref_map:
  input: 
    codon_align_aa_ref_seq = rules.clean_alignment.output.ref_seq,
    ref_seq_aa = rules.get_aa_ref_seq.output.ref_seq_aa
  output:
    align = temp("{run_id}/data/interim/{protein}_aa_codonalign_withrefseq"+"."+ config["codon_align"]["output_format"]),
    ref_map = '{run_id}/data/interim/{protein}_' + config['aa_ref_seq']['database'] + '.fasta.map'
  params:
    script_dir = config['script_dir']
  shell:
    "mafft --add {input.ref_seq_aa} --quiet --mapout {input.codon_align_aa_ref_seq} >  {output.align}"
    
## remove reference sequence from alignment
rule remove_refseq:
  input: rules.codon_align.output.codon_align_nt
  output: '{run_id}/data/interim/{protein}/codonalign/nt/codon_align_without_refseq.fasta'
  params:
    refseq_name=config['ref_seq']['query'],
    script_dir = config['script_dir']
  #log: 'logs/{protein}/remove_refseq.log'
  shell: 
    'Rscript {params.script_dir}/prep/sequence_data/remove_refseq.R -i {input} -o {output} -r {params.refseq_name}'
  # Alternatively delete the fasta entry with the refseq_name using python, not checked
  #run:
#    from Bio import SeqIO
#    with open(output[0],'w') as f:
#      for seq in SeqIO.parse(input, 'fasta'):
#        if seq.id != params.refseq_name:
#          SeqIO.write(seq, f, "fasta")
  
  
rule make_SAAV_table:
  input:
	# alignment
    alignment = rules.clean_alignment.output.alignment,
    ref_map = rules.get_ref_map.output.ref_map
  params:
    script_dir = config['script_dir']
  output:
    '{run_id}/data/processed/{protein}_{outcome_type}_{saav_cutoff}_SAAV_table.csv'
  shell:
    "Rscript {params.script_dir}/prep/sequence_data/make_SAAV_table.R -a {input.alignment} -r {input.ref_map} -d {wildcards.saav_cutoff} -o {output} --outcome_type {wildcards.outcome_type}" 

# TODO coreceptor and clinical data should be deleted as it is not necessary
checkpoint create_datasets:
  input:
    hla = rules.clean_hla_table.output,
    SAAV = rules.make_SAAV_table.output,
    coreceptor = rules.clean_coreceptor_data.output, 
    clinical = rules.clean_clinical_data.output
  params:
    dir = '{run_id}/data/processed/datasets/saav_freq_{saav_cutoff}/cd4_cutoff{cd4}/{protein}/{outcome_type}/{model_version}/', 
    script_dir = config["script_dir"]
  output: 
    directory('{run_id}/data/processed/datasets/saav_freq_{saav_cutoff}/cd4_cutoff{cd4}/{protein}/{outcome_type}/{model_version}/')
  shell: 
     "Rscript {params.script_dir}/prep/create_datasets.R --hla {input.hla} -s {input.SAAV} --coreceptor {input.coreceptor} --clinical {input.clinical} -d {params.dir} --cd4_cutoff {wildcards.cd4}"

rule create_random_dataset: 
  input:
    hla = rules.clean_hla_table.output,
    clinical = rules.clean_clinical_data.output
  params: 
    script_dir = config["script_dir"], 
    rand_nr = config['modeling']['rand_nr']
  output: 
    "{run_id}/data/processed/rand/seed_{seed}/cd4_cutoff{cd4}/df.rds"
  shell: 
    "Rscript {params.script_dir}/prep/create_rand_dataset.R --hla {input.hla} -c {input.clinical} -s {wildcards.seed} --cd4_cutoff {wildcards.cd4} -o {output} --rand_nr {params.rand_nr}"
      

rule select_variables:
  input: '{run_id}/data/processed/datasets/saav_freq_{saav_cutoff}/cd4_cutoff{cd4}/{protein}/{outcome_type}/{model_version}/{saav}/full_data/df.rds',
  params: 
    fixed = lambda wildcards: config['modeling']['fixed'][wildcards.model_version], 
    random = config['modeling']['random'], 
    outcome = config['modeling']['outcome'],
    hla = lambda wildcards: config['hla'][wildcards.model_version],
    script_dir=config['script_dir'],
    id_column = 'saav_id'
  output: 
    '{run_id}/data/processed/datasets/saav_freq_{saav_cutoff}/cd4_cutoff{cd4}/{protein}/{outcome_type}/{model_version}/{saav}/full_data/df_full.rds'
  #log:
  # 'logs/select_variables/cd4_cutoff{cd4}/{protein}/{outcome_type}/{model_version}/{saav}/full_data/log.txt'
  run:
    fixed = ",".join(params.fixed)
    hla = ",".join(params.hla)
    
    exec = 'Rscript {params.script_dir}/modeling/select_variables.R'
    exec += ' -i {input} -o {output}'
    exec += ' -f ' + fixed
    exec += ' -r {params.random}'
    exec += ' -y {params.outcome}'
    exec += ' --id  {params.id_column}'
    exec += ' --hla ' + hla
    
    shell(exec)
    

rule run_raxml: 
  input: rules.remove_refseq.output
  output: '{run_id}/data/processed/{protein}/phylo/raxml/RAxML_bestTree.{run_id}'
  #output: '{run_id}/data/processed/{protein}/phylo/raxml/RAxML_bestTree.{protein}'
  params: 
    model=config["raxml"]["model"],
    run_id = RUN_ID, 
    dir = '{run_id}/data/processed/{protein}/phylo/raxml/',
    seed = config["seed"]
  log: 'logs/{run_id}/{protein}/phylo/raxml.log' 
  benchmark:
    "benchmarks/{run_id}/{protein}/raxml.txt"
  shell:
    'raxmlHPC -s {input} -w $(pwd)/{params.dir} -n {params.run_id} -m {params.model} -p {params.seed} --silent'
    
