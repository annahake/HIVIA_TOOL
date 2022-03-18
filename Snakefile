# Copyright 2019 Anna Hake.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms. 

###
# import libraries
###

import glob
import os.path
import numpy
from pathlib import Path
import pandas as pd


###
# email report
###
onsuccess:
    body_text = "Nothing to report\n{log}"
    if config['notify']:
        shell('mail -s "[Snakemake] SUCCESS - HIVIA" {} <<< "{}"'.format(config['notify_email'], body_text))

onerror:
    if config['notify']:
      shell('mail -s "[Snakemake] ERROR - HIVIA" {} < {{log}}'.format(config['notify_email']))

###
# define variables
###
PROTEINS = config["proteins"]
FOLDS = numpy.arange(config["modeling"]["kfold"]) + 1
SEEDS = config['modeling']["seed"]
MODELS = config["modeling"]["model"]
CD4 = config['cd4_cutoff']
RANDS = config['modeling']['rand']
MODEL_VERSIONS = config["model_versions"]
OUTCOME_TYPE = config['outcome_type']
SAAV_CUTOFF = config["saav_freq_cutoff"]
HLA_CUTOFF = config["hla_freq_cutoff"]
RUN_ID = config['run_id']
wildcard_constraints:
  protein='('+'|'.join(PROTEINS)+')',
  model='('+'|'.join(MODELS)+')',
  rand= '('+'|'.join(RANDS)+')',
  seed ='('+'|'.join(SEEDS)+')',
  model_version = '('+'|'.join(MODEL_VERSIONS)+')', 
  saav_cutoff = '('+'|'.join(SAAV_CUTOFF)+')', 
  hla_cutoff = '('+'|'.join(HLA_CUTOFF)+')'


###
# functions
###  

  
  
###
# load rules (before target rules s.t. rules object can be used for dependencies. 
###
if config['include']['prep']:
  include: "rules/prep.smk"
if config['include']['final_models']:
  include:'rules/final_models.smk'
if config['include']['CV_model']:
  include: "rules/CV_model.smk"



localrules: all
###
# target rules
# if more complex, create all_input variable which is build with if statements or in a function
###
rule all:
  input:
    expand(rules.run_raxml.output, protein=PROTEINS, run_id = RUN_ID),
    expand(rules.CV_aggregate.output, protein=PROTEINS, seed=SEEDS, cd4=CD4, model= MODELS, model_version=MODEL_VERSIONS, outcome_type = OUTCOME_TYPE, saav_cutoff = SAAV_CUTOFF, hla_cutoff = HLA_CUTOFF, run_id = RUN_ID) if config['include']['CV_model'] else [],
    ## build of final models with full training data
    expand(rules.aggregate.output, protein=PROTEINS, seed=SEEDS, cd4=CD4, model= MODELS, model_version=MODEL_VERSIONS, outcome_type = OUTCOME_TYPE, saav_cutoff = SAAV_CUTOFF, hla_cutoff = HLA_CUTOFF, run_id = RUN_ID) if config['include']['final_models'] else [],
    expand(rules.create_saav_list.output,  protein=PROTEINS, run_id = RUN_ID, outcome_type = OUTCOME_TYPE, model_version=MODEL_VERSIONS, cd4=CD4, saav_cutoff = SAAV_CUTOFF, hla_cutoff = HLA_CUTOFF, seed=SEEDS) if config['include']['final_models'] else [],




