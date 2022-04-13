- NGS reads were mapped with a customized version of MinVar
    Huber, Michael, et al. "MinVar: a rapid and versatile tool for HIV-1 drug resistance genotyping by deep sequencing." Journal of virological methods 240 (2017): 7-13.
- Used implementation: Biocontainer via Docker

- Customizations: 
    o adjusted to full genome reference sequences
    o replaced question marks in ref seqs by hand
    o removed sampling of reads
    o stopped pipeline after mapping stage
    o mapping call adjusted to consider paired end data
- Subtype-specific reference is chosen by sequence similarity according to the
  sequences provided in 'subtype_references.fasta'
  and then iteratively turned into the sample consensus


Output Files:
-----------
mapped_reads/
- bam: mapping result
- fasta: consensus sample that was used for mapping
--------
subtype_references.fasta
- subtype-specific reference sequences (used refs: predominantly C, some Gs, few Bs)
