# snakeSTR


Handling of hemyzigous and haploid markers was not implemented 

# example

snakemake --profile profile/slurm/  --configfile config/default.yaml --printshellcmds --latency-wait 30   --config "parameter_set"="default" "reference"="input/reference/mzib.min_150.pseudohap2.1_HiC.purged.fasta" "bam_suffix"=".mzib.hic.purged.mkdup"
