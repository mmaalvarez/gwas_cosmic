#!/bin/bash

## load modules, and set root dir names and partition names, depending on cluster ('executor')

hostname=`echo $HOSTNAME | cut -d"." -f1 | cut -d"-" -f1`

if [[ "$hostname" == "lsflogin" || "$hostname" == "hpc" ]]; then

	module load nextflow/24.04.2-with-plugins
	module load singularity/4.1.1
	module load graphviz/12.1.1
	export root_dir="/nas/weka.gel.zone"
	export work_dir=""$root_dir"/re_gecip/cancer_pan/malvarez/gwas_cosmic"
	export executor="lsf"
	export partition_fast_short="short"
	export partition_slow_long="medium"
	export partition_slowest_unlimited="long"
	export project_group="-P re_gecip_cancer_pan"
	export TMPDIR="/re_scratch/re_gecip/cancer_pan/malvarez"

elif [[ "$hostname" == "fsupeksvr" ]]; then

	conda activate nextflow
	export root_dir="/g"
	export work_dir=$PWD
	export executor="slurm"
	export partition_fast_short="normal_prio"
	export partition_slow_long="normal_prio_lon"
	export partition_slowest_unlimited="normal_prio_unli"
	export project_group=""

elif [[ "$hostname" == "irblogin01" ]]; then

	module load R/4.2.1-foss-2022a Nextflow Anaconda3
	export root_dir="/data/gds"
	export work_dir="$PWD"
	export partition_fast_short="irb_cpu_iclk"
	export partition_slow_long="?"
	export partition_slowest_unlimited="?"
	export project_group=""

else
	echo "ERROR: HOSTNAME is not known: '`echo $HOSTNAME | cut -d"." -f1`'"
fi


mkdir -p log/

nextflow -log $PWD/log/nextflow.log run main.nf -resume

