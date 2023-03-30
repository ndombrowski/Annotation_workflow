#!/bin/sh
#SBATCH --partition=LMEMX
##SBATCH --nodelist=no88
#SBATCH --nodes=1	# require 1 nodes
#SBATCH --ntasks-per-node=30  # (by default, "ntasks"="cpus")
#SBATCH --mem-per-cpu=1GB
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

# Modules:
module load parallel/19.0122
module load interproscan/5.48.83

#prepare GNU parallel
source `which env_parallel.bash`

# Executable commands :
parallel -j3 'i={}; interproscan.sh --cpu 10 -i $i -d IPRscan/ -T IPRscan/temp --iprlookup --goterms -appl TIGRFAM,SFLD,SUPERFAMILY,Gene3D,Hamap,Coils,ProSiteProfiles,SMART,CDD,PRINTS,ProSitePatterns,Pfam,MobiDBLite,PIRSF,TMHMM,SignalP_GRAM_NEGATIVE,SignalP_GRAM_POSITIVE' ::: split_faa/File*.faa

