#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=storeBestAsFits
#SBATCH --mail-type=END
#SBATCH --output=slurm_out/slurm_%x.out
#SBATCH --mail-user=ntf229@nyu.edu

module purge

cd /home/ntf229/containers

singularity exec --overlay overlay-15GB-500K.ext3:ro \
	    /scratch/work/public/singularity/cuda11.0-cudnn8-devel-ubuntu18.04.sif \
/bin/bash -c "source /ext3/env.sh; 
python /home/ntf229/nihao2/bin/storeBestAsFits.py \
--ageSmooth=False --SF=False --tauClear=2.5 \
--numPhotons=1e9 --pixels=2000 \
--dustFraction=0.1 --maxTemp=16000 --SSP=FSPS_Chabrier"



