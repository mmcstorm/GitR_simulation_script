#!/bin/env bash
#SBATCH -J Sim_study
#SBATCH -N 2
#SBATCH --mem=512MB
#SBATCH --output=sim_study_%J.out
#SBATCH --error=sim_study_%J.err

module load statistical/R/4.0.2/gcc.8.3.1

Rscript ORDINAL_MainSimulationScriptSLURM.R $i $j