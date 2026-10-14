#!/bin/bash
#SBATCH -D /work/thompsok

# ----------------------------------------------------------------------
# slurm arguments
# ----------------------------------------------------------------------

#SBATCH -J relate_NLCD_CONUS_to_BBS_routes
#SBATCH -t 0-10:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --array=1-49

# ----------------------------------------------------------------------
# setup job output/error reports
# ----------------------------------------------------------------------

#SBATCH -o /work/%u/Beta_Div_Birds_Env/Land_Cover/%x-%j-%a_log.txt


# ----------------------------------------------------------------------
# load required modules
# ----------------------------------------------------------------------
module load foss/2022b R/4.2.2

# ----------------------------------------------------------------------
# execute task
# ----------------------------------------------------------------------

# run code
code=/gpfs1/data/iupdate/Rcode2/Beta_div/05_Relating_NLCD_and_CONUS_to_BBS_routes.R

Rscript "$code" 