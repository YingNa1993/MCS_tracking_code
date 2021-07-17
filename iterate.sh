#!/bin/bash
#SBATCH --job-name=0529          # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --ntasks=1             # Specify max. number of tasks to be invoked
#SBATCH --cpus-per-task=1      # Specify number of CPUs per task
#SBATCH --time=08:00:00        # Set a limit on the total run time
#SBATCH --account=bm0982       # Charge resources on this project account
#SBATCH --output=track0529.o%j    # File name for standard output
#SBATCH --error=track0529.o%j     # File name for standard error output
#shopt -s expand_aliases
#source ~/.bashrc

set -ex

   PFAD=`pwd`


for year in $( seq 2001 2007)
do 

FILENAMESTRING='gpmmm_cn_'${year}'04'
cd $PFAD
${PFAD}/irt_objects_release.x ${FILENAMESTRING}
${PFAD}/irt_tracks_release.x ${FILENAMESTRING}
${PFAD}/irt_tracklinks_release.x ${FILENAMESTRING}
cd /disk1/nay/PE1/4_iterative_raincell_tracking/3_3_olrpermonth/gpm_1/output/
sort -n -k2 irt_tracks_nohead_output_${FILENAMESTRING}.txt > irt_tracks_sorted_${FILENAMESTRING}.txt
${PFAD}/irt_trackmask_release.x ${FILENAMESTRING}
#cdo -r -f nc copy irt_objects_mask_${FILENAMESTRING}.srv irt_objects_mask_${FILENAMESTRING}.nc
cdo -r -f nc copy irt_tracks_mask_ccs_${FILENAMESTRING}.srv irt_tracks_mask_ccs_${FILENAMESTRING}.nc
cdo -r -f nc copy irt_tracks_mask_${FILENAMESTRING}.srv irt_tracks_mask_${FILENAMESTRING}.nc
rm irt_objects_mask_${FILENAMESTRING}.srv
rm irt_tracks_mask_ccs_${FILENAMESTRING}.srv
rm irt_tracks_mask_${FILENAMESTRING}.srv

FILENAMESTRING='gpmmm_cn_'${year}'06'
cd $PFAD
${PFAD}/irt_objects_release.x ${FILENAMESTRING}
${PFAD}/irt_tracks_release.x ${FILENAMESTRING}
${PFAD}/irt_tracklinks_release.x ${FILENAMESTRING}
cd /disk1/nay/PE1/4_iterative_raincell_tracking/3_3_olrpermonth/gpm_1/output/
sort -n -k2 irt_tracks_nohead_output_${FILENAMESTRING}.txt > irt_tracks_sorted_${FILENAMESTRING}.txt
${PFAD}/irt_trackmask_release.x ${FILENAMESTRING}
#cdo -r -f nc copy irt_objects_mask_${FILENAMESTRING}.srv irt_objects_mask_${FILENAMESTRING}.nc
cdo -r -f nc copy irt_tracks_mask_ccs_${FILENAMESTRING}.srv irt_tracks_mask_ccs_${FILENAMESTRING}.nc
cdo -r -f nc copy irt_tracks_mask_${FILENAMESTRING}.srv irt_tracks_mask_${FILENAMESTRING}.nc
rm irt_objects_mask_${FILENAMESTRING}.srv
rm irt_tracks_mask_ccs_${FILENAMESTRING}.srv
rm irt_tracks_mask_${FILENAMESTRING}.srv

done 

exit
