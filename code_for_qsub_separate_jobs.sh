#!/bin/bash

#  code_for_qsub_separate_jobs.sh
#  
#
#  Created by Jonah Zimmerman on 12/24/14.
#

# look for .out and .error files and delete if they exist

# open a file to write the submission script

#n=$1
#C_rdot=$2
#C_nuc_rate=$3
#C_death_rate=$4
#C_erf_death_rate=$5
#save_filename=$6

if [ -e "qsubrun$1.out" ]; then
rm "qsubrun$1.out"
fi

if [ -e "bubble_growth_run$1.submit" ]; then
rm "bubble_growth_run$1.submit"
fi

if [ -e "qsubrun$1.error" ]; then
rm "qsubrun$1.error"
fi

cat > bubble_growth_run$1.submit << EOF
#!/bin/bash
#$ -o qsubrun$1.out
#$ -e qsubrun$1.error
#$ -cwd
#$ -S /bin/bash

module load matlab
matlab -nodesktop -r "bubble_growth('$2','$3','$4','$5','$6')"
EOF

qsub "bubble_growth_run$1.submit"
