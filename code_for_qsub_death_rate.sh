#!/bin/bash

#  code_for_qsub_C.sh
#  
#
#  Created by Jonah Zimmerman on 12/8/14.
#

# look for .out and .error files and delet if they exist

# open a file to write the submission script

if [ -e qsubrun.out ]; then
rm qsubrun.out
fi

if [ -e qsubrun.error ]; then
rm qsubrun.error
fi

cat > bubble_growth_death_run.submit << EOF
#!/bin/bash
#$ -o qsubrun.out
#$ -e qsubrun.error
#$ -cwd
#$ -S /bin/bash

module load matlab
matlab -nodesktop < eval_death_rate.m
EOF

qsub bubble_growth_death_run.submit