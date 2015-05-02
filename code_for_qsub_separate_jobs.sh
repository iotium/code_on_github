#!/bin/bash

# look for .out and .error files and delete if they exist

# open a file to write the submission script

# test_no=$1

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
matlab -nodesktop -nosplash -nojvm -r "bubble_growth('$2','$3','$4','$5'); quit"
EOF

qsub "bubble_growth_run$1.submit"
