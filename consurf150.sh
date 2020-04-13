#!/bin/bash

# A simple script to take the 150 consuef values instead of the 300s
#
source /home/programs/anaconda/linux-5.3.6/init.sh
dbdir=/scratch/consurf_db
blastdir=/home/programs/ncbi-blast/ncbi-blast-2.9.0_linux
hmmerdir=/home/programs/hmmer-3.1b2/linux
cdhitdir=/home/programs/cd-hit-v4.6.8-2017-0621
mafftdir=/home/programs/mafft-7.294/
rate4sitedir=/home/programs/rate4site-3.0.0/src/rate4site/
prottestdir=/home/programs/prottest-3.4.2
scripts=consurf_scripts

# Use mapsci to produce an alignment
echo "Aligning the final sequences"
$mafftdir/bin/mafft-linsi --quiet --localpair --maxiterate 1000 --reorder --clustalout --namelength 30 ./accepted150.fasta >| ./postalignment150.aln

# Get the best protein matrix
echo "Running Prottest"
java -jar $prottestdir/prottest-3.4.2.jar -i postalignment150.aln -JTT -LG -MtREV -Dayhoff -WAG -CpREV -S 1 -threads 2 >| prottest150.out
best_model=`grep 'Best model according to BIC:' prottest.out | awk  '{print $6}'`
if [ "$best_model" == "JTT" ] ; then
    rate_model="-Mj"
elif [ "$best_model" == "LG" ] ; then
    rate_model="-Ml"
elif [ "$best_model" == "MtREV" ] ; then
    rate_model="-Mr"
elif [ "$best_model" == "Dayhoff" ] ; then
    rate_model="-Md"
elif [ "$best_model" == "WAG" ] ; then
    rate_model="-Mw"
elif [ "$best_model" == "CpREV" ] ; then
    rate_model="-MC"
else 
    rate_model="-Mj"
fi

# Run the rate4site to ge the sconsurf scores
echo "Running rate4site and grading the scores"
$rate4sitedir/rate4site_doublerep -ib -a 'PDB_ATOM' -s ./postalignment150.aln -zn $rate_model -bn -l ./r4s150.log -o ./r4s150.res  -x r4s150.txt >| r4s150.out

# Turn those scores into grades
PYTHONPATH=. python3 ../$scripts/r4s_to_grades.py r4s150.res initial150.grades
