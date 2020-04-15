#!/bin/bash

# A simple script to take the 150sp consuef values instead of the 300s
#
source /home/programs/anaconda/linux-5.3.6/init.sh
dbdir=/scratch/consurf_db
blastdir=/home/programs/ncbi-blast/ncbi-blast-2.9.0_linux
hmmerdir=/home/programs/hmmer-3.1b2/linux
cdhitdir=/home/programs/cd-hit-v4.6.8-2017-0621
mafftdir=/home/programs/mafft-7.294/
rate4sitedir=/home/programs/rate4site-3.0.0/src/rate4site/
prottestdir=/home/programs/prottest-3.4.2
consscripts=../consurf_scripts
critiscripts=../critires_scripts

# Use mapsci to produce an alignment
echo "Aligning the final sequences"
$mafftdir/bin/mafft-linsi --quiet --localpair --maxiterate 1000 --reorder --clustalout --namelength 30 ./accepted150sp.fasta >| ./postalignment150sp.aln

# Get the best protein matrix
echo "Running Prottest"
java -jar $prottestdir/prottest-3.4.2.jar -i postalignment150sp.aln -JTT -LG -MtREV -Dayhoff -WAG -CpREV -S 1 -threads 2 >| prottest150sp.out
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
$rate4sitedir/rate4site_doublerep -ib -a 'PDB_ATOM' -s ./postalignment150sp.aln -zn $rate_model -bn -l ./r4s150sp.log -o ./r4s150sp.res  -x r4s150sp.txt >| r4s150sp.out

# Turn those scores into grades
PYTHONPATH=. python3 $consscripts/r4s_to_grades.py r4s150sp.res initial150sp.grades

# Pull all the data togeather
python $critiscripts/get_consurf_home.py initial150sp.grades wild_consurf150sp.txt
ln -s wild_consurf150sp.txt wild_consurf.txt
python3 $critiscripts/assemble_data.py >| assemble.txt

# Get the stable/unstabla and results in the PDB numbering scheme
/bin/rm -rf results_ambnum.txt results.txt
$critiscripts/set_numbers.sh
python3 $critiscripts/find_stable_unstable.py >| results_ambnum.txt
PYTHONPATH=. python3 $critiscripts/print_results.py >| results.txt

rm wild_consurf.txt
mv assemble.txt assemble150sp.txt
mv results_ambnum.txt results_ambnum150sp.txt
mv results.txt results150sp.txt
