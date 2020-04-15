#!/bin/bash

# A simple script to take the 300top consuef values instead of the 300s
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
$mafftdir/bin/mafft-linsi --quiet --localpair --maxiterate 1000 --reorder --clustalout --namelength 30 ./accepted300top.fasta >| ./postalignment300top.aln

# Get the best protein matrix
echo "Running Prottest"
java -jar $prottestdir/prottest-3.4.2.jar -i postalignment300top.aln -JTT -LG -MtREV -Dayhoff -WAG -CpREV -S 1 -threads 2 >| prottest300top.out
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
$rate4sitedir/rate4site_doublerep -ib -a 'PDB_ATOM' -s ./postalignment300top.aln -zn $rate_model -bn -l ./r4s300top.log -o ./r4s300top.res  -x r4s300top.txt >| r4s300top.out

# Turn those scores into grades
PYTHONPATH=. python3 $consscripts/r4s_to_grades.py r4s300top.res initial300top.grades

# Save the consurf300 data
mv assemble.txt assemble300.txt
mv initial.grades initial300.grades
mv wild_consurf.txt wild_consurf300.txt
mv results_ambnum.txt results_ambnum300.txt
mv results.txt results300.txt

# Pull all the data togeather
python $critiscripts/get_consurf_home.py initial300top.grades wild_consurf300top.txt
ln -s wild_consurf300top.txt wild_consurf.txt
python3 $critiscripts/assemble_data.py >| assemble.txt

# Get the stable/unstabla and results in the PDB numbering scheme
/bin/rm -rf results_ambnum.txt results.txt
$critiscripts/set_numbers.sh
python3 $critiscripts/find_stable_unstable.py >| results_ambnum.txt
PYTHONPATH=. python3 $critiscripts/print_results.py >| results.txt

rm wild_consurf.txt
mv assemble.txt assemble300top.txt
mv results_ambnum.txt results_ambnum300top.txt
mv results.txt results300top.txt
