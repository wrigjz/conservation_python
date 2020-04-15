#!/bin/bash

#
# A simple script to replciate Consurf for a provided single chain pdb file called wild_min.pdb
# with a chain id of X, this is a post amber minimized pdb file
# setup anaconda environment
source /home/programs/anaconda/linux-5.3.6/init.sh
dbdir=/scratch/consurf_db
blastdir=/home/programs/ncbi-blast/ncbi-blast-2.9.0_linux
hmmerdir=/home/programs/hmmer-3.1b2/linux
cdhitdir=/home/programs/cd-hit-v4.6.8-2017-0621
mafftdir=/home/programs/mafft-7.294/
rate4sitedir=/home/programs/rate4site-3.0.0/src/rate4site/
prottestdir=/home/programs/prottest-3.4.2
scripts=consurf_scripts

# Remove output from previous runs
echo "Creating Fasta file"
/bin/rm -rf uniref90_list.txt prealignment.fasta postalignment.aln accepted.fasta
/bin/rm -rf homologues.fasta r4s_pdb.py initial.grades r4s.res prottest.out cdhit.log r4s.out

# generate the fasta file
python3 ../$scripts/mk_fasta.py wild_mini.pdb >| wild.fasta

# Jackhmmer the blast database looking for homologues
echo "Jackhmmering the Uniref90 DB"
$hmmerdir/binaries/jackhmmer -E 0.0001 --domE 0.0001 --incE 0.0001 -N 1 \
        -o wild_hmmer.out -A uniref90_list.txt wild.fasta $dbdir/uniref90.fasta

# Retrieve the sequences that Jackhmmer found
echo "Reformating the sequences from the Uniref90 DB"
$hmmerdir/binaries/esl-reformat fasta uniref90_list.txt >| homologues.fasta

# Run cd-hit to cluster everything to remove duplicates at the 95% level
echo "Clustering using cdhit and selecting the sequences"
$cdhitdir/cd-hit -i ./homologues.fasta -o ./cdhit.out -c 0.95 >| cdhit.log

echo "Rejecting some sequences"
python3 ../$scripts/select_seqs.py wild.fasta cdhit.out

# Here you get to choose either the accepted150 or accepted300 files for aligning
/bin/ln -s accepted300sp.fasta accepted.fasta

# Use mapsci to produce an alignment
echo "Aligning the final sequences"
$mafftdir/bin/mafft-linsi --quiet --localpair --maxiterate 1000 --reorder --clustalout --namelength 30 ./accepted.fasta >| ./postalignment.aln

# Get the best protein matrix
echo "Running Prottest"
java -jar $prottestdir/prottest-3.4.2.jar -i postalignment.aln -JTT -LG -MtREV -Dayhoff -WAG -CpREV -S 1 -threads 2 >| prottest.out
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
$rate4sitedir/rate4site_doublerep -ib -a 'PDB_ATOM' -s ./postalignment.aln -zn $rate_model -bn -l ./r4s.log -o ./r4s.res  -x r4s.txt >| r4s.out

# Turn those scores into grades
PYTHONPATH=. python3 ../$scripts/r4s_to_grades.py r4s.res initial.grades

# Save the consurf300sp data
mv assemble.txt assemble300sp.txt
mv initial.grades initial300sp.grades
mv wild_consurf.txt wild_consurf300sp.txt
mv results_ambnum.txt results_ambnum300sp.txt
mv results.txt results300sp.txt
