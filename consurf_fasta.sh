#!/bin/bash
###################################################################################################
## Jon Wright, IBMS, Academia Sinica, Taipei, 11529, Taiwan
## These files are licensed under the GLP ver 3, essentially you have the right
## to copy, modify and distribute this script but all modifications must be offered
## back to the original authors
###################################################################################################
#
# A simple script to replciate Consurf for a provided single chain fasta file

# Usage consurf_home file.pdb

if [ "$#" -ne 1 ]; then
    echo "Please give a (and only one) fasta format file"
    exit 1
fi

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

# Extract the title line from the fasta file
longtitle=`grep ^> $1`
title="${longtitle:1}"

# Remove output from previous runs
/bin/rm -rf uniref90_list.txt prealignment.fasta postalignment.aln accepted.fasta uniref.tmp 
/bin/rm -rf frequency.aln consurf_home.grades frequency.txt
/bin/rm -rf homologues.fasta r4s_pdb.py initial.grades r4s.res prottest.out cdhit.log r4s.out

# Jackhmmer the blast database looking for homologues
echo "Jackhmmering the Uniref90 DB"
$hmmerdir/binaries/jackhmmer -E 0.0001 --domE 0.0001 --incE 0.0001 -N 1 \
        -o $1\_hmmer.out -A uniref90_list.txt $1 $dbdir/uniref90.fasta

# Remove the target entry - probably not needed but good to do anyway
grep -v "$title" uniref90_list.txt >| uniref.tmp

# Retrieve the sequences that Jackhmmer found
echo "Reformating the sequences from the Uniref90 DB"
$hmmerdir/binaries/esl-reformat fasta uniref.tmp >| homologues.fasta

# Run cd-hit to cluster everything to remove duplicates at the 95% level
echo "Clustering using cdhit and selecting the sequences"
$cdhitdir/cd-hit -i ./homologues.fasta -o ./cdhit.out -c 0.95 >| cdhit.log

echo "Rejecting some sequences"
python3 ../$scripts/select_seqs.py $1.fasta cdhit.out

# Use mapsci to produce an alignment
echo "Aligning the final sequences"
$mafftdir/bin/mafft-linsi --quiet --localpair --maxiterate 1000 --reorder --clustalout --namelength 30 \
         ./accepted.fasta >| ./postalignment.aln
$mafftdir/bin/mafft-linsi --quiet --localpair --maxiterate 1000 --reorder --namelength 30 \
         ./accepted.fasta >| frequency.aln
python3 ../$scripts/get_frequency.py frequency.aln >| frequency.txt

# Get the best protein matrix
echo "Running Prottest"
java -jar $prottestdir/prottest-3.4.2.jar -i postalignment.aln -JTT -LG -MtREV -Dayhoff -WAG \
        -CpREV -S 1 -threads 2 >| prottest.out
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
$rate4sitedir/rate4site_doublerep -ib -a "$title" -s ./postalignment.aln -zn $rate_model -bn \
     -l ./r4s.log -o ./r4s.res  -x r4s.txt >| r4s.out

# Turn those scores into grades
PYTHONPATH=. python3 ../$scripts/r4s_to_grades_fasta.py r4s.res initial.grades
paste initial.grades frequency.txt > consurf_home.grades
