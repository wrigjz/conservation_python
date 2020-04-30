#!/bin/bash
###################################################################################################
## Jon Wright, IBMS, Academia Sinica, Taipei, 11529, Taiwan
## These files are licensed under the GLP ver 3, essentially you have the right
## to copy, modify and distribute this script but all modifications must be offered
## back to the original authors
###################################################################################################
#
# A simple script to replciate Consurf it looks at using the SEQRES records

# Usage consurf_home fileX.pdb
# Although here we always use original.pdb but we then have to use a pdb file from 
# THe atom numbers to be able to match with critires later

if [ "$#" -ne 1 ]; then
    echo "Please give just one pdb file - the one we will take the ATOM records from"
    exit 1
fi

# Check if we are on the cluster
if [ -z "$PBS_NODEFILE" ] ; then
    echo "Non cluster run"
    threads="1"
else
    echo "Cluster run"
    threads=`wc -l < $PBS_NODEFILE`
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
consscripts=consurf_scripts

# Remove output from previous runs
/bin/rm -rf uniref90_list.txt prealignment.fasta postalignment.aln accepted.fasta uniref.tmp 
/bin/rm -rf frequency.aln consurf_home.grades frequency.txt seqres.fasta
/bin/rm -rf homologs.fasta r4s_pdb.py initial.grades r4s.res prottest.out cdhit.log r4s.out

# Work out the chain ids
if [ ! -e "original_chain.txt" ] ; then
    echo "No orignial_chain.txt file to get the chain id from"
else
    chain=`cat original_chain.txt`
fi 

# generate the fasta file from the given pdb file
echo "Creating Fasta file"
python3 ../$consscripts/mk_fasta.py $1  >| cons.fasta
python3 ../$consscripts/mk_fasta_from_seqres.py original.pdb ${chain^}  >| seqres.fasta

# Jackhmmer the blast database looking for homologs
echo "Jackhmmering the Uniref90 DB"
$hmmerdir/binaries/jackhmmer -E 0.0001 --domE 0.0001 --incE 0.0001 -N 1 \
        -o cons_hmmer.out -A uniref90_list.txt seqres.fasta $dbdir/uniref90.fasta

# Remove the PDB_ATOM / given fasta entry - probably not needed but good to do anyway
grep -v PDB_ATOM uniref90_list.txt >| uniref.tmp

# Retrieve the sequences that Jackhmmer found
echo "Reformating the sequences from the Uniref90 DB"
$hmmerdir/binaries/esl-reformat fasta uniref.tmp >| homologs.fasta

# Run cd-hit to cluster everything to remove duplicates at the 95% level
echo "Clustering using cdhit and selecting the sequences"
$cdhitdir/cd-hit -i ./homologs.fasta -o ./cdhit.out -c 0.95 >| cdhit.log

echo "Rejecting some sequences"
python3 ../$consscripts/select_seqs.py seqres.fasta cdhit.out

# Use mapsci to produce an alignment
echo "Aligning the final sequences"
$mafftdir/bin/mafft-linsi --quiet --localpair --maxiterate 1000 --reorder --clustalout \
  --thread $threads --namelength 30 ./accepted.fasta >| ./postalignment.aln
$mafftdir/bin/mafft-linsi --quiet --localpair --maxiterate 1000 --reorder \
  --thread $threads --namelength 30 ./accepted.fasta >| frequency.aln

# Calculate the residue frequencies for homologs aligned to the inital given sequence
python3 ../$consscripts/get_frequency.py frequency.aln >| frequency.txt

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

# Run the rate4site to get the consurf scores - sometimes this fails and if so we then run
# the older version which seems to do better but has less options and is far slower
echo "Running rate4site and grading the scores"
$rate4sitedir/rate4site_doublerep -ib -a 'PDB_ATOM' -s ./postalignment.aln -zn $rate_model -bn \
       -l ./r4s.log -o ./r4s.res  -x r4s.txt >| r4s.out

# Check if rate4site ran okay, if not then we run the older version which lacks "LG" so if protest
# recommended that we need to change it to "JTT"
if [ $? -ne 0 ] ; then
    echo "R4S 3.0 failed so we'll try the older version"
    if [ "$best_model" == "LG" ] ; then
        rate_model="-Mj"
    fi
    $rate4sitedir/rate4site.old_slow -ib -a 'PDB_ATOM' -s ./postalignment.aln -zn $rate_model -bn \
       -l ./r4s.log -o ./r4s.res  -x r4s.txt >| r4s.out
fi

# Turn those scores into grades
PYTHONPATH=. python3 ../$consscripts/r4s_to_grades.py r4s.res seqres.grades N
paste seqres frequency.txt >| seqres_home.grades
