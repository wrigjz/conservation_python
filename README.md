# conservation_python
#########################################################################################################################
## Jon Wright, IBMS, Academia Sinica, Taipei, 11529, Taiwan
## These files are licensed under the GLP ver 3, essentially you have the right
## to copy, modify and distribute this script but all modifications must be offered
## back to the original authors
#########################################################################################################################
#
This is an implimentation of the consurf system
This code has been tested on Centos 7.x and Ubuntu 22.04

Requirements:
python3 with BioPython (such as anaconda), uniref90 db, HMMER, CD-Hit, Mapsci,
Protest and rate4site

It takes a single chain pdb file (often provided from critires) or a single chain fasta file
It uses the 300 top acceptable homologs

Each CS run must be done in a separate directory with a fasta or pdb file c$as the input

To run it if you have a pdb coordinate file as:
   $consurf_scripts/consurf_home.sh file.pdb

To run it if you have a fasta file as:
   $consurf_scripts/consurf_home.sh file.fasta

If you want to use the SEQRES records from the pdb file instead of the ATOM records then use:
   $consurf_scripts/consurf_seqres.sh file.pdb

Note the numbering will run from 1->X however
