#########################################################################################################################
## Jon Wright, IBMS, Academia Sinica, Taipei, 11529, Taiwan
## These files are licensed under the GLP ver 3, essentially you have the right
## to copy, modify and distribute this script but all modifications must be offered
## back to the original authors
#########################################################################################################################
#
This is an implimentation of the consurf system
It takes a single chain pdb file (often provided from critires) or a single chain fasta file
It uses 150 evenly sampled acceptable homologs 

To run you need a directory structure such as

/top/consurf_scripts
    /pdbidchain    e.g. 1fnfa
    /pdbidchain/file.pdb    - a single chain pdb file
    or
    /pdbidchain/file.fasta  - a single chain fasta file

Usage: consurf_home.sh file
       consurf_seqres file.pdb

To run it if you have a pdb file as:
   $consurf_scripts/consurf_home.sh file.pdb

To run it if you have a fasta file as:
   $consurf_scripts/consurf_home.sh file.fasta

If you want to use the SEQRES records from the pdb file instead of the ATOM records then:
   $consurf_scripts/consurf_seqres.sh file.fasta
Note the numbering will run from 1->X however

