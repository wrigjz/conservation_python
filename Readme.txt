#########################################################################################################################
## Jon Wright, IBMS, Academia Sinica, Taipei, 11529, Taiwan
## These files are licensed under the GLP ver 3, essentially you have the right
## to copy, modify and distribute this script but all modifications must be offered
## back to the original authors
#########################################################################################################################
#
This is an implimentation of the consurf system
It takes a single chain pdb file (often provided from critires) or a single chain fasta file
It uses 150 evenly sampled acceptable homologues 

To run you need a directory structure such as

/top/consurf_scripts
    /pdbidchain    e.g. 1fnfa
    /pdbidchain/file.pdb    - a single chain pdb file
    or
    /pdbidchain/file.fasta  - a single chain fasta file

Usage: consurf_home.sh file P|F
    where P is for pdb files and F is for fasta files

To run it if you have a pdb file as:
   $consurf_scripts/consurf_home.sh file.pdb P

To run it if you have a fasta file as:
   $consurf_scripts/consurf_home.sh file.fasta F


