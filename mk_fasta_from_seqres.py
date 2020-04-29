#!/usr/bin/python3
###################################################################################################
## Jon Wright, IBMS, Academia Sinica, Taipei, 11529, Taiwan
## These files are licensed under the GLP ver 3, essentially you have the right
## to copy, modify and distribute this script but all modifications must be offered
## back to the original authors
###################################################################################################
#
# This take a PDB file and then creates a file with a fasta format for a single chain
# It takes the sequence from the SEQRES recorss not from the ATOM ones
#
# Usage:
# python3 get_fasta_seqres.py file.pdb CHAINID
# e.g. python3 get_fasta_seqres.py original.pdb A

import sys
#import os

# Look up table for converting 3 letter AA codes to single letter codes
RESIDUETAB = {
    'ALA' : 'A',
    'CYS' : 'C',
    'ASP' : 'D',
    'GLU' : 'E',
    'PHE' : 'F',
    'GLY' : 'G',
    'HIS' : 'H',
    'ILE' : 'I',
    'LYS' : 'K',
    'LEU' : 'L',
    'MET' : 'M',
    'ASN' : 'N',
    'PRO' : 'P',
    'GLN' : 'Q',
    'ARG' : 'R',
    'SER' : 'S',
    'THR' : 'T',
    'VAL' : 'V',
    'TYR' : 'Y',
    'TRP' : 'W',
    'CYX' : 'C',
    'HSP' : 'H',
    'HSD' : 'H',
    'HSE' : 'H',
    'HIP' : 'H',
    'HID' : 'H',
    'HIE' : 'H',
    'ACE' : '',
    'NME' : '',
}

# Get the output name from the command line
if len(sys.argv) != 3:
    print("Needs a argument with the PDB file and chain ID to extract from\n")
    sys.exit(0)

TEMP = sys.argv[2]
CHAIN = TEMP.upper()

# Open the file,
INFILE = open(sys.argv[1], "r")
print(">PDB_ATOM")
# Loop over each line of the input file splitting it into different words
for TMLINE in INFILE:
    if TMLINE[0:6] == "SEQRES":
        SEQRES, LINE, CHAININ, COUNTER, *RESIDUE = [x.strip() for x in TMLINE.split()]
        if CHAININ == CHAIN:  # Only fo this for our chain
            residue_number = len(RESIDUE)
            for i in range(0, residue_number):
                #print(RESIDUE[i])
                if RESIDUE[i] in RESIDUETAB:
                    single = str(RESIDUETAB.get(RESIDUE[i]))
                    print(single, end="")
                else:
                    print(".", end="")
