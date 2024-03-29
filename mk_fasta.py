#!/usr/bin/python3
###################################################################################################
## Jon Wright, IBMS, Academia Sinica, Taipei, 11529, Taiwan
## These files are licensed under the GLP ver 3, essentially you have the right
## to copy, modify and distribute this script but all modifications must be offered
## back to the original authors
###################################################################################################
#
# This take a PDB file and then creates a file with a fasta format for a single chain
# It will ouput ACE/NME as X
# and write another file that allows for the final consurf numbering
# to be related back to the initial PDB file

import sys
#import os

# Set initial values for a few variables
PREVIOUS_RESNO = "-999"
PREVIOUS_RESNO_SHORT = "-999"

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
    'CYM' : 'C',
    'HSP' : 'H',
    'HSD' : 'H',
    'HSE' : 'H',
    'HIP' : 'H',
    'HID' : 'H',
    'HIE' : 'H',
    'ACE' : '',
    'NME' : '',
    ' ZN' : '',
}

# Get the output name from the command line
if len(sys.argv) == 1:
    print("Needs a argument with the list of PDB files to extract the chain from\n")
    sys.exit(0)

# Open the file,
INFILE = open(sys.argv[1], "r")
OUTFILE = open("r4s_pdb.py", "w")
print(">PDB_ATOM")
INDEX = 0
R4SINDEX = 0
OUTFILE.write("R4S_2_PDB={\n")
# Loop over each line of the input file splitting it into different columns
for TMLINE in INFILE:
    if TMLINE[0:4] == "ATOM" or TMLINE[0:6] == "HETATM":
    # Save the info we need for later
        resname = TMLINE[17:20]
        atom_name = TMLINE[12:16]
        chain = TMLINE[21:22]
        RESID_LONG = str(TMLINE[22:27])
        # This block looks for a change in the residue number  , and prints out the 1 letter codes
        if PREVIOUS_RESNO != RESID_LONG:
            INDEX += 1
            # Make a note of the residue number for the next line we look at
            PREVIOUS_RESNO = RESID_LONG
            # Convert from 3 letter code to single letter code and print it out
            if resname in RESIDUETAB:
                SINGLE = str(RESIDUETAB.get(resname))
                print(SINGLE, end="")
            else:
                print(".", end="")
            if SINGLE != "":
                R4SINDEX += 1
                outline = "    '" + "{:>d}".format(R4SINDEX) + "' : '" + "{:>d}".format(INDEX) \
                       + "' ,\n"
                OUTFILE.write(outline)
print("") # print a newline to end the sequence output
OUTFILE.write("}")
