#!/usr/bin/python3
###################################################################################################
## Jon Wright, IBMS, Academia Sinica, Taipei, 11529, Taiwan
## These files are licensed under the GLP ver 3, essentially you have the right
## to copy, modify and distribute this script but all modifications must be offered
## back to the original authors
###################################################################################################
#
# Read the rete4site result files, then
# Split the scores into grades, this is done by takeing the (already normalised) scores
# and swetting the most negative score as equal to-4.5 bin WIDTHs, we then assign
# all the numbers into equal bin WIDTHs -4.5->-3.5 ... -0.5->+0.5
# with the last (consurf_grade 9) being open ended so +3.5->infinity

import sys
import os

# Read in the original PDB numbering scheme if the file exists
if os.path.exists('r4s_pdb.py'):
    from r4s_pdb import R4S_2_PDB
    WAS_A_PDB = 1
else:
    WAS_A_PDB = 0

INFILE = open(sys.argv[1], "r")
OUTFILE = open(sys.argv[2], "w")

if sys.argv[3] == "N":
    WAS_A_PDB = 0

RESIDUETAB = {
    'A' : 'ALA',
    'C' : 'CYS',
    'D' : 'ASP',
    'E' : 'GLU',
    'F' : 'PHE',
    'G' : 'GLY',
    'H' : 'HIS',
    'I' : 'ILE',
    'K' : 'LYS',
    'L' : 'LEU',
    'M' : 'MET',
    'N' : 'ASN',
    'P' : 'PRO',
    'Q' : 'GLN',
    'R' : 'ARG',
    'S' : 'SER',
    'T' : 'THR',
    'V' : 'VAL',
    'Y' : 'TYR',
    'W' : 'TRP',
}

LOWEST = 0
for LINE in INFILE:
    if LINE[0:1] == "#" or LINE == "\n": # look for lines without leading # and are not blank
        continue
    resnum, resname, in3, *junk = [x.strip() for x in LINE.split()]
    score = float(in3)
    if score < LOWEST:
        LOWEST = score
INFILE.seek(0)

# Bin WIDTH is set so that the LOWEST number is 4.5 WIDTHs from 0
WIDTH = LOWEST / -4.5
#print(LOWEST, WIDTH)

# Write out the output file headers
OUTLINE = "  # SEQ 3LETT PDB COLOUR  SCORE \n"
OUTFILE.write(OUTLINE)

#Set up the 9 bins ranges - not bin 1 doesn't have an upper range
for LINE in INFILE:
    if LINE[0:1] == "#" or LINE == "\n": # look for lines without leading # and are not blank
        continue
    resnum, resname, in3, *junk = [x.strip() for x in LINE.split()] # Read in Res Name/Num/Score
    score = float(in3)
    if WAS_A_PDB == 1:
        original = (R4S_2_PDB.get(resnum)) # get the 'original' resnumber from the initial pdb file
    else:
        original = resnum # We had a fasta file so no original number
    threeletter = (RESIDUETAB.get(resname)) # Get the three letter code

    # now take each score and test it against the coloured bins for consurf grades
    # if it matches then print it out with the grade
    for i in range(0, 9):
        lower = LOWEST + (i * WIDTH)
        upper = lower + WIDTH
        binnumber = 9 - i
        if lower <= score < upper: # find bins for residues
            #print(resnum,resname,threeletter,original,score,binnumber)
            OUTLINE = "{:>4s}  ".format(resnum) + "{:1s}  ".format(resname) + \
                      "{:3s} ".format(threeletter) + "{:>4s}    ".format(original) + \
                      "{:d}   ".format(binnumber) + "{:>6.3f} ".format(score) + "\n"
            OUTFILE.write(OUTLINE)
        if binnumber == 1 and score >= upper: # Catch the ones above 1
            #print(resnum,resname,threeletter,original,score,binnumber)
            OUTLINE = "{:>4s}  ".format(resnum) + "{:1s}  ".format(resname) + \
                      "{:3s} ".format(threeletter) + "{:>4s}    ".format(original) + \
                      "{:d}   ".format(binnumber) + "{:>6.3f} ".format(score) + "\n"
            OUTFILE.write(OUTLINE)
