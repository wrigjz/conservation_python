#!/usr/bin/python
#########################################################################################################################
## and Jon Wright, IBMS, Academia Sinica, Taipei, 11529, Taiwan
## These files are licensed under the GLP ver 3, essentially you have the right
## to copy, modify and distribute this script but all modifications must be offered
## back to the original authors
#########################################################################################################################
# Read the rete4site result files, then
# Split the scores into grades, this is done by takeing the (already normalised) scores
# and swetting the most negative score as equal to-4.5 bin widths, we then assign 
# all the numbers into equal bin widths -4.5->-3.5 ... -0.5->+0.5
# with the last (consurf_grade 9) being open ended so +3.5->infinity

import sys
import os

INFILE = open(sys.argv[1],"r")
OUTFILE = open(sys.argv[2],"w")

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

lowest = 0
for LINE in INFILE:
    if LINE[0:1] == "#" or LINE == "\n": # look for lines without leading # and are not blank
        continue
    resnum, resname, in3, *junk = [x.strip() for x in LINE.split()]
    score = float(in3)
    if score < lowest: lowest = score
INFILE.seek(0)

# Bin width is set so that the lowest number is 4.5 widths from 0
width = lowest / -4.5
#print(lowest, width)

# Read in the original PDB numbering scheme
from r4s_pdb import R4S_2_PDB 

# Write out the output file headers
outline="  # SEQ 3LETT PDB COLOUR  SCORE\n"
OUTFILE.write(outline)

#Set up the 9 bins ranges - not bin 1 doesn't have an upper range
for LINE in INFILE:
    if LINE[0:1] == "#" or LINE == "\n": # look for lines without leading # and are not blank
        continue
    resnum, resname, in3, *junk = [x.strip() for x in LINE.split()] # Read in Res Name/Num/Score
    score = float(in3)
    original = (R4S_2_PDB.get(resnum)) # get the 'original' resnumber from the initial pdb file for CS
    threeletter = (RESIDUETAB.get(resname)) # Get the three letter code

    # now take each score and test it against the coloured bins for consurf grades
    # if it matches then print it out with the grade
    for i in range(0,9):
        lower = lowest + (i * width)
        upper = lower + width
        binnumber = 9 - i
        if score >= lower and score < upper: # find bins for residues
            #print(resnum,resname,threeletter,original,score,binnumber)
            outline="{:>4s}  ".format(resnum) + "{:1s}  ".format(resname) + "{:3s} ".format(threeletter) + \
                    "{:>4s}    ".format(original) + "{:d}   ".format(binnumber) + "{:>6.3f} ".format(score) + "\n"
            OUTFILE.write(outline)
        if binnumber == 1 and score >= upper: # Catch the ones above 1
            #print(resnum,resname,threeletter,original,score,binnumber)
            outline="{:>4s}  ".format(resnum) + "{:1s}  ".format(resname) + "{:3s} ".format(threeletter) + \
                    "{:>4s}    ".format(original) + "{:d}   ".format(binnumber) + "{:>6.3f} ".format(score) + "\n"
                    

            OUTFILE.write(outline)

