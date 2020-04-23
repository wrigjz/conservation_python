#!/usr/bin/python3
###################################################################################################
## and Jon Wright, IBMS, Academia Sinica, Taipei, 11529, Taiwan
## These files are licensed under the GLP ver 3, essentially you have the right
## to copy, modify and distribute this script but all modifications must be offered
## back to the original authors
###################################################################################################


# This script reads a target input (arg1) and the aligned sequences
# It then reports how many residues are aligned with each residue in the PDB_ATOM line
# It also keeps a count of how common each residue type is at each point and prints them out
# in order of most common first

# usage:
# $mafftdir/bin/mafft-linsi --quiet --localpair --maxiterate 1000 --reorder --namelength 30  \
#        unaligned.fast.fasta >| frequency.fasta
# python3 get_frequency.py frequency.fasta >| frequency.txt

import sys
#import os

if len(sys.argv) != 2:
    print("Please give the aligned fasta file.")
    exit()

INFILE = open(sys.argv[1], "r")

# Intial array setup
TITLE = []
SEQUENCE = []
POSITION = []
ALA = []
CYS = []
ASP = []
GLU = []
PHE = []
GLY = []
HIS = []
ILE = []
LYS = []
LEU = []
MET = []
ASN = []
PRO = []
GLN = []
ARG = []
SER = []
THR = []
VAL = []
TRP = []
TYR = []
NAME = ["" for i in range(0, 20)]
TEMP = ["" for i in range(0, 20)]
INDEX = -1

# Setup single latter codes for the tuple later
NAME[0] = "A"
NAME[1] = "C"
NAME[2] = "D"
NAME[3] = "E"
NAME[4] = "F"
NAME[5] = "G"
NAME[6] = "H"
NAME[7] = "I"
NAME[8] = "K"
NAME[9] = "L"
NAME[10] = "M"
NAME[11] = "N"
NAME[12] = "P"
NAME[13] = "Q"
NAME[14] = "R"
NAME[15] = "S"
NAME[16] = "T"
NAME[17] = "V"
NAME[18] = "W"
NAME[19] = "Y"

# Read in the fasta file
for LINE in INFILE:
    LINE = LINE.rstrip("\n") # remove the newline
    if LINE[0:1] == ">":
        TITLE.append(0)
        SEQUENCE.append("")
        POSITION.append("N")  # N for None - if there are amino acids it will become an integer
        INDEX += 1    # increment the index
        TITLE[INDEX], *junk = [x.strip() for x in LINE.split("/")]
    else:
        SEQUENCE[INDEX] = SEQUENCE[INDEX] + LINE
INFILE.close()

# Set up the arrays for holding the residue counts
# Each position in the target sequence needs to be able to have a variable for each of the
# 20 amino acids
ALIGNED_LENGTH = len(SEQUENCE[0])
for i in range(0, ALIGNED_LENGTH):
    ALA.append(0)
    CYS.append(0)
    ASP.append(0)
    GLU.append(0)
    PHE.append(0)
    GLY.append(0)
    HIS.append(0)
    ILE.append(0)
    LYS.append(0)
    LEU.append(0)
    MET.append(0)
    ASN.append(0)
    PRO.append(0)
    GLN.append(0)
    ARG.append(0)
    SER.append(0)
    THR.append(0)
    VAL.append(0)
    TRP.append(0)
    TYR.append(0)

# At this point look for the non "-" positions in the target sequence and count how many other
# sequences are not "-" too
#print(SEQUENCE[0])
for i in range(0, ALIGNED_LENGTH):     # Loop over each residues in the target sequence
    if SEQUENCE[0][i] != "-" and SEQUENCE[0][i] != "X":
        #print(SEQUENCE[0][i])
        POSITION[i] = 0
        for j in range(1, INDEX+1):    # Loop over all the other aligned homolog sequences
            if SEQUENCE[j][i] != "-":
                POSITION[i] += 1      # Keeps track of the number of other seqs with an AA
                                      # at this place
            if SEQUENCE[j][i] == "A":
                ALA[i] += 1
            if SEQUENCE[j][i] == "C":
                CYS[i] += 1
            if SEQUENCE[j][i] == "D":
                ASP[i] += 1
            if SEQUENCE[j][i] == "E":
                GLU[i] += 1
            if SEQUENCE[j][i] == "F":
                PHE[i] += 1
            if SEQUENCE[j][i] == "G":
                GLY[i] += 1
            if SEQUENCE[j][i] == "H":
                HIS[i] += 1
            if SEQUENCE[j][i] == "I":
                ILE[i] += 1
            if SEQUENCE[j][i] == "K":
                LYS[i] += 1
            if SEQUENCE[j][i] == "L":
                LEU[i] += 1
            if SEQUENCE[j][i] == "M":
                MET[i] += 1
            if SEQUENCE[j][i] == "N":
                ASN[i] += 1
            if SEQUENCE[j][i] == "P":
                PRO[i] += 1
            if SEQUENCE[j][i] == "Q":
                GLN[i] += 1
            if SEQUENCE[j][i] == "R":
                ARG[i] += 1
            if SEQUENCE[j][i] == "S":
                SER[i] += 1
            if SEQUENCE[j][i] == "T":
                THR[i] += 1
            if SEQUENCE[j][i] == "V":
                VAL[i] += 1
            if SEQUENCE[j][i] == "W":
                TRP[i] += 1
            if SEQUENCE[j][i] == "Y":
                TYR[i] += 1

# Print out the resuls for the non "-" positions
print("Seq: Res Fnd / Aln  Residues in order of being observed")
for i in range(0, ALIGNED_LENGTH):     # Loop over each residues in the target sequence
    if POSITION[i] != "N":
        TEMP[0] = ALA[i]
        TEMP[1] = CYS[i]
        TEMP[2] = ASP[i]
        TEMP[3] = GLU[i]
        TEMP[4] = PHE[i]
        TEMP[5] = GLY[i]
        TEMP[6] = HIS[i]
        TEMP[7] = ILE[i]
        TEMP[8] = LYS[i]
        TEMP[9] = LEU[i]
        TEMP[10] = MET[i]
        TEMP[11] = ASN[i]
        TEMP[12] = PRO[i]
        TEMP[13] = GLN[i]
        TEMP[14] = ARG[i]
        TEMP[15] = SER[i]
        TEMP[16] = THR[i]
        TEMP[17] = VAL[i]
        TEMP[18] = TRP[i]
        TEMP[19] = TYR[i]
        # Form the single letter codes and their associated values into a tuple so we can
        # remove 0's and sort it from most frequent to least
        MERGED_LIST_ORIG = tuple(zip(NAME, TEMP))
        MERGED_LIST = list(filter(lambda a: a[1] != 0, MERGED_LIST_ORIG)) # Remove 0's  from list
        MAX_ARRAY = sorted(MERGED_LIST, key=lambda a: a[1], reverse=True) # Do the sorting

        # Finally print the results
        print("Seq:", "{:>2}".format(SEQUENCE[0][i]), "{:>4}".format(POSITION[i]), "/",\
              "{:>4}".format(INDEX), end=" ")
        res_list = ""
        array_len = len(MAX_ARRAY)
        for j in range(0, array_len):
            res_list = res_list + MAX_ARRAY[j][0] + ", "
        print(res_list)
