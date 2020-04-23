#!/usr/bin/python3
#########################################################################################################################
## and Jon Wright, IBMS, Academia Sinica, Taipei, 11529, Taiwan
## These files are licensed under the GLP ver 3, essentially you have the right
## to copy, modify and distribute this script but all modifications must be offered
## back to the original authors
#########################################################################################################################

#$mafftdir/bin/mafft-linsi --quiet --localpair --maxiterate 1000 --reorder --namelength 30  ./accepted.fasta >| prop.fasta

# This script reads a target input (arg1) and the aligned sequences
# It then reports how many residues are aligned with each residue in the PDB_ATOM line
# It also keeps a count of how common each residue type is at each point

# usage:
# python3 get_propensities.py prop.fasta

import sys
import os

if len(sys.argv) != 2:
    print("Please give the aligned fasta file.")
    exit()

INFILE = open(sys.argv[1], "r")

# Intial setup
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
name = ["" for i in range(0,20)]
temp = ["" for i in range(0,20)]
INDEX = -1

# Setup single latter codes for the tuple later
name[0] = "A"
name[1] = "C"
name[2] = "D"
name[3] = "E"
name[4] = "F"
name[5] = "G"
name[6] = "H"
name[7] = "I"
name[8] = "K"
name[9] = "L"
name[10] = "M"
name[11] = "N"
name[12] = "P"
name[13] = "Q"
name[14] = "R"
name[15] = "S"
name[16] = "T"
name[17] = "V"
name[18] = "W"
name[19] = "Y"

# Read in the fasta file
for LINE in INFILE:
    LINE = LINE.rstrip("\n") # remove the newline
    if LINE[0:1] == ">":
        TITLE.append(0)
        SEQUENCE.append("")
        POSITION.append("N")
        INDEX += 1    # increment the index
        TITLE[INDEX], *junk = [x.strip() for x in LINE.split("/")]
    else:
        SEQUENCE[INDEX] = SEQUENCE[INDEX] + LINE
INFILE.close()

# Set upi the arrayes for holding the residue counts
aligned_length=len(SEQUENCE[0])
for i in range(0,aligned_length):
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

#Here we need to set up the indexes for each position

# At this point look for the non "-" positions and count how many other sequences are not "-" too
#print(SEQUENCE[0])
for i in range(0,aligned_length):     # Loop over each residues in the target sequence
    if SEQUENCE[0][i] != "-" and SEQUENCE[0][i] != "X":
        #print(SEQUENCE[0][i])
        POSITION[i] = 0
        for j in range(1,INDEX+1):    # Loop over all the other sequences
            if SEQUENCE[j][i] != "-":
                POSITION[i] += 1
            if SEQUENCE[j][i] == "A": ALA[i] += 1
            if SEQUENCE[j][i] == "C": CYS[i] += 1
            if SEQUENCE[j][i] == "D": ASP[i] += 1
            if SEQUENCE[j][i] == "E": GLU[i] += 1
            if SEQUENCE[j][i] == "F": PHE[i] += 1
            if SEQUENCE[j][i] == "G": GLY[i] += 1
            if SEQUENCE[j][i] == "H": HIS[i] += 1
            if SEQUENCE[j][i] == "I": ILE[i] += 1
            if SEQUENCE[j][i] == "K": LYS[i] += 1
            if SEQUENCE[j][i] == "L": LEU[i] += 1
            if SEQUENCE[j][i] == "M": MET[i] += 1
            if SEQUENCE[j][i] == "N": ASN[i] += 1
            if SEQUENCE[j][i] == "P": PRO[i] += 1
            if SEQUENCE[j][i] == "Q": GLN[i] += 1
            if SEQUENCE[j][i] == "R": ARG[i] += 1
            if SEQUENCE[j][i] == "S": SER[i] += 1
            if SEQUENCE[j][i] == "T": THR[i] += 1
            if SEQUENCE[j][i] == "V": VAL[i] += 1
            if SEQUENCE[j][i] == "W": TRP[i] += 1
            if SEQUENCE[j][i] == "Y": TYR[i] += 1

# Print out the resuls for the non "-" positions
print("Seq: Res Fnd / Aln  Residues in order of being observed")
for i in range(0,aligned_length):     # Loop over each residues in the target sequence
    if POSITION[i] != "N":
        temp[0] = ALA[i]
        temp[1] = CYS[i]
        temp[2] = ASP[i]
        temp[3] = GLU[i]
        temp[4] = PHE[i]
        temp[5] = GLY[i]
        temp[6] = HIS[i]
        temp[7] = ILE[i]
        temp[8] = LYS[i]
        temp[9] = LEU[i]
        temp[10] = MET[i]
        temp[11] = ASN[i]
        temp[12] = PRO[i]
        temp[13] = GLN[i]
        temp[14] = ARG[i]
        temp[15] = SER[i]
        temp[16] = THR[i]
        temp[17] = VAL[i]
        temp[18] = TRP[i]
        temp[19] = TYR[i]
        # Form the single letter codes and their associated values into a tuple so we can remove 0's and sort it
        MERGED_LIST_ORIG = tuple(zip(name, temp))
        MERGED_LIST = list(filter(lambda a: a[1] != 0, MERGED_LIST_ORIG)) # Remove 0 values from list
        MAX_ARRAY = sorted(MERGED_LIST, key=lambda a: a[1], reverse=True)

        print("Seq:", "{:>2}".format(SEQUENCE[0][i]), "{:>4}".format(POSITION[i]), "/", "{:>4}".format(INDEX), end=" ")
        res_list = ""
        array_len=len(MAX_ARRAY)
        for j in range(0,array_len):
            res_list = res_list + MAX_ARRAY[j][0] + ", "
        print(res_list)
