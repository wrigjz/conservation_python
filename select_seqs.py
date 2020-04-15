#!/usr/bin/python3

# This script reads a target input (arg1) and a cdhit fasta output file (arg2)
# it then does a pairwise alignment, it accepts alignmensts that are between
# 35-95% seqid with the target  and have a length 60% or over compared with
# with the target
# It also checks if the same homologue is used twice and if there is an overlap
# of > 10% the shorter sequence is rejected
# It then prints out the top 150 hits, if there are less than 150 it prints all
# of thosee too to accepted150top.fasta
# It then prints out the top 300 hits, if there are less than 300 it prints all
# of thosee too to accepted300top.fasta
# It then prints out an even spread of 150 hits to the accetped150sp.fasta file
# It then prints out an even spread of 300 hits to the accetped300sp.fasta file

# usage:
# python3 select_seqs.py reference.fasta cdhit.out

import sys
import shutil
from Bio import pairwise2

if len(sys.argv) != 3:
    print("Please give reference fasta file and the homologues fasta file.")
    exit()

TARGETFILE = open(sys.argv[1], "r")
INFILE = open(sys.argv[2], "r")
PREALIGN150TOP = open("accepted150top.fasta", "w")
PREALIGN300TOP = open("accepted300top.fasta", "w")
PREALIGN150SP = open("accepted150sp.fasta", "w")
PREALIGN300SP = open("accepted300sp.fasta", "w")
REJECTFILE = open("rejected.fasta", "w")

# Intial setup is for 100000 sequences
TITLE = []
SHORT_TITLE = []
SEQUENCE = []
REJECT = []  # A number other than 0 is the percentage for rejection
             # Acceptable means it's not a reject
INDEX = -1

# Read in the target seqeunce to TITLE[0] and SEQUENCE[0]
for LINE in TARGETFILE:
    LINE = LINE.rstrip("\n") # remove the newline
    if LINE[0:1] == ">":
        TITLE.append(0)
        SHORT_TITLE.append(0)
        SEQUENCE.append("")
        REJECT.append("Acceptable")
        INDEX += 1    # increment the index
        TITLE[INDEX] = LINE
        SHORT_TITLE[INDEX] = LINE
    else:
        SEQUENCE[INDEX] = SEQUENCE[INDEX] + LINE
TARGETFILE.close()
TARGET_LEN = len(SEQUENCE[0])

# Now read in the cdhit output and save to TITLE[1] ++ and SEQUENCE[1] ++ and SHORT_TITLE[1] ++
for LINE in INFILE:
    LINE = LINE.rstrip("\n") # remove the newline
    if LINE[0:1] == ">":
        TITLE.append(0)
        SHORT_TITLE.append(0)
        SEQUENCE.append("")
        REJECT.append("Acceptable")
        INDEX += 1    # increment the index
        TITLE[INDEX] = LINE
        SHORT_TITLE[INDEX], *junk = [x.strip() for x in LINE.split("/")]
    else:
        SEQUENCE[INDEX] = SEQUENCE[INDEX] + LINE
INFILE.close()

# Now start to sort out the sequences we want to keep
# Loop over the lot and remove those that don't meet the criteria against the target sequence
for i in range(1, INDEX+1):
    alignments = pairwise2.align.globalxx(SEQUENCE[0], SEQUENCE[i])
    score = alignments[0][2]
    percentage = (score/TARGET_LEN) * 100
    overlap = (len(SEQUENCE[i]) / TARGET_LEN) *100
    percentage_digits = round(percentage, 2)
    overlap_digits = round(overlap, 2)
    if percentage <= float(35): # Reject not enough seqid
        REJECT[i] = "Too Low: " + str(percentage_digits)
    elif percentage >= float(95): # Reject too much seqid
        REJECT[i] = "Too High: " + str(percentage_digits)
    elif overlap <= float(60): # Reject for lack of overlap
        REJECT[i] = "Too Short: " + str(overlap_digits)

# At this point we need to check the short Uniref names if there is a match with another
# Get the two lenghts, align shorter with longer, get percentage
# if > 10% we reject the shorter
#>UniRef90_P13726/37-243 [subseq from] Tissue factor n=43 Tax=Catarrhini TaxID=9526 RepID=TF_HUMAN
#>UniRef90_A0A2I3REM1/21-227 [subseq from] Coagulation factor III, tissue factor n=1
# We use local algnment becasue we are looking for overgaps not a full-algnment
# We set the mismatch and gap penalties to -10 because they would not exist in an overlap
# percentage >= 10 means too much overkap so we reject the shorted chain
for i in range(1, INDEX+1):
    for j in range(i+1, INDEX+1):
        if SHORT_TITLE[i] == SHORT_TITLE[j]:
            length1 = len(SEQUENCE[i])
            length2 = len(SEQUENCE[j])
            # Now do the choosing only pick #2 is if is longer than #1
            if length2 > length1:
                longer = length2
                shorter_seq = i
                alignments = pairwise2.align.localms(SEQUENCE[i], SEQUENCE[j], 1, -10, -10, -1)
            else:
                longer = length1
                shorter_seq = j
                alignments = pairwise2.align.localms(SEQUENCE[j], SEQUENCE[i], 1, -10, -10, -1)
            score = alignments[0][2]
            percentage = (score/longer) * 100
            percentage_digits = round(percentage, 2)
            if percentage >= 10:
                REJECT[shorter_seq] = "Overlaps with itself elsewhere: " + str(percentage_digits)

# At this point we need to create a new list of those that are acceptable
# so we can loop over these later looking for whatever sets we want to make
# we make this as a tuple
MERGED_LIST_ORIG = tuple(zip(TITLE, SEQUENCE, REJECT))
# List only Acceptable values
MERGED_LIST_ACC = list(filter(lambda a: a[2] == "Acceptable", MERGED_LIST_ORIG))
# List only non Acceptable values
MERGED_LIST_REJ = list(filter(lambda a: a[2] != "Acceptable", MERGED_LIST_ORIG))
REMAINING = len(MERGED_LIST_ACC)
REJECTED = len(MERGED_LIST_REJ)

# If there are not enough homologues, min of 5, we stop at this point
if REMAINING < 4:
    print("Not enough homologues for conservation analysis: ", INDEX, REMAINING)
    exit()

# Write out the target SEQUENCE first to each of the accepted files
PREALIGN150TOP.write(MERGED_LIST_ACC[0][0])
PREALIGN150TOP.write("\n")
PREALIGN150TOP.write(MERGED_LIST_ACC[0][1])
PREALIGN150TOP.write("\n")
PREALIGN300TOP.write(MERGED_LIST_ACC[0][0])
PREALIGN300TOP.write("\n")
PREALIGN300TOP.write(MERGED_LIST_ACC[0][1])
PREALIGN300TOP.write("\n")
PREALIGN150SP.write(MERGED_LIST_ACC[0][0])
PREALIGN150SP.write("\n")
PREALIGN150SP.write(MERGED_LIST_ACC[0][1])
PREALIGN150SP.write("\n")
PREALIGN300SP.write(MERGED_LIST_ACC[0][0])
PREALIGN300SP.write("\n")
PREALIGN300SP.write(MERGED_LIST_ACC[0][1])
PREALIGN300SP.write("\n")

# If there are < 150 remaining, write them all out
if REMAINING <= 150:
    for i in range(1, REMAINING):
        PREALIGN150TOP.write(MERGED_LIST_ACC[i][0])
        PREALIGN150TOP.write("\n")
        PREALIGN150TOP.write(MERGED_LIST_ACC[i][1])
        PREALIGN150TOP.write("\n")
    # Close the other files and copy this file to all the other accepted files
    PREALIGN150TOP.close()
    PREALIGN150SP.close()
    PREALIGN300TOP.close()
    PREALIGN300SP.close()
    shutil.copyfile("accepted150top.fasta", "accepted150sp.fasta")
    shutil.copyfile("accepted150top.fasta", "accepted300top.fasta")
    shutil.copyfile("accepted150top.fasta", "accepted300sp.fasta")

# Write out the top 150 file if there are more than 150 accepted homologue
if REMAINING > 150:
    for i in range(1, 151): # Write out the top 150
        PREALIGN150TOP.write(MERGED_LIST_ACC[i][0])
        PREALIGN150TOP.write("\n")
        PREALIGN150TOP.write(MERGED_LIST_ACC[i][1])
        PREALIGN150TOP.write("\n")
    # Here we try to write out 150 that evenly sample the remaining homologues
    INTERVAL = int(REMAINING /150) # Calculate the interval for the 150 samples
    for i in range(1, REMAINING, INTERVAL):
        PREALIGN150SP.write(MERGED_LIST_ACC[i][0])
        PREALIGN150SP.write("\n")
        PREALIGN150SP.write(MERGED_LIST_ACC[i][1])
        PREALIGN150SP.write("\n")
    PREALIGN150TOP.close()
    PREALIGN150SP.close()

# If the number of accepted homologues is between 150 and 300 we just write out all of them
# and copy that to the 300 even spread file
if 150 < REMAINING <= 300:
    for i in range(1, REMAINING):
        PREALIGN300TOP.write(MERGED_LIST_ACC[i][0])
        PREALIGN300TOP.write("\n")
        PREALIGN300TOP.write(MERGED_LIST_ACC[i][1])
        PREALIGN300TOP.write("\n")
    # Close the 300 even spread file and copy this file to it
    PREALIGN300TOP.close()
    PREALIGN300SP.close()
    shutil.copyfile("accepted300top.fasta", "accepted300sp.fasta")

# Write out the top 300 file if there are more than 300 accetped homologues
if REMAINING > 300:
    for i in range(1, 301):  # Write out the top 300
        PREALIGN300TOP.write(MERGED_LIST_ACC[i][0])
        PREALIGN300TOP.write("\n")
        PREALIGN300TOP.write(MERGED_LIST_ACC[i][1])
        PREALIGN300TOP.write("\n")
    # Here we try to write out 300 that evenly sample the REMAINING homologues
    INTERVAL = int(REMAINING / 300) # Calculate the interval for the 300 samples
    for i in range(1, REMAINING, INTERVAL):
        PREALIGN300SP.write(MERGED_LIST_ACC[i][0])
        PREALIGN300SP.write("\n")
        PREALIGN300SP.write(MERGED_LIST_ACC[i][1])
        PREALIGN300SP.write("\n")
    PREALIGN300TOP.close()
    PREALIGN300SP.close()

# Loop back over the rejected sequences and print them out with a reason
for i in range(0, REJECTED):
    REJECTFILE.write(MERGED_LIST_REJ[i][0])
    REJECTFILE.write("\n")
    REJECTFILE.write(MERGED_LIST_REJ[i][2])
    REJECTFILE.write("\n")
    REJECTFILE.write(MERGED_LIST_REJ[i][1])
    REJECTFILE.write("\n")
REJECTFILE.close()
