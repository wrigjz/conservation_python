#!/usr/bin/python3
###################################################################################################
## Jon Wright, IBMS, Academia Sinica, Taipei, 11529, Taiwan
## These files are licensed under the GLP ver 3, essentially you have the right
## to copy, modify and distribute this script but all modifications must be offered
## back to the original authors
##################################################################################################
#
# This script reads a target input (arg1) and a cdhit fasta output file (arg2)
# it then does a pairwise alignment, it accepts alignmensts that are between
# 35-95% seqid with the target  and have a length 60% or over compared with
# with the target
# It also checks if the same homologue is used twice and if there is an overlap
# of > 10% the shorter sequence is rejected
# It then prints out 300 evenly sampled hits, if there are less than 300 it prints all
# of thosee too to accepted.fasta

# usage:
# python3 select_seqs.py reference.fasta cdhit.out

import sys
#import shutil
from Bio import pairwise2

if len(sys.argv) != 3:
    print("Please give reference fasta file and the homologs fasta file.")
    exit()

TARGETFILE = open(sys.argv[1], "r")
INFILE = open(sys.argv[2], "r")
PREALIGN = open("accepted.fasta", "w")
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
FAST_OPTION = 0
for i in range(1, INDEX):
    if REJECT[i] == "Acceptable": # if it's not an acceptbale one then don't check
        for j in range(i+1, INDEX+1):
            if REJECT[j] == "Acceptable": # If it's not acceptable one then don't check
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
                        REJECT[shorter_seq] = "Overlaps with itself elsewhere: " \
                            + str(percentage_digits)
        if REJECT[i] == "Acceptable": # This one is still acceptable so increment the counter
            FAST_OPTION += 1
        if FAST_OPTION == 300: # We have found enough, lets stop checking and go to writing out
            break

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

# If there are not enough homologs, min of 5, we stop at this point
if REMAINING < 6:
    print("Not enough homologs for conservation analysis: ", INDEX, REMAINING)
    exit()

# Write out the target SEQUENCE first to the accepted file
PREALIGN.write(MERGED_LIST_ACC[0][0])
PREALIGN.write("\n")
PREALIGN.write(MERGED_LIST_ACC[0][1])
PREALIGN.write("\n")

# If there are < 300 remaining, write them all out, start at 1 because the target is  0
if REMAINING <= 300:
    for i in range(1, REMAINING):
        PREALIGN.write(MERGED_LIST_ACC[i][0])
        PREALIGN.write("\n")
        PREALIGN.write(MERGED_LIST_ACC[i][1])
        PREALIGN.write("\n")

# Write out the top 300 file if there are more than 300 accepted homologue
if REMAINING > 300:
    for i in range(1, 301):
        PREALIGN.write(MERGED_LIST_ACC[i][0])
        PREALIGN.write("\n")
        PREALIGN.write(MERGED_LIST_ACC[i][1])
        PREALIGN.write("\n")
PREALIGN.close() # Close the prealign file

# Loop back over the rejected sequences and print them out with a reason
for i in range(0, REJECTED):
    REJECTFILE.write(MERGED_LIST_REJ[i][0])
    REJECTFILE.write("\n")
    REJECTFILE.write(MERGED_LIST_REJ[i][2])
    REJECTFILE.write("\n")
    REJECTFILE.write(MERGED_LIST_REJ[i][1])
    REJECTFILE.write("\n")

# Now write to the reject file any that were accetable but we already had enough
for i in range(301, REMAINING):
    REJECTFILE.write(MERGED_LIST_ACC[i][0])
    REJECTFILE.write("\n")
    REJECTFILE.write("Was acceptable but we already had 300\n")
    REJECTFILE.write(MERGED_LIST_ACC[i][1])
    REJECTFILE.write("\n")
REJECTFILE.close() # Close the reject file
