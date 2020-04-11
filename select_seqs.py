#!/usr/bin/python3

# This script read a target input (arg1) and a cdhit fasta output file (arg2)
# it then does a pairwise alignment, it accepts alignmensts that are between 35-95% seqid with the target
# and have a length 60% or over compared with with the target
# It also checks if the same homologue is used twice and if there is an overlap of > 10% the shorter
# sequence is rejected
# It then prints out the top 150 hits, if there are less than 150 left it prints all of thosee too
# if there are more than 150 left then it picks 150 of those left-overs randomly

# usage:
# python3 select_seqs.py reference.fasta homologues.fasta [150|300]

import os
import sys
import random
import decimal
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

if len(sys.argv) < 4:
    print ("Please give reference fasta file, the homologues fasta file and then either 150 or 300 depending on if you want just the top 150 matches or the top 150 + another 150 randomly choosen matches")
    exit()
if sys.argv[3] == "150":
    print("Taking the top 150 matches")
elif sys.argv[3] == "300":
    print("Taking the top 150 and then another 150 random ones from the other matches")
else:
    print ("Please give reference fasta file, the homologues fasta file and then either 150 or 300 depending on if you want just the top 150 matches or the top 150 + another 150 randomly choosen matches")
    exit()

TARGETFILE=open(sys.argv[1],"r")
INFILE=open(sys.argv[2],"r")
METHOD=sys.argv[3]
PREALIGN=open("accepted.fasta","w")
REJECT=open("rejected.fasta","w")

# Intial setup is for 100000 sequences
title       = []
short_title = []
sequence    = []
reject      = []  # Done = Accepted and written, Reconsider = available for reconsideration, 
                                              # A number other than 0 is the percentage for rejection
                                              # Acceptable means it's not a reject
index       = -1
total_aligns_written = 0
reconsider = 0

# Read in the target seqeunce to title[0] and sequence[0]
for LINE in TARGETFILE:
    LINE = LINE.rstrip("\n") # remove the newline
    if LINE[0:1] == ">":
        title.append(0)
        short_title.append(0)
        sequence.append("")
        reject.append("acceptable")
        index += 1    # increment the index
        title[index]       = LINE
        short_title[index] = LINE
    else:
        sequence[index] = sequence[index] + LINE
TARGETFILE.close()
target_len = len(sequence[0])

# Write out the target sequence first
PREALIGN.write(title[0])
PREALIGN.write("\n")
PREALIGN.write(sequence[0])
PREALIGN.write("\n")

# Now read in the cdhit output and save to title[1] ++ and sequence[1] ++ and short_title[1] ++
for LINE in INFILE:
    LINE = LINE.rstrip("\n") # remove the newline
    if LINE[0:1] == ">":
        title.append(0)
        short_title.append(0)
        sequence.append("")
        reject.append("Acceptable")
        index += 1    # increment the index
        title[index] = LINE
        short_title[index], *junk = [x.strip() for x in LINE.split("/")]
    else:
        sequence[index] = sequence[index] + LINE
INFILE.close()

# Now sort out the sequences we want to keep
# Loop over the lot and remove those that don't meet the criteria 
for i in range(1,index+1):
    #if total_aligns_written < int(150):  # We only write 150 out maximum
        alignments = pairwise2.align.globalxx(sequence[0], sequence[i])
        score = alignments[0][2]
        percentage = (score/target_len) * 100
        overlap = (len(sequence[i] ) / target_len) *100
        percentage_digits = round(percentage,2)
        overlap_digits = round(overlap,2)
        if percentage <= float(35):   # Reject not enough seqid
            reject[i] = "Too Low: " + str(percentage_digits)
        elif percentage >= float(95): # Reject too much seqid
            reject[i] = "Too High: " + str(percentage_digits)
        elif overlap <= float(60):    # Reject for lack of overlap
            reject[i] = "Too Short: " + str(overlap_digits)

# At this point we need to check the short Uniref anmesm if there is a match with another
# Get the two lenghts, align shorter with longer, get percentage
# if > 10% we reject the shorter
#>UniRef90_P13726/37-243 [subseq from] Tissue factor n=43 Tax=Catarrhini TaxID=9526 RepID=TF_HUMAN
#>UniRef90_A0A2I3REM1/21-227 [subseq from] Coagulation factor III, tissue factor n=1 Tax=Pan troglodytes TaxID=9598 RepID=A0A2I3REM1_
# We use local algnment becasue we are looking for overgaps not a full-algnment
# We set the mismatch and gap penalties to -10 because they would not exist in an overlap
# percentage >= 10 means too much overkap so we reject the shorted chain
for i in range(1,index+1):
    for j in range(i+1,index+1):
        if short_title[i] == short_title[j]:
            length1=len(sequence[i])
            length2=len(sequence[j])
            # Now do the choosing only pick #2 is if is longer than #1
            if length2 > length1: 
                longer = length2
                shorter_seq = i
                alignments = pairwise2.align.localms(sequence[i], sequence[j],1,-10,-10,-1)
            else:
                longer = length1
                shorter_seq = j
                alignments = pairwise2.align.localms(sequence[j], sequence[i],1,-10,-10,-1)
            score = alignments[0][2]
            percentage = (score/longer) * 100
            percentage_digits = round(percentage,2)
            if percentage >= 10:
                reject[shorter_seq] = "Overlaps with itself elsewhere: " + str(percentage_digits)

# Loop back over the sequences and start to print them out
# We always do this and print out the top 150 sequencces - or until we run out of sequences
for i in range(1,index+1):  
    # if it's acceptable and we have written <= 150 so far write this one as well
    if reject[i] == "Acceptable" and total_aligns_written < 150:
        PREALIGN.write(title[i])
        PREALIGN.write("\n")
        PREALIGN.write(sequence[i])
        PREALIGN.write("\n")
        reject[i] = "Done"    # Mark this one as being DONE
        total_aligns_written += 1
    # If it's a reject write to the reject file and give a reason
    elif reject[i] != "Acceptable":    
        note = "REJECTED: reason is: " + str(reject[i])
        REJECT.write(title[i])
        REJECT.write("\n")
        REJECT.write(note)
        REJECT.write("\n")
        REJECT.write(sequence[i])
        REJECT.write("\n")
    # If we have written >150 and this one is also acceptable mark it for later consideration
    elif total_aligns_written >= 150:
        reject[i] = "Reconsider"
        reconsider += 1
    else: # Anything else gets rejected - really should not see this though
        note = "REJECTED: not sure why:" + str(reject[i])
        REJECT.write(title[i])
        REJECT.write("\n")
        REJECT.write(note)
        REJECT.write("\n")
        REJECT.write(sequence[i])
        REJECT.write("\n")


# At this point if asked we will randomly choose 150 from the reconsider group, unless there are < 150 of them
# in which case we'll take them all

if METHOD == "300": # Only use this if reconsider number is > 0
   if reconsider > 0 and reconsider <= 150:  # less than 150, but over 0 loop over the entire set looking for "R" and use them all
      print("Using all")
      for i in range(1,index+1):
           if reject[i] == "Reconsider":
               PREALIGN.write(title[i])
               PREALIGN.write("\n")
               PREALIGN.write(sequence[i])
               PREALIGN.write("\n")
               reject[i] = "Done"  # Mark this one as DONE

# If there are more than 150, then we will pick ones at random, check if they are Done, Reconsider or rejectss (a number != 0)
# if reconsider we use them and then set them to done "Done"
   elif reconsider > 150 :
       print("Reconsidering")
       reconsider_loop = 1
       while reconsider_loop <= 150:
            picked = random.randint(1,index+1)
            if reject[picked] == "Reconsider":
               PREALIGN.write(title[picked])
               PREALIGN.write("\n")
               PREALIGN.write(sequence[picked])
               PREALIGN.write("\n")
               reject[picked] = "Done"  # Mark this one as DONE
               reconsider_loop += 1

PREALIGN.close()
REJECT.close()
