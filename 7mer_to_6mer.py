import argparse
import operator
import math



arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("inputfile")

args=arg_parser.parse_args()

kmers={}
amounts={}

a=0

with open(args.inputfile) as f:
    next(f)
    for line in f:
        sequence=line.split()[0]
        #print sequence
        energyA=float(line.split()[1])
        energyB=float(line.split()[2])
        kmer1=sequence[0:-1]
        kmer2=sequence[1:]
        if kmer1 in kmers:
            kmers[kmer1]+=energyA + energyB
            amounts[kmer1]+=2.0
        else:
            kmers[kmer1]=energyA + energyB
            amounts[kmer1]=2.0
        if kmer2 in kmers:
            kmers[kmer2]+=energyA + energyB
            amounts[kmer2]+=2.0
        else:
            kmers[kmer2]=energyA + energyB
            amounts[kmer2]=2.0

for kmer in kmers:
    #print kmer  +"\t" + str(kmers[kmer]) + "\t" + str(amounts[kmer]) + "\t" + str(kmers[kmer]/amounts[kmer])
    print kmer  +"\t" + str(kmers[kmer]/amounts[kmer])

