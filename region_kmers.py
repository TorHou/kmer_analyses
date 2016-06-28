import argparse
import operator
import math
from Bio import SeqIO



arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("inputfile")
arg_parser.add_argument("bnsfile")
arg_parser.add_argument("--k", type=int, default=6)
arg_parser.add_argument("--window", type=int, default=2)
arg_parser.add_argument("--bns_column", type=int, default=5)

args=arg_parser.parse_args()
k = args.k
window = args.window

threshold = 0
kmer_affinities={}
kmers_many={}
filetype = ""

with open(args.bnsfile) as f:
    for line in f:
        if not line.startswith("["):
            kmer_affinities[line.split()[0]]=float(line.split()[args.bns_column])


fasta_sequences = SeqIO.parse(open(args.inputfile),'fasta')
for fasta in fasta_sequences:
    sequence = str(fasta.seq).upper()
    lstart = len(sequence)-k
    max_affinity = -1000
    affinity_window_max = 0
    affinity_sum = 0
    affinity_gm = 1.0
    for i in range(0,lstart+1):
        kmer=sequence[i:i+k]
        if kmer_affinities[kmer] > max_affinity:
            max_affinity = kmer_affinities[kmer]
        affinity_sum = affinity_sum + kmer_affinities[kmer]
        affinity_gm = float(affinity_gm) * float(kmer_affinities[kmer])
        affinity_window = 0
        if i+k+window-1<=len(sequence):
            for j in range(0,window):
                kmer=sequence[i+j:i+k+j]
                affinity_window = affinity_window + kmer_affinities[kmer]
        if affinity_window > affinity_window_max:
            affinity_window_max = affinity_window
    
            
    
    affinity_sum_norm = affinity_sum / len(sequence)
    affinity_gm_norm= math.pow(affinity_gm,1.0/float(len(sequence)))
    affinity_window_max_norm= affinity_window_max / window

    print fasta.id[2:] +"\t"+str(max_affinity)+"\t"+str(affinity_window_max)+"\t"+str(affinity_sum_norm)+"\t"+str(affinity_gm_norm)


