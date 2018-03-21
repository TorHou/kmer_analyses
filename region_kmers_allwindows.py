import argparse
import operator
import math
import numpy as np
from Bio import SeqIO



arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("inputfile")
arg_parser.add_argument("bnsfile")
arg_parser.add_argument("--k", type=int, default=6)
#arg_parser.add_argument("--window", type=int, default=1)
arg_parser.add_argument("--bns_column", type=int, default=5)

args=arg_parser.parse_args()
k = args.k
#window = args.window

threshold = 0
kmer_affinities={}
kmers_many={}
filetype = ""

with open(args.bnsfile) as f:
    for line in f:
        if not line.startswith("["):
            kmer_affinities[line.split()[0]]=float(line.split()[args.bns_column])

maxwindow = 100
minwindow = 1


fasta_sequences = SeqIO.parse(open(args.inputfile),'fasta')
for fasta in fasta_sequences:
    sequence = str(fasta.seq).upper()
    lstart = len(sequence)-k
    #max_affinity = -1000
    affinity_window_max = 0
    affinity_sum = 0
    affinity_gm = 1.0
    affinity_window_max = np.array([0.0]*maxwindow)
    for i in range(0,lstart+1):
        kmer=sequence[i:i+k]
        #if kmer_affinities[kmer] > max_affinity:
        #    max_affinity = kmer_affinities[kmer]
        affinity_sum = affinity_sum + kmer_affinities[kmer]
        affinity_gm = float(affinity_gm) * float(kmer_affinities[kmer])
        affinity_window = np.arange(minwindow,maxwindow+1)
        for idx, window in enumerate(range(minwindow,maxwindow+1)):
            if i+k+window-1<=len(sequence):
                elements = np.array([0.0]*window)
                for j in range(0,window):
                    kmer=sequence[i+j:i+k+j]
                    elements[j] = kmer_affinities[kmer]
                    #affinity_window[window] = affinity_window[window] * kmer_affinities[kmer]
                product = math.pow(np.prod(elements),1.0/float(window))
                if affinity_window_max[idx] < product: 
                    affinity_window_max[idx] = product
           
            #if affinity_window > affinity_window_max:
            #    affinity_window_max = affinity_window
    
            
    
    affinity_sum_norm = affinity_sum / len(sequence)
    if len(sequence)-k+1 > 0 :
        affinity_gm_norm= math.pow(affinity_gm,1.0/float(len(sequence)-k+1))
    #affinity_window_max_norm= affinity_window_max / window
    window_energies = ""
    #for window in range(minwindow, maxwindow+1):
    #    window_energies += str(affinity_window_max[window]) + "\t"
    #window_energies = "\t".join(`affinity` for affinity in  affinity_window_max)
    window_energies = "\t".join(map(str,affinity_window_max))
        

    print fasta.id +"\t"+window_energies+"\t"+str(affinity_sum_norm)+"\t"+str(affinity_gm_norm)


