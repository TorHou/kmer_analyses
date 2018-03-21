import argparse
from Bio import SeqIO
import operator
from math import sqrt

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("inputfile")
arg_parser.add_argument("--k", type=int, default=7)
arg_parser.add_argument("--overlap", type=bool, default=False)

args=arg_parser.parse_args()
k = args.k

alphabet_size = 4
possible_kmers = alphabet_size**k

kmers={}

fasta_sequences = SeqIO.parse(open(args.inputfile),'fasta')
sumks = 0.0
positions = 0

for fasta in fasta_sequences:
    seq = fasta.seq.upper()
    #print seq 
    lstart = len(seq)-k
    for i in range(0,lstart+1):
        positions += 1
        kmer=seq[i:i+k]
        if str(kmer) in kmers:
            if args.overlap:
                [score, last] = kmers[str(kmer)]
                kmers[str(kmer)] = [score+1, positions]
                sumks += 1
            else:
                [score, last] = kmers[str(kmer)]
                if positions-last > k:
                    kmers[str(kmer)] = [score+1, positions]
                    sumks += 1
        else:
            kmers[str(kmer)] = [1, positions]
            sumks += 1
            
expectedRandom = sumks/possible_kmers
                
kmers_red ={}
for kmer in kmers:
    kmers_red[kmer] = kmers[kmer][0]
    
sorted_x = sorted(kmers_red.items(), key=operator.itemgetter(1), reverse=True)
print "Overall number of kmers: " + str(sumks)
#print "\n".join(str(i)[0]+": "+str(i/sumks) for i in sorted_x[0:10])
print "\n".join(str(i[0])+": "+str(i[1])+"\tStdRes: "+ str((i[1]-expectedRandom)/sqrt(expectedRandom)) for i in sorted_x[0:50])   


            


