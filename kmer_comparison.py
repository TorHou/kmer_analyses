import argparse
from Bio import SeqIO
import operator
from math import sqrt
from scipy.stats import binom

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("inputfile1")
arg_parser.add_argument("inputfile2")
arg_parser.add_argument("--k", type=int, default=7)

args=arg_parser.parse_args()
k = args.k

alphabet_size = 4
possible_kmers = 4**k

kmers={}

fasta_sequences = SeqIO.parse(open(args.inputfile1),'fasta')
sumks = 0.0
for fasta in fasta_sequences:
    seq = fasta.seq.upper()
    #print seq 
    lstart = len(seq)-k
    for i in range(0,lstart+1):
        kmer=seq[i:i+k]
        sumks += 1
        if str(kmer) in kmers:
            kmers[str(kmer)] += 1
        else:
            kmers[str(kmer)] = 1
            
kmer_freq = {}
for kmer in kmers:
    kmer_freq[kmer] = kmers[kmer]/sumks
#print kmer_freq

fasta_sequences = SeqIO.parse(open(args.inputfile2),'fasta')
sumks2 = 0.0
kmers2={}
for fasta in fasta_sequences:
    seq = fasta.seq.upper()
    lstart = len(seq)-k
    for i in range(0,lstart+1):
        kmer=seq[i:i+k]
        sumks2 += 1
        if str(kmer) in kmers2:
            kmers2[str(kmer)] += 1
        else:
            kmers2[str(kmer)] = 1

kmer_sig = {}
for kmer in kmers2:
    if not kmer in kmers:
        kmers[kmer]=0.5
    fr2 = kmers2[kmer]/sumks2
    fr1 = kmers[kmer]/sumks

    kmer_sig[kmer] = [kmers[kmer],kmers2[kmer],fr1/fr2]
#print "binom: " + str(binom.cdf(3,30,0.1))

kmer_sig_tuples = [(k, v[0],v[1],v[2]) for k, v in kmer_sig.iteritems()]

kmer_sig_sorted = sorted(kmer_sig_tuples, key=(operator.itemgetter(3)), reverse=False)
print "counts1\tcounts2\tsignificance"
for v,kmer in enumerate(kmer_sig_sorted):
    if v < 40:
        print kmer
    
    
#expectedRandom = sumks/possible_kmers
                
#sorted_x = sorted(kmers.items(), key=operator.itemgetter(1), reverse=True)
#print "Overall number of kmers: " + str(sumks)
#print "\n".join(str(i)[0]+": "+str(i/sumks) for i in sorted_x[0:10])
#print "\n".join(str(i[0])+": "+str(i[1])+"\tStdRes: "+ str((i[1]-expectedRandom)/sqrt(expectedRandom)) for i in sorted_x[0:100])   


            


