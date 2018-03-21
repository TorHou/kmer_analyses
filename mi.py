import argparse


arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("inputfile")
arg_parser.add_argument("--k", type=int, default=7)

args=arg_parser.parse_args()
k = args.k

threshold = 0
kmers={}
kmers_many={}
filetype = ""

with open(args.inputfile) as infile:
    for line in infile:
        if line.startswith('>'): 
            filetype = "fasta"
            break
        else:
            filetype = "sam"
            

with open(args.inputfile) as infile:
    if filetype == "fasta":
        whole_line = ""
        for line in infile:
            line = line.rstrip().upper()
            if not line.startswith('>'):
                whole_line = whole_line + line
            else:
                lstart = len(whole_line)-k
                for i in range(0,lstart+1):
                    kmer=whole_line[i:i+k]
                    if kmer in kmers:
                        kmers[kmer] += 1
                    else:
                        kmers[kmer] = 1
                whole_line = ""
    elif filetype == "sam":
        for line in infile:
            split_line = line.split()
            seq = split_line[9]
            lstart = len(seq)-k
            for i in range(0,lstart+1):
                kmer=seq[i:i+k]
                if kmer in kmers:
                    kmers[kmer] += 1
                else:
                    kmers[kmer] = 0
            
    
                    
                
                
for kmer in kmers:             
    if kmers[kmer] > threshold:
        print kmer + "\t" + str(kmers[kmer])
            


