import os, sys, argparse, time
import numpy as np
from Bio import SeqIO

## Private import ##
from defs3 import *

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',type=str,required=True)
parser.add_argument('-o','--output_dir',type=str,required=True)
parser.add_argument('-m','--miss_pair',default=0.22) #0.22
parser.add_argument('-mx','--max_len',default=5000)
parser.add_argument('-mn','--min_len',default=1000)
parser.add_argument('-r','--repetitions',default=3)
args = parser.parse_args()

arch = args.input
ext = arch.split(".")[-1]

Miss_pair = args.miss_pair

# Diretórios com os Ids
output = args.output_dir

try:
    os.mkdir(output)
except:
    pass

# Arquivo de saída com os consensus
outfile_consensus = open("%s/consensus_%s"%(output,arch.split("/")[-1]),"w")

stats = open("%s/Stats_genomeAss.txt"%(output),"w")
stats.write("\n------------------------------------\n")

sum_seqs = 0
reads_count = 0
mean_time_read = 0
count_t = 0 #Contagem de reads totais
rep_count = {}

initt = time.time()
for rec in SeqIO.parse(arch,ext):
    count_t += 1
    rec.id = (rec.id).split()[0]
    rec.id = (rec.id).split("/")[0]

    if len(rec.seq) > int(args.min_len)*int(args.repetitions):
        print("-------------------------------------------------------------------------------")
        bool_r = False

        stat = strat_align(rec,Miss_pair,output,outfile_consensus,args.min_len,args.max_len,int(args.repetitions)) 

        if stat != None:
            seq_l = stat[0]
            bool_r = stat[1]
            time_read = stat[2]
            reps = stat[3]
            print("\n")
            
            if bool_r == True and not os.stat("%s/%s/consensus.fasta"%(output,rec.id)).st_size == 0:
                sum_seqs += seq_l
                reads_count += 1
                mean_time_read += time_read
                try:
                    rep_count[reps] += 1
                except:
                    rep_count[reps] = 1


outfile_consensus.close()

end_time = time.time()

stats.write("Mean_len: %s\nReads_count: %s\nmean_time_read: %s\nN_repetições: %s\nTotal_reads: %s\nT_time: %s"%(sum_seqs/reads_count,reads_count,mean_time_read/reads_count,rep_count,count_t,end_time-initt))
stats.close()

