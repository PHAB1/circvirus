import os, sys, argparse, time, shutil
import numpy as np
from Bio import SeqIO

#Config file
from config import *

parser = argparse.ArgumentParser()
parser.add_argument('-i','--reads',type=str,required=True)
parser.add_argument('-r','--reference',type=str,required=True)
parser.add_argument('-d','--directory',type=str,required=True)
parser.add_argument('-o','--output_dir',type=str,required=True)
parser.add_argument('-mr','--repetitions',default=3)
parser.add_argument('-mq','--min_qual',default=3)
parser.add_argument('-mf','--min_freq',default=0.01)
args = parser.parse_args()

reads = args.reads
reference = args.reference
directory = args.directory
output = args.output_dir
min_rep = int(args.repetitions)
min_qual = int(args.min_qual)
min_freq = int(args.min_freq)

dic = {}
for rec in SeqIO.parse(reads,"fasta"):
    if rec.id == "query":
        continue

    try:
        dic[(rec.id).split("_")[0]].append(rec.seq)
    except KeyError:
        dic[(rec.id).split("_")[0]] = [rec.seq]

dic_del_key = []
for item in dic.items():
    id = item[0] 
    seqs = item[1]

    file_ = open("%s/%s/msa2.fasta"%(directory,id),'w')
    try:
        shutil.rmtree("%s/%s/variant_vcf"%(directory,id))
    except NotADirectoryError:
        pass
    except FileNotFoundError:
        pass

    count = 0
    for seq in seqs:
        file_.write(">%s_%s\n%s\n"%(id,count,seq)) 
        count+=1

    if count < min_rep:
        dic_del_key.append(id)
        try:
            shutil.rmtree("%s/%s/variant_vcf"%(directory,id))
        except:
            pass

for id_ in dic_del_key:
    del dic[id_]


file_.close()

out_file = open("%s/vars"%(output),"w")
out_file_table = open("%s/vars_table"%(output),"w")

for item in dic.items():
    id = item[0]

    output_dir = ("%s/%s"%(directory,id))
    reads = ("%s/msa2.fasta"%(output_dir))
    os.system("%s -a -x map-ont %s %s > %s/minimapv.sam"%(minimap_path,reference,reads,output_dir))
    os.system("%s view -bS %s/minimapv.sam > %s/minimapv.bam"%(samtools_path,output_dir,output_dir))
    os.system("%s sort %s/minimapv.bam -o %s/sortedv.bam"%(samtools_path,output_dir,output_dir))
    os.system("%s index %s/sortedv.bam"%(samtools_path,output_dir))
    
    bam = ("%s/sortedv.bam"%(output_dir))
    os.system("%s -f %s -i %s -d -o %s/variant_vcf"%(medaka_var_path,reference,bam,output_dir))

    

    for line in open("%s/%s/%s"%(output_dir,"variant_vcf","round_1.vcf"),"r"):
        out_file.write(line)

out_file.close()

total_org = 0
pos = {}
for line in open("%s/vars"%(output),"r"):
    if line[0:4] == "##CL":
        id = "".join("".join(line.split(";")[0:-1]).split("/")[-2]) 

    if line[0:6] == "CHROM":
        total_org += 1
    
    if line[0:5] == "query":
        line = line.split()
        if line[6] == "PASS" and float(line[5]) > min_qual:
            p = line[1]
            m0 = line[3]
            m1 = line[4]

            try:
                pos[m0,p,m1].append(id)
            except KeyError:
                pos[m0,p,m1] = [id]

for item in pos.items():
    out_file_table.write("%s%s%s\t%s\n\n"%(item[0][0],item[0][1],item[0][2],item[1]))















