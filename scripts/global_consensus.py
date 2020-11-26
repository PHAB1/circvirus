import os, sys, argparse, time
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Align import AlignInfo
from Bio.Seq import Seq
import numpy as np

#Config file
from config import *

import dot
from defs3 import *

time1 = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',type=str,required=True)
parser.add_argument('-d','--dir',type=str,required=True)
parser.add_argument('-o','--output_dir',type=str,required=True)
parser.add_argument('-mc','--min_cover',default=50)
parser.add_argument('-mi','--min_id',default=0.8)
args = parser.parse_args()

arch = args.input
con_dir = args.dir
output_dir = args.output_dir
min_cover = float(args.min_cover)
min_id = float(args.min_id)

# Creating output directory
try:
    os.mkdir(output_dir)
except:
    pass

genome_lin = open("%s/genomes_linearized.fasta"%(output_dir),"w")
stats = open("%s/stats"%(output_dir),"w")

# Joining all genomes + refseqs dictionary(Id:[seq,start,end])
seqs_dic = {}
i = 0

# Seqs_dic = {Id:[seq,start,end} 
for rec in SeqIO.parse(arch,"fasta"):
    seqs_dic[(rec.id).split("/")[0]] = [rec.seq,i,i+len(rec.seq)] 
    i += len(rec.seq)


count = 1
seq_lenght = 0
seq_lenght_sub = 0
keys_dic = [key for key in seqs_dic.keys()]
for key in keys_dic:
    time_iter = time.time()
    seq_len = []

    try:
        seqs_dic[key]
    except:
        continue
    
    #Update seqs_dictionary
    i = 0
    for rec in seqs_dic.items():
        len_seq_i = len(rec[1][0])

        rec[1][1] = i
        rec[1][2] = i+len_seq_i
        i += len_seq_i
 
    #Update seq
    seq = ""
    for item in seqs_dic.items():
        seq += "".join(item[1][0])

    seq_lenght = len(seqs_dic[key][0])
    seq_lenght_sub = len(seqs_dic[key][0])

    output_dir_count = "%s/%s"%(output_dir,count)

    try:
        os.mkdir("%s/%s"%(output_dir,count))
    except:
        pass
    
    #seqy
    seqy = str(seqs_dic[key][0][50:600])

    # Dotplot
    lenghts,lenghts_inv,list_lenghts,list_lenghts_inv = dot.dot_(seq,seqy,0)

    # removing interference #
    # Invert
    for point in out_line(list_lenghts_inv):
        try:
            lenghts_inv.pop(point)
        except:
            pass

    # Normal
    for point in out_line(list_lenghts):
        try:
            lenghts.pop(point)
        except:
            pass

    lenghts = low_val(lenghts)
    lenghts_inv = low_val(lenghts_inv)

    lenghts.update(lenghts_inv)

    f = open("%s/%s/genomes.fasta"%(output_dir,count),"w")
    t = open("%s/%s/msa_transition.fasta"%(output_dir,count),"w")
    q = open("%s/%s/query_transition.fasta"%(output_dir,count),"w")    
    
    seq_dic_est = {} #{Id:[corrected_seq,sseq,qseq,start,qstart]}
    
    #for pt in lenghts find correct seq and ~pt similar
    bool_ = False
    for pt in lenghts.items():

        for item in seqs_dic.items(): # Find correspondent coord.
            if item[1][2] > pt[1][0]:
                
                A = False
                for fas in SeqIO.parse("%s/%s/msa.fasta"%(con_dir,"".join(item[0][0:].split("/")[0])[0:]),"fasta"): #Concatamers fasta directory 0:-1 -> Correção no HPV 

                    max_dif = 500
                    if seq_lenght - max_dif < 0:
                        max_dif = 100

                    if len(fas.seq) > seq_lenght + max_dif or len(fas.seq) < seq_lenght - max_dif:
                        continue

                    #Inv Seqs
                    if pt[0] not in lenghts_inv:
                        est_seq = fas.seq[pt[0]-item[1][1]:]+fas.seq[0:pt[0]-item[1][1]]

                        # Update seqs_dic
                        A = True
                        seq_dic_est[fas.id] = [est_seq]
                        
                        if bool_ == False:
                            seq_dic_est["query"] = [Seq(str((item[1][0])[pt[0]-item[1][1]:]+(item[1][0])[0:pt[0]-item[1][1]]))]
                            bool_ = True
                    
                    # Normal Seqs
                    if pt[0] in lenghts_inv:
                        est_seq = (fas.seq[pt[0]-item[1][1]:]+fas.seq[0:pt[0]-item[1][1]]).reverse_complement()
                        
                        # Update seqs_dic_est
                        A = True
                        seq_dic_est[fas.id] = [est_seq]
                        
                
                # Update Seq_dictionary
                if A == True:
                    del seqs_dic[item[0]]

                break
    
    trans_lenght = 800
    for item in seq_dic_est.items():
        f.write(">%s\n%s\n"%(item[0],item[1][0]))
        t.write(">%s\n%s\n"%(item[0],item[1][0][0:trans_lenght])) 

        if item[0] == "query":
            q.write(">%s\n%s\n"%(item[0],item[1][0][0:trans_lenght])) 

    f.close()
    t.close()
    q.close()
   
    # Alignment -> Transitions Segments
    os.system("makeblastdb -in %s/msa_transition.fasta -dbtype nucl -parse_seqids"%(output_dir_count))
    os.system("blastn -db %s/msa_transition.fasta -query %s/query_transition.fasta -outfmt '6 sseqid sstart send slen qstart qend qlen evalue length nident mismatch gaps sseq qseq qseqid' -max_target_seqs 100000000 -out %s/blastn_transition.aln"%(output_dir_count,output_dir_count,output_dir_count))
    
    max_ref_len = []
    for curr_aln in open("%s/blastn_transition.aln"%(output_dir_count),"r"):
        curr_aln = curr_aln.split()
        seq_dic_est[curr_aln[0]].append(str(curr_aln[12]))
        seq_dic_est[curr_aln[0]].append(str(curr_aln[13]))
        seq_dic_est[curr_aln[0]].append(int(curr_aln[1]))
        seq_dic_est[curr_aln[0]].append(int(curr_aln[4]))

        max_ref_len.append(int(curr_aln[4]))

    max_ref_len = max(max_ref_len)
    
    reps_dic = {}
    query = open("%s/query.fasta"%(output_dir_count),"w")
    con_corr = open("%s/%s/genomes_corrected.fasta"%(output_dir,count),"w")
    bool_=True
    for item in seq_dic_est.items():
        # Blastn low score sequences correction
        try:
            item[1][0]
            item[1][1]
        except IndexError:
            continue

        # Find starting point with respect to the reference
        pos_s1=1
        pos_q1=1
        for s1,q1 in zip(item[1][1],item[1][2]):
            if s1 != "-":
                pos_s1 += 1
            
            if q1 != "-":
                pos_q1 += 1
            
            if pos_q1+int(item[1][4]) == max_ref_len:
                final_start_pt = pos_s1 + int(item[1][3]) - 1
                s1 = item[1][0]

                reps_dic[item[0]] = s1[final_start_pt:]+s1[:final_start_pt] 

                item[1][0] = s1[final_start_pt:]+s1[:final_start_pt]

                if item[0] == "query" and bool_:
                    query_len = float(len(s1[final_start_pt:]+s1[:final_start_pt]))
                    query.write(">%s\n%s\n"%(item[0],s1[final_start_pt:]+s1[:final_start_pt]))
                    bool_ = False
                
                break

    for rep in reps_dic.items():
        con_corr.write(">%s\n%s\n"%(rep[0],rep[1]))

    query.close()
    con_corr.close()
               

    '''
    # Alignment
    in_file = "%s/%s/msa_transition.fasta"%(output_dir,count)
    outfile = "%s/%s/msa_transition.aln"%(output_dir,count)
    
    clustalw_cline = MuscleCommandline(input = in_file,out = outfile,clw=False)
    clustalw_cline()

    c_align = AlignIO.read(outfile, "fasta")

    summary_align = AlignInfo.SummaryInfo(c_align)
    consensus = summary_align.dumb_consensus()
    my_pssm = summary_align.pos_specific_score_matrix(consensus)

    # Lower "-" in alignment
    pos_list = [0,1000]
    pos = 0
    for col in my_pssm:
        if col["-"] < pos_list[1]:
            pos_list[0] = pos
            pos_list[1] = col["-"]
        
        pos+=1

    # Acurate start position in genomes.fasta
    seq_start_pos_dic = {}
    for fas in c_align:
        real_pos = 0
        for pos in range(len(fas.seq)): 
            if pos == pos_list[0]:
                seq_start_pos_dic[fas.id] = real_pos
                break

            if fas.seq[pos] != "-":
                real_pos += 1

    # Correct sequences 
    cover = 0
    con_corr = open("%s/%s/genomes_corrected.fasta"%(output_dir,count),"w")
    query = open("%s/query.fasta"%(output_dir_count),"w")
    bool_=True
    for fas in SeqIO.parse("%s/%s/genomes.fasta"%(output_dir,count),"fasta"):
        cover += 1
        try:
            con_corr.write(">%s\n%s\n"%(fas.id,fas.seq[seq_start_pos_dic[fas.id]:]+fas.seq[:seq_start_pos_dic[fas.id]]))
        except:
            continue
        if bool_:
            query.write(">%s\n%s\n"%(fas.id,fas.seq[seq_start_pos_dic[fas.id]:]+fas.seq[:seq_start_pos_dic[fas.id]]))
            bool_ = False

    query.close()
    '''

    #Blastn alignment
    os.system("makeblastdb -in %s/genomes_corrected.fasta -dbtype nucl -parse_seqids"%(output_dir_count))
    os.system("blastn -db %s/genomes_corrected.fasta -query %s/query.fasta -outfmt '6 sseqid sstart send slen qstart qend qlen evalue length nident mismatch gaps sseq qseq qseqid' -max_target_seqs 100000 -out %s/blastn.aln"%(output_dir_count,output_dir_count,output_dir_count))
    
    reps_dic = {}
    con_corr = open("%s/%s/genomes_corrected.fasta"%(output_dir,count),"w")
    for curr_align in open("%s/blastn.aln"%(output_dir_count),"r"):
        curr_align = curr_align.split()
        seq_id = curr_align[0]
        nident = float(curr_align[9])
        sstart = int(curr_align[1])
        send = int(curr_align[2])
        
        if send-sstart < 0:
            nident = 0

        reps_dic[seq_id] = [seq_dic_est[seq_id][0],nident]
        

        seq_len.append(len(seq_dic_est[seq_id][0]))

    
    cover = 0
    for rep in reps_dic.items():

        if rep[1][1]/query_len < min_id:
            continue
            
        con_corr.write(">%s\n%s\n"%(rep[0],rep[1][0]))
        cover += 1


    con_corr.close()
    
    soma = sum(seq_len)
    dp = np.std(seq_len)

    if cover < min_cover:
        continue

    stats.write("\n----------------------\n")
    stats.write("%s\ncobertura: %s\nTempo_iter: %s\nMédia_len: %s\ndp: %s"%(count,cover,time.time()-time_iter,soma/cover,dp))

    '''
    # Verify lenght and %id
    blastn_aln_fileR = open("%s/blastn.aln"%(output_dir_count),"r")
    blastn_list = []
    for curr_align in blastn_aln_fileR:
        curr_align = curr_align.split()
        curr_align_pid = float(curr_align[9])/float(curr_align[6])
        len_dif = abs(float(curr_align[6])-float(curr_align[3]))

        if curr_align_pid > 0.9 and len_dif < 100.0:
            print(curr_align_pid, len_dif)
            blastn_list.append("%s\n"%(" ".join(curr_align)))

    
    blastn_aln_fileR.close()
    blastn_aln_fileW = open("%s/blastn.aln"%(output_dir_count),"w")
    
    for curr_align in blastn_list:
        blastn_aln_fileW.write(curr_align)

    blastn_aln_fileW.close()
    '''

    #Convert Blastn to M5
    input_conv = open("%s/blastn.aln"%(output_dir_count),"r")
    output_conv = open("%s/blastn.m5"%(output_dir_count),"w")
    conv_B_M5(input_conv,output_conv)

    #PBDAGON Consensus
    exec_cmd("pbdagcon -t 0 -c 1 -m 1 %s/blastn.m5"%(output_dir_count),"%s/consensus.fasta"%(output_dir_count))
    
    #Racon+medaka
    exec_cmd("%s -a -x map-ont %s/consensus.fasta %s/genomes_corrected.fasta"%(minimap_path,output_dir_count,output_dir_count) ,"%s/minimap.sam"%(output_dir_count))
    exec_cmd("%s view -bS %s/minimap.sam"%(samtools_path,output_dir_count), "%s/minimap.bam"%(output_dir_count))
    exec_cmd("%s sort %s/minimap.bam -o %s/sorted.bam"%(samtools_path,output_dir_count,output_dir_count),-3)
    exec_cmd("%s index %s/sorted.bam"%(samtools_path,output_dir_count),-3)
    exec_cmd("%s %s/genomes_corrected.fasta %s/minimap.sam %s/consensus.fasta"%(racon_path,output_dir_count,output_dir_count,output_dir_count), "%s/consensus_racon_%sx.fasta"%(output_dir_count,cover))
    exec_cmd("%s -i %s/genomes_corrected.fasta -d %s/consensus_racon_%sx.fasta -o %s"%(medaka_consensus_path,output_dir_count,output_dir_count,cover,output_dir_count),-3,"Global consensus %s Finished!"%(count))
   
    for seq in SeqIO.parse(open("%s/consensus.fasta"%(output_dir_count)),"fasta"):
        seq.id = count
        genome_lin.write(">%s\n%s\n"%(seq.id,seq.seq))

    count+= 1
    
    '''
    os.system("sudo %s -i %s/consensus_racons.fasta -d %s/genomes_corrected.fasta -o %s/consensus_medaka.fasta"%(medaka_consensus_path,output_dir_count,output_dir_count,output_dir_count))
    '''

    '''
    [plt.plot(val[0],val[1],"+",color="red") for val in lenghts.values()]
    [plt.plot(val[0],val[1],"+",color="blue") for val in lenghts_inv.values()]
    plt.show()
    '''
stats.close()
print(time.time()-time1) 
genome_lin.close()
