from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO
from Bio.Align import AlignInfo
import os, shutil, traceback,subprocess
import numpy as np
import statistics as st
import time

# Config file #
from config import *

# Private defs
import dot
from blastnTom5 import *

# Call subprocess
def exec_cmd(comand,out,message=""):
    
    try:
        out = open(out,"w")
    except ValueError:
        if out == -3:
            pass
        else:
            print(error)
            sys.exit(-1)

    list_of_args = comand.split()

    try:
        subprocess.check_call(list_of_args,stdout=out)
        print(message)
    except ImportError:
        print("unable to import subprocess")
        sys.exit(-1)

    except OSError:
        print("1) may be some executive file is missing in %s. 2) Check program paths." % str(list_of_args))
        sys.exit(-1)

    except subprocess.CalledProcessError:
        print("error in executing the commands")
        sys.exit(-1)


nucs = "A,T,C,G"

def seq_inversor(seq):
    seq_inv = ""
    for nuc in seq:
        if nuc == "A":
            seq_inv = seq_inv + "T"
        
        if nuc == "T":
            seq_inv = seq_inv + "A"
        
        if nuc == "C":
            seq_inv = seq_inv + "G"
        
        if nuc == "G":
            seq_inv = seq_inv + "C"
    
    return(seq_inv[::-1])


def compare(seq,seq_test,count,Miss_pair):
    coordy = []
    for len_seq in range(len(seq)):
        miss = 0
        
        try:
            for len_seq_test in range(len(seq_test)):
                if seq[len_seq + len_seq_test] != seq_test[len_seq_test]:
                    miss +=1

                if miss/int(len(seq_test)) > Miss_pair:
                    break
                    

        except IndexError:            
            break

        if miss/int(len(seq_test)) <= Miss_pair:
            coordy.append(len_seq)
            count += 1
    
    return(coordy)


# return focus points in a key dictionary
def low_val(lenghts):
    list_del = []
    for i in lenghts.keys():
        for j in lenghts.keys():
            if abs(i-j) < 100 and abs(i-j) > 0:
                if lenghts[i][1] < lenghts[j][1]:
                    list_del.append(j)
                else:
                    list_del.append(i)


    for del_ in list_del:
        try:
            lenghts.pop(del_) 
        except:
            pass

    return lenghts

# Points out line 
def out_line(list_lenghts):
    list_pop = []
    for a in list_lenghts:
        c = 0
        for i in range(a-20,a+20):
            c += list_lenghts.count(i)
            if c > 5:
                break

        if c < 5:
            try:
                list_pop.append(a)
            except:
                pass

    return(list_pop)


# Pegar os pts na origem da seq original + create query file
def origin_pt(seq,outfile,lenghts_inv):
    try:
        align = AlignIO.read(outfile, "clustal")
    except: 
        shutil.rmtree(output_dir)
        return(print(""))
    

    dic_pos_gaps = {}

    for seq in align:
        gaps_count = 0
        for nuc in seq:
            if nuc == "-":
                gaps_count += 1 

            if nuc in nucs:
                break
        
        dic_pos_gaps[seq.description] = gaps_count

    jump = (max(dic_pos_gaps.values()))

    dic_pts = []
    dic_pts_inv = []
    for seq in align: 
        nuc_count = 0
        nuc_gaps_count = 0
        for nuc in seq:
            nuc_gaps_count += 1
            if nuc in nucs:
                nuc_count += 1

            if nuc_gaps_count == jump:
                break
        
        if int(seq.description)-20 in lenghts_inv:
            dic_pts_inv.append(int(seq.description)-nuc_count)
        
        if int(seq.description)-20 not in lenghts_inv:
            dic_pts.append(int(seq.description)+nuc_count)
            
    return(sorted(dic_pts),sorted(dic_pts_inv))

# Consensus
def cons_func(my_pssm):
    consensus = []
    for line in my_pssm:
        list_maf = []
        nuc_mf = max(line.values()) # Max Frequency
        
        for item in line:
            if line[item] == nuc_mf:
                list_maf.append(item)
        
        if len(list_maf)==1:
            if "-" in list_maf:
                consensus.append("")
            else:
                consensus.append(list_maf[0])
        else:
            if "-" in list_maf:
                list_maf.remove("-")
            
            consensus.append(list_maf[0])
            
    return("".join(consensus))

# Nanopolish steps
def racon(reads, output_dir, output, single_consensus):
    # racon
    os.system("sudo %s -a -x map-ont %s %s > %s/minimap.sam"%(minimap_path,single_consensus,reads,output_dir)) 
    os.system("sudo %s view -bS %s/minimap.sam > %s/minimap.bam"%(samtools_path,output_dir,output_dir))
    os.system("sudo %s sort %s/minimap.bam -o %s/sorted.bam"%(samtools_path,output_dir,output_dir))
    os.system("sudo %s index %s/sorted.bam"%(samtools_path,output_dir))
    os.system("sudo %s %s/msa.fasta %s/minimap.sam %s/consensus.fasta > %s/consensus_racon.fasta"%(racon_path,output_dir,output_dir,output_dir,output_dir))   
    

def alter_readdb(output_dir,reads):
    readdb_orig, readdb_new = "%s.index.readdb"%(reads), "%s.index.readdb_"%(output_dir)
    output_write = open(readdb_new,"w")
    for line in open(readdb_orig,"r"):
        save1 = line.split()[0]
        try:
            save2 = line.split()[1]
        except IndexError:
            output_write.write("%s\t%s\n"%(save1,save2))
            continue
            
        output_write.write(line)
    
    output_write.close()
    
    os.remove(readdb_orig)
    os.rename(readdb_new,readdb_orig)

    
def strat_align(rec,Miss_pair,output,outfile_consensus,min_len,max_len,min_reps):
    start = time.time()
    output_dir = "%s/%s"%(output,rec.id)    
    single_consensus = "%s/consensus.fasta"%(output_dir)
    
    seq = str(rec.seq)
    
    id_f = rec.id
    # Output directory
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass

    print("Arquivo de saída: %s"%(output_dir))

    # Get seq from File
    seq = seq[50:]

    seqy = seq[:800]

    # Compare seq_test with Seq
    lenghts,lenghts_inv,list_lenghts,list_lenghts_inv = dot.dot_(seq,seqy,Miss_pair)
    
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


    # lowest value from each line (dotplot)
    lenghts = low_val(lenghts)
    lenghts_inv = low_val(lenghts_inv)

    lenghts.update(lenghts_inv)

    # Align transition sequences
    msa_transition = open("%s/msa_transition.txt"%(output_dir),"w")

    tran_seq = []
    tran_start = 20
    tran_end = 300
    
    for tran_point in lenghts.keys(): # P/ cada pt de transição
        try:
            seq[tran_point]
        except:
            continue

        if tran_point == 0:
            msa_transition.write(">%s\n%s\n\n"%(tran_point,seq[tran_point:tran_point+tran_end]))

        if tran_point != 0 and tran_point not in lenghts_inv:
            msa_transition.write(">%s\n%s\n\n"%(tran_point-tran_start,seq[tran_point-tran_start:tran_point+tran_end]))
        
        if tran_point in lenghts_inv:
            msa_transition.write(">%s\n%s\n\n"%(tran_point+tran_start,seq_inversor(seq[tran_point-tran_end:tran_point+tran_start])))
             
        if tran_point in range(len(seq)-tran_start,len(seq)) and tran_point in lenghts_inv:
            msa_transition.write(">%s\n%s\n\n"%(tran_point,seq_inversor(seq[tran_point-tran_end:tran_point])))
        
    msa_transition.close()

    in_file = "%s/msa_transition.txt"%(output_dir)
    outfile = "%s/msa_transition.aln"%(output_dir)

    clustalw_cline = ClustalwCommandline("clustalw",infile = in_file,outfile = outfile)
    
    try:
        clustalw_cline()
    except:
        shutil.rmtree(output_dir)
        return(print("few sequences for alignment"))


    # pts in the original sequence
    try: 
        origin_pts, origin_pts_inv = origin_pt(seq,outfile,lenghts_inv)
    except:
        shutil.rmtree(output_dir)
        return(print("few sequences for alignment"))

    origin_pts_tot = origin_pts + origin_pts_inv
    print("Normal: %s\nInvertida: %s"%(origin_pts, origin_pts_inv))

    final_align_txt = open("%s/msa.fasta"%(output_dir),"w")
    reads = ("%s/msa.fasta"%(output_dir))

    # Genome linarization + Repetitions add to file contigs
    test_list, seqs_lenght = [], []
    num = 0
    rec_id_save = rec.id
    
    query = open("%s/query.fasta"%(output_dir),"w")
    bool_ = True
    for pt in sorted(origin_pts_tot):
        test_list.append(pt)
        if len(test_list) == 3:
            test_list.pop(0)
        
        if len(test_list) == 2:
            rec.description = ""
            rec.id = rec_id_save
            
            if test_list[0] in origin_pts and test_list[1] in origin_pts:
                rec.id = "%s_%s"%(rec.id,num)

                if rec.seq != "":
                    try:
                        SeqIO.write(rec[test_list[0]+50:test_list[1]+50],final_align_txt,"fasta")
                        if bool_ == True:
                            SeqIO.write(rec[test_list[0]+50:test_list[1]+50],query,"fasta")
                            bool_ = False
                    except:
                        pass

                seqs_lenght.append(len(seq[test_list[0]:test_list[1]]))
                num += 1
                
 
            if test_list[0] in origin_pts_inv and test_list[1] in origin_pts_inv:
                reverse = rec.reverse_complement()
                reverse.description = ""
                reverse.id = "%s%s"%(rec.id,num)
                
                if reverse.seq != "":
                    try: 
                        SeqIO.write(reverse[test_list[0]+50:test_list[1]+50],final_align_txt,"fasta") 
                        if bool_ == True:
                            SeqIO.write(rec[test_list[0]+50:test_list[1]+50],query,"fasta")
                            bool_ = False
                    except:
                        pass

                seqs_lenght.append(len(seq[test_list[0]:test_list[1]]))
                num += 1 
                

            else:
                pass # Jogar em outro arquivo depois
    query.close()

    rec.id = rec_id_save
    # Finish program if low sequences
    if num < min_reps:
        shutil.rmtree(output_dir)
        return(print("< min repetitons"))
    
    # Compare lenghts
    bool_ = True
    dif_len = False
    for rep_l in seqs_lenght:
        if bool_:
            rep_0 = rep_l
            bool_ = False
            continue

        max_dif = 300
        if rep_l-max_dif < 0:
            max_dif = 100
        if rep_l+max_dif < rep_0 or rep_l-max_dif > rep_0:
            dif_len = True 

    if st.mean(seqs_lenght) < float(min_len) or st.mean(seqs_lenght) > float(max_len) or dif_len == True:
        shutil.rmtree(output_dir)
        return(print("Bad lenght"))


    final_align_txt.close()
    
    # Blastn alignment
    os.system("makeblastdb -in %s/msa.fasta -dbtype nucl -parse_seqids"%(output_dir))
    os.system("blastn -db %s/msa.fasta -query %s/query.fasta -outfmt '6 sseqid sstart send slen qstart qend qlen evalue length nident mismatch gaps sseq qseq qseqid' -out %s/blastn.aln"%(output_dir,output_dir,output_dir))
    
    # Converter blastn to M5
    input_conv = open("%s/blastn.aln"%(output_dir),"r")
    output_conv = open("%s/blastn.m5"%(output_dir),"w")
    conv_B_M5(input_conv,output_conv)

    #PBDAGON Consensus
    exec_cmd("pbdagcon -t 0 -c 1 -m 1 %s/blastn.m5"%(output_dir), "%s/consensus.fasta"%(output_dir))

    for rec in SeqIO.parse("%s/consensus.fasta"%(output_dir),"fasta"):
        outfile_consensus.write(">%s\n%s\n"%(id_f,rec.seq)) 
    
    return(len(seq),True,time.time()-start,num)

