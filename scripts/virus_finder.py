import os, argparse
from Bio import SeqIO

from config import *
from defs3 import *

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',type=str,required=True)
parser.add_argument('-o','--output',type=str,required=True)

args = parser.parse_args()
output = args.output

## Kaiju files ##
kaiju_input_file = args.input
seq_tax_known = open("%sknown.fasta"%(output),"w")
seq_tax_unknown = open("%sunknown.fasta"%(output),"w") 
kaiju_output_file = "%skaiju.out"%(output)

## Diamond files ##
diamond_input_file = "%sknown.fasta"%(output)
diamond_output_file = "%smatches.m8"%(output)
diamond_table = "%sdiamond.out"%(output)
diamond_table_tax = "%sdiamond.names.out"%(output)

## Kaiju Def ##
def kaiju_match(kaiju_db_dir,kaiju_input_file,kaiju_output_file):
    os.system("sudo %s -t %s/nodes.dmp -f %s -i %s -o %s"%(kaiju_path,kaiju_db_dir,kaijudb,kaiju_input_file,kaiju_output_file))

    dict_matchs = {}
    # Create dic[Id_seq] = [seq,tax_idKaiju,tax_idDiamond]
    for rec in SeqIO.parse(kaiju_input_file,"fasta"):
        rec.id = "".join((rec.id).split("/")[0])
        dict_matchs[rec.id] = [str(rec.seq)]

    for line in open(kaiju_output_file,"r"):
        line = line.split()
        seq_tax_known.write(">%s\n%s\n"%(line[1],dict_matchs[line[1]][0])) if line[0] == "C" else seq_tax_unknown.write(">%s\n%s\n"%(line[1],dict_matchs[line[1]][0]))
    seq_tax_known.close()

    '''
    # Análise gráfica
    os.system("sudo ~/kaiju/bin/kaiju2krona -t %s/nodes.dmp -n %s/names.dmp -i %s -o %skaiju.out.krona"%(kaiju_db_dir,kaiju_db_dir,kaiju_output_file,output))
    os.system("sudo ktImportText -o %skaiju.out.html %skaiju.out.krona"%(output,output))

    # Taxon
    os.system("sudo ~/kaiju/bin/kaiju-addTaxonNames -t %s/nodes.dmp -n %s/names.dmp -i %s -o %skaiju.names.out"%(kaiju_db_dir,kaiju_db_dir,kaiju_output_file, output))
    '''
    print("----------------------------------------------------------------")
    return(dict_matchs)

## Diamond Def ##
def diamond_match(diamond_db,diamond_input_file,diamond_output_file):
    os.system("sudo %s blastx -d %s -q %s -f 102 -o %s"%(diamond_path,diamond_db,diamond_input_file,diamond_output_file))

## Run ##
dict_matchs = kaiju_match(kaiju_db_dir,kaiju_input_file,kaiju_output_file)
diamond_match(diamond_db,diamond_input_file,diamond_output_file)


# Construct table
for line in open(diamond_output_file,"r"):
    line = line.split()
    dict_matchs[line[0]].append(line[1])

diamond_table_w = open(diamond_table,"w")
for item in dict_matchs.items():
    try:
        if item[1][1] != "0":
            diamond_table_w.write("C\t%s\t%s\n"%(item[0],item[1][1]))
        if item[1][1] == "0":
            diamond_table_w.write("U\t%s\t%s\n"%(item[0],0))
    except IndexError:
        pass

diamond_table_w.close()

# Análise gráfica
os.system("sudo ~/kaiju/bin/kaiju2krona -t %s/nodes.dmp -n %s/names.dmp -i %s -o %sdiamond.out.krona"%(kaiju_db_dir,kaiju_db_dir,diamond_table,output))
os.system("sudo ktImportText -o %sdiamond.out.html %sdiamond.out.krona"%(output,output))

# Taxonomia
os.system("sudo ~/kaiju/bin/kaiju-addTaxonNames -t %s/nodes.dmp -n %s/names.dmp -i %s -o %sdiamond.names.out"%(kaiju_db_dir,kaiju_db_dir,diamond_table, output))

seq_tax_known = open("%sknown_2.fasta"%(output),"w")
for line in open(diamond_table_tax,"r"):
    line = line.split()
    seq_tax_known.write(">%s|%s\n%s\n"%(line[1]," ".join(line[3:]),dict_matchs[line[1]][0])) if line[0] == "C" else seq_tax_unknown.write(">%s\n%s\n"%(line[1],dict_matchs[line[1]][0]))

seq_tax_known.close()
seq_tax_unknown.close()
