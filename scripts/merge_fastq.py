import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',type=str,required=True)
parser.add_argument('-o','--output',type=str,required=True)
args = parser.parse_args()

paths = args.input.split(",")
output = args.output

output_file = open("%s"%(output),"w")

for path in paths:
    ext = path.split(".")[-1]
    for rec in SeqIO.parse(path,ext):
        output_file.write(str(">%s\n"%(rec.id)))
        output_file.write(str("%s\n"%(rec.seq)))

output_file.close()
