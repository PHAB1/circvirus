# circvirus

# Description

# Requeriments

- Python 3.0 or >

- MUSCLE

- Clustalx

- Biopython

- PBDAGCON

- racon and requeriments 

- Kaiju and requeriments 

- Diamond and requeriments 

- Medaka and requeriments

- Blastn


# Use

Genome assemble:
The output contains all genomes in the sought specifications, you may have to make changes to the config file (pathway for minimap2, samtools and racon programs).
 
	First use: 
	
		1. cd ~/circvirus/scripts/
		2. python3 setup.py build_ext --inplace


	python3 ~/circvirus/scripts/genomeassemble.py -i <sequences> -o <output> ## Lembrar de ver se o fast5 será necessário ##
	
		<sequences>
			-i --input
			input file in FASTA/FASTQ 
		
		<output>
			-o --output
			output directory

	options:
		-m --miss_pair => Default = 0.22
		The error rate for comparing k-mers 
		
		-mx --max_len => Default = 3000
		Max genome lenght

		-mn --min_len => Default = 1200
		Min genome length				 

		-r --repetitions => Default = 3
		Min number of repetitions. as higher as the number of repeats, the better the quality of the resulting genome.


Genome match: Genome classification using Kaiju and Diamond. you may have to make changes to the config file (pathway for Kaiju input and Diamond reference database). To set up a reference database for DIAMOND see https://github.com/bbuchfink/diamond.

	python3 ~/circvirus/scripts/virusmatch.py -i <sequences> -o <output>
		
		<sequences>
			-i --input 
			input file in FASTA/FASTQ format
		<output>
			-o --output
			output directory and files (you can choose to add a name to the files to differentiate them. E.g. ~/dir/run1)










