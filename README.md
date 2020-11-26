# circvirus

# Requeriments

- Python 3.0 or >

- Kaiju (virus_finder only)

- Diamond (virus_finder only)

- racon

- Medaka

clustalw pbdagcon build-essential cmake ncbi-blast+ clustalw gcc make libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev numpy statsmodels==0.9.0 HTSeq pandas pybedtools==0.7.10 cython muscle biopython==1.74 -> Can be installed using => python3 ~/circvirus/install.py (Requires pip3).


# Use

Concatamer Consensus:
The output contains all genomes in the sought specifications (in directories) and a file with all the consensus, you may have to make changes to the config file (pathway for minimap2, samtools and racon programs). => config file in ~/circvirus/scripts/config.py
 
	First use: 
	
		1. cd ~/circvirus/scripts/
		2. python3 setup.py build_ext --inplace


	python3 ~/circvirus/scripts/concatamer_consensus.py -i <sequences> -o <output>
	
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

Global Consensus:
Use concatamer consensus output file as input to generate Global consensus of all viruses found in the sample, the output contains the alignment of each consensus in genomes_corrected.fasta and consensus in consensus.fasta file within directories containing information from each cluster of sequences. additionally, global consensus is available in the file genomes_linearized.fasta

	python3 ~/circvirus/scripts/global_consenus.py -i <sequences> -d <directory> -o <output>
	
		<sequences>
			-i --input
			Concatamer consensus output file in FASTA format 
		
		<directoy>
			-d --dir
			Concatamer consensus directory output
			
		<output>
			-o --output
			output directory
			
		options:
			-mc --min_cover => defalt = 50 
			Minimun coverage for global consensus
		
			-mi --min_id
			Minimun identity between concatamer consensus sequence and reads
		


Virus Finder: Genome classification using Kaiju and Diamond. you may have to make changes to the config file (pathway for Kaiju input and Diamond reference database (e.g. nr.dmnd)). To set up a reference database for DIAMOND see https://github.com/bbuchfink/diamond.

	python3 ~/circvirus/scripts/virus_finder.py -i <sequences> -o <output>
		
		<sequences>
			-i --input 
			input file in FASTA/FASTQ format
			
		<output>
			-o --output
			output directory and files (you can choose to add a name to the files to differentiate them. E.g. ~/dir/run1)










