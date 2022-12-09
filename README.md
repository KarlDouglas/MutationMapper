# NGS_pipeline

MutationMapper is a program designed to quantify the effect of targeted mutagenesis. Eg. error prone PCR or in-vivo approaches
- The program takes two FASTQ files as input (adaptor trimmed) and returns
    - A scatter plot illustrating the distribution of point mutations
    - A bar plot displaying the ration of bases that were mutated
    - A series of bar plots showing the ratio of what the mutated bases were substituted with 
    - A .txt with information used to make the plots (sequence depth, nucleotides counted, mutations found, mutation frequency, mutations per read(0-10 mutations))
- Bowtie2 is used to map the reads, the query template and additional settings can be adjusted in the run2.sh file


For Use on KU's erda service:
1: Go to the erda service site (https://erda.dk/wsgi-bin/jupyter.py) and start MODI, select the SLURM and open the terminal
2: Go to the modi_mount directory and clone my github repository: cd modi_mount && git clone https://github.com/KarlDouglas/MutationMapper.git
3: Give permission to run the two shell scripts: chmod +x run.sh && chmod +x run2.sh
4: Run the script: sbatch run.sh
thats it! The script will make a directory named data where the result will be deposited

note: The script takes a while to run on a full dataset. if you want to test the pipeline oyu can copy a subset of the data eg. 100.000 lines using command: sed -n '1,100000p' "infile" > "outfile"


Note on development.
If you want to tweak the parameters of the data analysis you mainly want to focus on the mian.py file as the functions are actually written here. The other python files simply calls on different functions from this file. 
