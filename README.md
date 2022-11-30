# NGS_pipeline

To use the NGS pipeline run following commands

1: Go to the erda service site (https://erda.dk/wsgi-bin/jupyter.py) and start MODI, select the SLURM and open the terminal
2: Go to the modi_mount directory and clone my github repository: cd modi_mount && git clone https://github.com/KarlDouglas/NGS_pipeline.git
3: Give permission to run the two shell scripts: chmod +x run.sh && chmod +x run2.sh
4: Run the script: sbatch run.sh
thats it! The script will make a directory named data where the result will be deposited

note: The script takes a while to run on a full dataset. if you want to test the pipeline oyu can copy a subset of the data eg. 100.000 lines using command: sed -n '1,100000p' "infile" > "outfile"

