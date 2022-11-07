cd ..
fastq_masker -q 30 -r N -i /home/ptc872_ku_dk/modi_mount/unzip/X204SC22080130-Z01-F001/raw_data/CY002/CY002_FDDP220275354-1a_H7VYHDRX2_L1_1.fq -o masked1.fq
fastq_masker -q 30 -r N -i /home/ptc872_ku_dk/modi_mount/unzip/X204SC22080130-Z01-F001/raw_data/CY002/CY002_FDDP220275354-1a_H7VYHDRX2_L1_2.fq -o masked2.fq

cutadapt -b file:barcodes.fasta -m 20 -o {name}.1.fq masked1.fq
cutadapt -b file:barcodes.fasta -m 20 -o {name}.2.fq masked2.fq

fastq_pair BC1.1.fq BC1.2.fq
fastq_pair BC2.1.fq BC2.2.fq
fastq_pair BC3.1.fq BC3.2.fq
fastq_pair BC4.1.fq BC4.2.fq
fastq_pair BC5.1.fq BC5.2.fq
fastq_pair BC6.1.fq BC6.2.fq
fastq_pair BC7.1.fq BC7.2.fq
fastq_pair BC8.1.fq BC8.2.fq
fastq_pair BC9.1.fq BC9.2.fq
fastq_pair BC10.1.fq BC10.2.fq

/home/ptc872_ku_dk/modi_mount/NGmerge/NGmerge -1 BC1.1.fq.paired.fq -2 BC1.2.fq.paired.fq -o BC1.merged.fastq
/home/ptc872_ku_dk/modi_mount/NGmerge/NGmerge -1 BC2.1.fq.paired.fq -2 BC2.2.fq.paired.fq -o BC2.merged.fastq
/home/ptc872_ku_dk/modi_mount/NGmerge/NGmerge -1 BC3.1.fq.paired.fq -2 BC3.2.fq.paired.fq -o BC3.merged.fastq
/home/ptc872_ku_dk/modi_mount/NGmerge/NGmerge -1 BC4.1.fq.paired.fq -2 BC4.2.fq.paired.fq -o BC4.merged.fastq
/home/ptc872_ku_dk/modi_mount/NGmerge/NGmerge -1 BC5.1.fq.paired.fq -2 BC5.2.fq.paired.fq -o BC5.merged.fastq
/home/ptc872_ku_dk/modi_mount/NGmerge/NGmerge -1 BC6.1.fq.paired.fq -2 BC6.2.fq.paired.fq -o BC6.merged.fastq
/home/ptc872_ku_dk/modi_mount/NGmerge/NGmerge -1 BC7.1.fq.paired.fq -2 BC7.2.fq.paired.fq -o BC7.merged.fastq
/home/ptc872_ku_dk/modi_mount/NGmerge/NGmerge -1 BC8.1.fq.paired.fq -2 BC8.2.fq.paired.fq -o BC8.merged.fastq
/home/ptc872_ku_dk/modi_mount/NGmerge/NGmerge -1 BC9.1.fq.paired.fq -2 BC9.2.fq.paired.fq -o BC9.merged.fastq
/home/ptc872_ku_dk/modi_mount/NGmerge/NGmerge -1 BC10.1.fq.paired.fq -2 BC10.2.fq.paired.fq -o BC10.merged.fastq
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450  -l 16 --suppress 1,2,3,6,7,8,9 BC1.merged.fastq > BC1.M.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC1.1.fq.single.fq > BC1.F.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC1.2.fq.single.fq > BC1.R.map

/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC2.merged.fastq > BC2.M.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC2.1.fq.single.fq > BC2.F.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC2.2.fq.single.fq > BC2.R.map

/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC3.merged.fastq > BC3.M.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC3.1.fq.single.fq > BC3.F.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC3.2.fq.single.fq > BC3.R.map

/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC4.merged.fastq > BC4.M.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC4.1.fq.single.fq > BC4.F.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC4.2.fq.single.fq > BC4.R.map

/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC5.merged.fastq > BC5.M.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC5.1.fq.single.fq > BC5.F.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC5.2.fq.single.fq > BC5.R.map

/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC6.merged.fastq > BC6.M.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC6.1.fq.single.fq > BC6.F.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC6.2.fq.single.fq > BC6.R.map

/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC7.merged.fastq > BC7.M.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC7.1.fq.single.fq > BC7.F.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC7.2.fq.single.fq > BC7.R.map

/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC8.merged.fastq > BC8.M.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC8.1.fq.single.fq > BC8.F.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC8.2.fq.single.fq > BC8.R.map

/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3 -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC9.merged.fastq > BC9.M.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC9.1.fq.single.fq > BC9.F.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC9.2.fq.single.fq > BC9.R.map

/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3 -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC10.merged.fastq > BC10.M.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC10.1.fq.single.fq > BC10.F.map
/home/ptc872_ku_dk/modi_mount/bowtie/bowtie -x  /home/ptc872_ku_dk/modi_mount/bowtie/indexes/leu2TAG.ebwt -3 4 -5 4 -n 3  -e 450 -l 16 --suppress 1,2,3,6,7,8,9 BC10.2.fq.single.fq > BC10.R.map

mv BC1.F.map Seq_analysis
mv BC1.R.map Seq_analysis
mv BC1.M.map Seq_analysis
mv BC2.F.map Seq_analysis
mv BC2.R.map Seq_analysis
mv BC2.M.map Seq_analysis
mv BC3.F.map Seq_analysis
mv BC3.R.map Seq_analysis
mv BC3.M.map Seq_analysis
mv BC4.F.map Seq_analysis
mv BC4.R.map Seq_analysis
mv BC4.M.map Seq_analysis
mv BC5.F.map Seq_analysis
mv BC5.R.map Seq_analysis
mv BC5.M.map Seq_analysis
mv BC6.F.map Seq_analysis
mv BC6.R.map Seq_analysis
mv BC6.M.map Seq_analysis
mv BC7.F.map Seq_analysis
mv BC7.R.map Seq_analysis
mv BC7.M.map Seq_analysis
mv BC8.F.map Seq_analysis
mv BC8.R.map Seq_analysis
mv BC8.M.map Seq_analysis
mv BC9.F.map Seq_analysis
mv BC9.R.map Seq_analysis
mv BC9.M.map Seq_analysis
mv BC10.F.map Seq_analysis
mv BC10.R.map Seq_analysis
mv BC10.M.map Seq_analysis

cd Seq_analysis

python testrun.py