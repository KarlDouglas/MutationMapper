import main
mask1 = main.mask("fastq_demofile.txt")
mask2 = main.mask("fastq_demofile2.txt")
paired = main.fastq_pair(mask1,mask2)
sort = main.sort_barcodes(paired)
merged = main.merge(sort)


#merged = main.merge(paired)

