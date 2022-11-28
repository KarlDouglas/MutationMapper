import main
mask1 = main.mask("fastq_demofile.txt")
mask2 = main.mask("fastq_demofile2.txt")
sort = main.sort_barcodes(mask1,mask2)
main.write_file(sort, "file2")
count = main.count_nucleotides("file1")
#depth = main.count_depth("file1")
#mutations = main.calculate_mutations(count, depth)
#print(mutations)