import main
mask1 = main.mask("CY001_1.fq") # inputfile 1
mask2 = main.mask("CY001_2.fq") # inputfile 2
sort = main.sort_barcodes(mask1,mask2) # sorts the data by barcode
index = 0
for i in sort:
    index += 1
    main.write_file(i, ("CY001BC"+str(index))+"_sorted") # writes a txt file for each barcode to give to bowtie2
