import main
mask1 = main.mask("CY001_1.fq")
mask2 = main.mask("CY001_2.fq")
sort = main.sort_barcodes(mask1,mask2)
index = 0
for i in sort:
    index += 1
    main.write_file(i, ("CY001BC"+str(index))+"_sorted")
