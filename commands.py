import main
mask1 = main.mask("CY002_1_1M.fq")
mask2 = main.mask("CY002_2_1M.fq")
sort = main.sort_barcodes(mask1,mask2)
index = 0
for i in sort:
    index += 1
    main.write_file(i, ("BC"+str(index))+"_sorted")