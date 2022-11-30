import main
import matplotlib.pyplot as plt
i = 0
while (i<= 11):
    i += 1
    BC = main.count_nucleotides("BC"+str(i)+"_aligned")
    BC_depth = main.count_depth("BC"+str(i)+"_aligned")
    BC_mutations = main.calculate_mutations(BC, BC_depth)
    BC_dist = main.plot_mutation_distribution(BC_mutations)
    plt.savefig("BC"+str(i)+"_dist.png")
    plt.close()
    BC_substitutions = main.plot_substitution_base_proberbility(BC)
    plt.savefig("BC"+str(i)+"_substitutions.png")
    plt.close()
    BC_mutations = main.plot_mutation_base_proberbility(BC)
    plt.savefig("BC"+str(i)+"_mutations.png")
    plt.close()
# Add function: plot mutations per read
# to CSV for each barcode: mutations per read 01,total amount of bases counted, total mutations,total mutation percentage, sequencing depth, mutation percentage per base