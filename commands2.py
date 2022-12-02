import main
import matplotlib.pyplot as plt
i = 0
while (i<= 11):
    i += 1
    BC = main.count_nucleotides("CY001BC"+str(i)+"_aligned")
    BC_depth = main.count_depth("CY001BC"+str(i)+"_aligned")
    BC_mutations = main.calculate_mutations(BC, BC_depth)
    BC_dist = main.plot_mutation_distribution(BC_mutations)
    plt.savefig("CY001BC"+str(i)+"_dist.png")
    plt.close()
    BC_substitutions = main.plot_substitution_base_proberbility(BC)
    plt.savefig("CY001BC"+str(i)+"_substitutions.png")
    plt.close()
    BC_dist = main.plot_mutation_base_proberbility(BC)
    plt.savefig("CY001BC"+str(i)+"_mutations.png")
    plt.close()
    total = main.total_mutations(BC)
    mutations = main.mutation_per_read("CY001BC"+str(i)+"_aligned")
    csv = main.CSV(BC_mutations,"CY001file"+str(i)+".txt",BC_depth,total,mutations)
