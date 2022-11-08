import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def sort_barcodes(fastq_file, barcode_file):
    "Takes a fasta file with barcodes and a fastq file with data. Returns a dataframe with barcode indexed columns"
    df = pd.DataFrame(pd.read_table(fastq_file, header=None).values.reshape(-1, 4),columns=['read_id', 'seq', '+', 'qual'])
    series_seq = pd.Series(df["seq"])
    with open(barcode_file) as BC:
        r = BC.read()
    df_sorted_by_BC = series_seq.str.extractall(r) # use this list as search groups df_BC["seq"].values
    cut_df = df_sorted_by_BC.str.replace(r'^GCAAAT', "", regex=True)  # Removes the barcode if matching to the regular expression, how can i use this on df?
    return print(cut_df)

def filter(fastq_file):
    "Docstring"
    df = pd.DataFrame(pd.read_table(fastq_file, header=None).values.reshape(-1, 4), columns = ['read_id', 'seq','+', 'qual'])
    df_qual = df["qual"].str.split(pat="",expand=True)
    df_seq = df["seq"].str.split(pat="",expand=True)
    phred_30_df = df_qual.replace(["!",'"',"#","$","%","&","'","(",")","*","+",",","-",".","/","0","1","2","3","4","5","6","7","8","9",":",";","<","=",">"],value = "N")
    masked_sequence_df = df_seq.mask(phred_30_df == "N","N")
    return (masked_sequence_df)

def fastq_pair(fastq_file1,fastq_file2):
    "Docstring"


def merge_txt_files(file1, file2, file3, name):
    "..."
    filenames = file1, file2, file3
    with open(name, "w") as outfile:
        for filename in filenames:
            with open(filename) as infile:
                contents = infile.read()
                outfile.write(contents)
    return name

def count_depth(file):
    "text"
    input_file = open(file)
    depth = input_file.readlines()
    return len(depth)

def count_nucleotides(f_reads):
    "Takes a txt file from a bowtie output with alignment position followed by the sequence seperated by a blank space, returns a dict of dicts with nucleotide position as outer key, nucleotide as inner key and number of nucleotides counted as value"
    dict_of_nucleotide_positions = {}
    input_file = open(f_reads)
    lines = input_file.readlines()
    for line in lines:
        elements = re.split("\s", line)  #splits alignment position from mutations
        alignment_position = elements[0]
        sequence = elements[1]
        for nucleotide in enumerate(sequence): #enumerate to keep track of the alignment position
            position = int(alignment_position)+nucleotide[0]
            if position not in dict_of_nucleotide_positions:
                dict_of_mutations = {"A": 0, "T": 0, "C": 0, "G": 0, "N":0}
                dict_of_nucleotide_positions[position] = dict_of_mutations  # Adds a dict with the nucleotide position as key
            outerkey = position
            innerkey = nucleotide[1]
            if innerkey == "N":
                dict_of_nucleotide_positions[outerkey][innerkey] += 0
            else:
                dict_of_nucleotide_positions[outerkey][innerkey] += 1
    return dict_of_nucleotide_positions

def calculate_mutations(dict, depth):
    "docstring"
    dict_of_mutation_percentages = {}
    items = dict.items()
    for item in items:
        position = item[0]
        dict_of_nucleotides = item[1]
        list_of_nucleotides = dict_of_nucleotides.values()
        wt = max(list_of_nucleotides)
        mutations = sum(list_of_nucleotides)-wt
        mutation_percentage = [mutations/sum(list_of_nucleotides)]
        if wt > 0.1*depth and position not in dict_of_mutation_percentages:
            dict_of_mutation_percentages[position] = mutation_percentage
    return dict_of_mutation_percentages

def total_mutations(dict):
    "docstring"
    items = dict.items()
    list_of_mutations = []
    list_of_wt = []
    for item in items:
        dict_of_nucleotides = item[1]
        list_of_nucleotides = dict_of_nucleotides.values()
        wt = max(list_of_nucleotides)
        mutation = sum(list_of_nucleotides)-wt
        list_of_wt.append(wt)
        list_of_mutations.append(mutation)
    return sum(list_of_wt), sum(list_of_mutations), sum(list_of_mutations)/sum(list_of_wt)

def plot_mutation_distribution(dict_of_mutation_percentages):
    "Docstring"
    lists = sorted(dict_of_mutation_percentages.items())
    x, y = zip(*lists)
    plt.scatter(x, y, s=4, color="black")
    return

def plot_substitution_base_proberbility(dict):
    "Docstring"
    dict_of_wt = {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0}
    A = { "T": 0, "C": 0, "G": 0}
    T = {"A": 0, "C": 0, "G": 0}
    G = {"T": 0, "C": 0, "A": 0}
    C = {"T": 0, "A": 0, "G": 0}
    for key, value in dict.items():
        wt = max(value, key=value.get)
        count = max(value.values())
        dict_of_wt[wt] += count
        if "A" == wt:
            substitution = list(value.values())
            A["T"] += substitution[1]
            A["C"] += substitution[2]
            A["G"] += substitution[3]
        if "T" == wt:
            substitution = list(value.values())
            T["A"] += substitution[0]
            T["C"] += substitution[2]
            T["G"] += substitution[3]
        if "C" == wt:
            substitution = list(value.values())
            C["A"] += substitution[0]
            C["G"] += substitution[3]
            C["T"] += substitution[1]
        if "G" == wt:
            substitution = list(value.values())
            G["A"] += substitution[0]
            G["T"] += substitution[1]
            G["C"] += substitution[2]
    AX = ["T","C","G"]
    TX = ["A","C","G"]
    CX = ["A","G","T"]
    GX = ["A","T","C"]
    ax1 = plt.subplot2grid((2, 2), (0, 0))
    ax2 = plt.subplot2grid((2, 2), (0, 1))
    ax3 = plt.subplot2grid((2, 2), (1, 0))
    ax4 = plt.subplot2grid((2, 2), (1, 1))
    A = {k: v / total for total in (sum(A.values()),) for k, v in A.items()}
    T = {k: v / total for total in (sum(T.values()),) for k, v in T.items()}
    C = {k: v / total for total in (sum(C.values()),) for k, v in C.items()}
    G = {k: v / total for total in (sum(G.values()),) for k, v in G.items()}
    ax1.bar(AX, A.values(), color = ["blue"])
    ax1.set_title('A')
    ax2.bar(TX, T.values(),color="green")
    ax2.set_title('T')
    ax3.bar(CX, C.values(),color="orange")
    ax3.set_title('C')
    ax4.bar(GX, G.values(),color="red")
    ax4.set_title('G')
    plt.tight_layout()
    return

def plot_mutation_base_proberbility(dict):
        "Docstring"
        dict_of_wt = {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0}
        A = {"T": 0, "C": 0, "G": 0}
        T = {"A": 0, "C": 0, "G": 0}
        G = {"T": 0, "C": 0, "A": 0}
        C = {"T": 0, "A": 0, "G": 0}
        for key, value in dict.items():
            wt = max(value, key=value.get)
            count = max(value.values())
            dict_of_wt[wt] += count
            if "A" == wt:
                substitution = list(value.values())
                A["T"] += substitution[1]
                A["C"] += substitution[2]
                A["G"] += substitution[3]
            if "T" == wt:
                substitution = list(value.values())
                T["A"] += substitution[0]
                T["C"] += substitution[2]
                T["G"] += substitution[3]
            if "C" == wt:
                substitution = list(value.values())
                C["A"] += substitution[0]
                C["G"] += substitution[3]
                C["T"] += substitution[1]
            if "G" == wt:
                substitution = list(value.values())
                G["A"] += substitution[0]
                G["T"] += substitution[1]
                G["C"] += substitution[2]
            if "N" == wt:
                None
        x1 = ["A", "T", "C", "G"]
        y1 = [sum(A.values()) / dict_of_wt["A"], sum(T.values()) / dict_of_wt["T"], sum(C.values()) / dict_of_wt["C"],sum(G.values()) / dict_of_wt["G"]]  # height of bar is number of substitutions divided by number of wt bases counted
        plt.bar(x1, y1, color=["blue", "green", "orange", "red"])
        plt.tight_layout()
        return
