import re
import matplotlib.pyplot as plt
import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


def mask(fastq_file):
    """Takes a fastq file and returns a pandas dataframe with the values, basecalls with a phredscore under 30 are masked with 'N' """
    df = pd.DataFrame(pd.read_table(fastq_file, header=None, comment=" ",).values.reshape(-1, 4), columns=['read_id', 'seq', '+', 'qual']).drop(columns="+")  # reads a fastq file into a dataframe where each entry is indexed at the column as read_id,seq and qual
    df_qual = df["qual"].str.split(pat="", expand=True)  # splitting the quality into each individual charecter to assess their value
    df_seq = df["seq"].str.split(pat="", expand=True).fillna("") # splitting the sequence into each individual charecters, filling the gaps of shorter sequences with empty strings
    phred_30_df = df_qual.replace(["!", '"', "#", "$", "%", "&", "'", "(", ")", "*", "+", ",", "-", ".", "/", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ":", ";", "<", "=", ">"], value = "N") #
    masked_sequence_df = df_seq.iloc[:, 7:].mask(phred_30_df == "N", "N")  # maks base calls with a phred score of 30 or lower, disregards the first 6 basecalls (the barcode)
    frames = [df_seq.iloc[:, :7], masked_sequence_df]  # joins the barcode and the masked sequence
    x = pd.concat(frames, axis=1)  # joins the barcode and the masked sequence
    df["seq"] = x.apply((lambda x: "".join(x)), axis=1)  # why are random N's added to the shorter sequences?
    return (df)


def sort_barcodes(df1,df2):
    """Takes a fasta file with barcodes and a fastq file with data. Returns a dataframe with barcode indexed columns"""
    series_seq = pd.Series(df1["seq"])
    series_seq2 = pd.Series(df2["seq"])
    df_sorted_by_BC1 = series_seq.str.extractall("(?P<BC1>^GCAAAA.+)").reset_index(drop=True)
    df_sorted_by_BC2 = series_seq.str.extractall("(?P<BC2>^GCAAAG.+)").reset_index(drop=True)
    df_sorted_by_BC3 = series_seq.str.extractall("(?P<BC3>^GCAAAC.+)").reset_index(drop=True)
    df_sorted_by_BC4 = series_seq.str.extractall("(?P<BC4>^GCAAAT.+)").reset_index(drop=True)
    df_sorted_by_BC5 = series_seq.str.extractall("(?P<BC5>^GCAAGA.+)").reset_index(drop=True)
    df_sorted_by_BC6 = series_seq.str.extractall("(?P<BC6>^GCAACA.+)").reset_index(drop=True)
    df_sorted_by_BC7 = series_seq.str.extractall("(?P<BC7>^GCAATA.+)").reset_index(drop=True)
    df_sorted_by_BC8 = series_seq.str.extractall("(?P<BC8>^GCAGAA.+)").reset_index(drop=True)
    df_sorted_by_BC9 = series_seq.str.extractall("(?P<BC9>^GCACAAA.+)").reset_index(drop=True)
    df_sorted_by_BC10 = series_seq.str.extractall("(?P<BC10>^GCATAA.+)").reset_index(drop=True)
    df_sorted_by_BC1_2 = series_seq2.str.extractall("(?P<BC1>^GCAAAA.+)").reset_index(drop=True)
    df_sorted_by_BC2_2 = series_seq2.str.extractall("(?P<BC2>^GCAAAG.+)").reset_index(drop=True)
    df_sorted_by_BC3_2 = series_seq2.str.extractall("(?P<BC3>^GCAAAC.+)").reset_index(drop=True)
    df_sorted_by_BC4_2 = series_seq2.str.extractall("(?P<BC4>^GCAAAT.+)").reset_index(drop=True)
    df_sorted_by_BC5_2 = series_seq2.str.extractall("(?P<BC5>^GCAAGA.+)").reset_index(drop=True)
    df_sorted_by_BC6_2 = series_seq2.str.extractall("(?P<BC6>^GCAACA.+)").reset_index(drop=True)
    df_sorted_by_BC7_2 = series_seq2.str.extractall("(?P<BC7>^GCAATA.+)").reset_index(drop=True)
    df_sorted_by_BC8_2 = series_seq2.str.extractall("(?P<BC8>^GCAGAA.+)").reset_index(drop=True)
    df_sorted_by_BC9_2 = series_seq2.str.extractall("(?P<BC9>^GCACAAA.+)").reset_index(drop=True)
    df_sorted_by_BC10_2 = series_seq2.str.extractall("(?P<BC10>^GCATAA.+)").reset_index(drop=True)
    list_of_dfs = [df_sorted_by_BC1,df_sorted_by_BC1_2,df_sorted_by_BC2,df_sorted_by_BC2_2,df_sorted_by_BC3,df_sorted_by_BC3_2, df_sorted_by_BC4,df_sorted_by_BC4_2,df_sorted_by_BC5,df_sorted_by_BC5_2,df_sorted_by_BC6,df_sorted_by_BC6_2,df_sorted_by_BC7,df_sorted_by_BC7_2,df_sorted_by_BC8,df_sorted_by_BC8_2,df_sorted_by_BC9,df_sorted_by_BC9_2,df_sorted_by_BC10,df_sorted_by_BC10_2]
    BC_df = pd.concat(list_of_dfs, axis=1).replace(r'^^[A-Z]{6}', "", regex=True) # Removes the barcode ie. the first six letters of the sequence
    combined_BC1 = pd.concat([df_sorted_by_BC1,df_sorted_by_BC1_2],axis=0).replace(r'^^[A-Z]{6}', "", regex=True)
    BC1 = combined_BC1["BC1"].values.tolist()# change to output all barcodes like this
    return (BC1)



def write_file(list_of_sequences, filename):
    """Takes a list of sequences and a filename as a string, returns a fastafile"""
    with open(filename, "w") as outfile:
        for i in list_of_sequences:
            outfile.write(i+"\n")
    return print(list_of_sequences)


def count_depth(file):
    "text"
    input_file = open(file)
    depth = input_file.readlines()
    return len(depth)


def count_nucleotides(f_reads):
    """Takes a txt file from a bowtie output with alignment position followed by the sequence seperated by a blank space, returns a dict of dicts with nucleotide position as outer key, nucleotide as inner key and number of nucleotides counted as value"""
    dict_of_nucleotide_positions = {}
    input_file = open(f_reads)
    lines = input_file.readlines()
    for line in lines:
        elements = re.split("\s", line)  #splits alignment position from mutations
        alignment_position = elements[3]
        sequence = elements[9]
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
    """docstring"""
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
    """docstring"""
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
    """Docstring"""
    lists = sorted(dict_of_mutation_percentages.items())
    x, y = zip(*lists)
    plt.scatter(x, y, s=4, color="black")
    return


def plot_substitution_base_proberbility(dict):
    """Docstring"""
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
        """Docstring"""
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
