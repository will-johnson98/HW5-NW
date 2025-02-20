# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/Tursiops_truncatus_BRD2.fa")

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    
    alignments = []
    
    score_mm = nw.align(hs_seq, mm_seq)[0]
    alignments.append(("Mouse", score_mm))
    
    score_gg = nw.align(hs_seq, gg_seq)[0]
    alignments.append(("Chicken", score_gg))
    
    score_br = nw.align(hs_seq, br_seq)[0]
    alignments.append(("Shoebill", score_br))
    
    score_tt = nw.align(hs_seq, tt_seq)[0]
    alignments.append(("Dolphin", score_tt))
    
    alignments.sort(key=lambda x: x[1], reverse=True)

    # Print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print("\nSpecies ordered by similarity to human BRD2:")
    for species, score in alignments:
        print(f"{species}")
    print("\nAlignment scores with human BRD2:")
    for species, score in alignments:
        print(f"{species}: {score:.2f}")
    

if __name__ == "__main__":
    main()
