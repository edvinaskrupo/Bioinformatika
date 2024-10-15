from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

# kodonu lentele
codons = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "_", "TAG": "_",
    "TGC": "C", "TGT": "C", "TGA": "_", "TGG": "W"
}

# nustatome start ir stop kodonus
START_CODON = "ATG"
STOP_CODONS = ["TAA", "TAG", "TGA"]

# ieskome ilgiausiu kodonu poru
def identify_longest_orf_pairs(seq):
    valid_pairs = []

    def find_longest_pairs(dna_seq):
        start_positions = [i for i in range(len(dna_seq)) if dna_seq[i:i+3] == START_CODON]
        for start in start_positions:
            farthest_stop = -1
            for stop_codon in STOP_CODONS:
                stop_index = dna_seq.find(stop_codon, start + 3)
                if stop_index != -1 and all(dna_seq[start + 3:stop_index].find(stop) == -1 for stop in STOP_CODONS):
                    if stop_index > farthest_stop:
                        farthest_stop = stop_index
            if farthest_stop != -1:
                valid_pairs.append((start, farthest_stop + 3))
                
    # analizuojame tiesioginę seką ir reverse komplementą
    find_longest_pairs(seq)
    find_longest_pairs(str(Seq(seq).reverse_complement()))
    
    # atfiltruojame trumpesnius nei 100bp kodonus
    return [pair for pair in valid_pairs if pair[1] - pair[0] >= 100]

# konvertuojame i baltymo seka
def translate_to_protein(seq, pairs):
    proteins = []
    for start, end in pairs:
        proteins.append("".join([codons[seq[i:i+3]] for i in range(start, end, 3)]))
    return proteins
