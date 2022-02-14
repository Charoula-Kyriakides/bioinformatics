import  numpy as np
import itertools
import re

DNA_nucleotides = list('AGCT')
complement = dict(zip(DNA_nucleotides, DNA_nucleotides[::-1]))

codons = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V', 'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V', 'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V', 'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V', 'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A', 'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A', 'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A', 'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D', 'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D', 'UAA': 'Stop', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E', 'UAG': 'Stop', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E', 'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G', 'UGC': 'C', 'CGC': 'R', 
'AGC': 'S', 'GGC': 'G', 'UGA': 'Stop', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G', 'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}

def read_FASTA(FASTA):
    """Given a fasta file, returns a dictionary of labels as keys and sequences as values"""
    f = open(FASTA, "r")
    txt = f.read()
    fasta = txt.split(">")
    dict_fasta = {}
    for i in fasta:
        seq =  i.split("\n")
        if len(seq) > 1:
            dict_fasta[seq[0 ]] = ''.join(seq[1:])
    return dict_fasta

def array_FASTA(FASTA):
    """Opens a FASTA file and returns the genes contained in an array form"""
    dict_FASTA = read_FASTA(FASTA).values()
    arr_FASTA = np.array([list(x) for x in dict_FASTA], dtype=object)
    return arr_FASTA

def frequency(seq):
    """Find the frequency of each Nucleotide in a sequence"""
    freq = {}
    for nuc in DNA_nucleotides:
        freq[nuc] = seq.count(nuc)
    return freq

def transcription(seq, file=False, label=None):
    """Trascribe DNA to RNA"""
    seq = seq.replace("T","U")
    if file:
        if label:
            f = open("RNA_seq_{}.txt".format(label),"w")
        else:
            label = ""
            f = open("RNA_seq_{}.txt".format(label),"w")
        f.write(seq)
        f.close()
        return f
    return seq

def reverse_complement(seq):
    """Produce the reverse complement of a DNA sequence"""
    reverse = ""
    for nuc in seq:
        reverse += complement[nuc]
    return reverse[::-1]

def reverse_complement_multiple(FASTA):
    """Apply the reverse_complement function to the values of a dictionary"""
    dict_FASTA = read_FASTA(FASTA)
    complements = {k:reverse_complement(v) for k, v in dict_FASTA.items()}
    return complements

def GC_content(seq):
    """Finds % of CG content in a single sequence"""
    return (seq.count("G") + seq.count("C"))/len(seq)*100

def GC_content_dict(FASTA):
    """Finds the sequence containing the greatest GC content.
    Returns the name of sequence and the % of GC as a tuple"""
    dict_fasta = read_FASTA(FASTA)
    n_GC = {}
    for i in dict_fasta:
        n_GC[i] = GC_content( dict_fasta[i] )
    return max(n_GC, key=n_GC.get), max(n_GC.values())
        
def Hemming_distance(seq):
    """Number of bases in which two strings differ"""
    seq = seq.split("\n")
    ar1 = np.array(list(seq[0]))
    ar2 = np.array(list(seq[1]))
    return len(ar1) - sum( ar1 == ar2) 

def Mendel_first_law(DD=0, DR =0, RR=0):
    """Find the probability that an allel which includes D 
    will be produced from a random selection of two allels
    """
    tot = DD + RR + DR

    #Probabilities that allel is RR
    p_RR = (RR / tot)*(RR-1)/(tot-1)
    p_mixed = (RR/tot)*(.5*DR)/(tot-1) + (.5*DR/tot)*RR/(tot-1)
    p_DR = (.5*DR/tot)*(.5*(DR-1)/(tot-1))

    # So probability that allel is DR or DD
    prop = 1 - ( p_RR + p_mixed + p_DR )
                    
    return (prop)

def permutation(n):
    """Finds all possible permutations of a list of length n.
    Input: integer n, 
    Output: integer, numpy array
    """
    ordered = np.arange(1,n+1)
    permutations = np.array(list( itertools.permutations(ordered) ))
    return len(permutations), permutations

def motif(file_name):
    """Find the positions of a sub sequence (motif) in a gene"""
    motif1 = open(file_name)
    seq = motif1.readline()
    sub_seq = motif1.readline().strip('\n')
    f_motif = open('motif', 'w')
    start = 0
    while True:
        start = seq.find(sub_seq, start, -1) + 1
        if start != 0:
            f_motif.write(str(start) + ' ')
        else:
            break
    f_motif.close()
    return f_motif

def count_nucleotides(FASTA):
    """Takes a FASTA file and returns the count of each Nucleotide on every row
    in a dictionary and list form"""
    genes = array_FASTA(FASTA)
    cons_li = []
    f = open('ancestor.txt', 'w')

    for nuc in DNA_nucleotides:
        counts = np.count_nonzero(genes==nuc, axis=0)
        f.write(nuc + ': ' + str(counts).strip('[').strip(']') + '\n')
        cons_li.append(counts)
    return cons_li, f

def common_ancestor(FASTA):
    """Takes a FASTA file and returns count of each Nucleotide on every row in a dictionary form,
    and a possible common ancestor gene"""
    cons_li, f = count_nucleotides(FASTA)
    ancestor = []
    for row in np.array(cons_li, dtype=object).T:
        ancestor.append(np.argmax(row))
    map_idx = dict (zip( [0,1,2,3], DNA_nucleotides))
    ancestor = list(map( map_idx.get, ancestor))
    f.write(''.join(ancestor))
    f.close()
    return f 

def overlap_graphs(FASTA,n=3):
    """Return a file of pairsof DNA strings which overlap by n number of bases;
        their prefix and suffix must match, but the whole string must not. """
    orig_seq = read_FASTA(FASTA)
    f = open('overlap.txt','w')
    for key1 in orig_seq:
        for key2 in orig_seq:
            if orig_seq[key1][-n:] == orig_seq[key2][:n] and key1 != key2:
                f.write(key1 + ' ' + key2 + '\n')
    f.close()
    return f

def decode_protein(RNA_string):
    """Takes and RNA in string format and returns the decoded protein"""
    triples = re.findall('...?', RNA_string)
    protein = ''.join(list( map(codons.get, triples)))
    return protein

def splice(FASTA):
    """Splices a DNA and returns the protein.
    Input: fasta file with the RNA first and introns later.
    Returns a string.
    """
    DNA_coll = list(read_FASTA(FASTA).values())
    string = DNA_coll[0]
    introns = DNA_coll[1:]
    for intron in introns:
        string = string.replace(intron,'')
    RNA = transcription(string, file=False)
    prot = decode_protein(RNA)
    return prot.replace('Stop','')

def spliced_motif(FASTA):
    """Takes fasta of two sequences: 1 main and 1 subsequence.
    returns the indexes of a spliced motif in a list form.
    """
    seqs = array_FASTA(FASTA)
    idx = -1
    indexes = []
    for base in seqs[1]:
        idx = seqs[0].index(base, idx+1)
        indexes.append(idx+1)
    return indexes