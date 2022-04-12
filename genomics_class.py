import  numpy as np
import pandas as pd


class Read_Genome():
    def __init__(self, FASTA):
        self.DNA_nucleotides = list('AGCT')
        self.complement = dict(zip(self.DNA_nucleotides, self.DNA_nucleotides[::-1]))
        self.codons = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V', 'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V', 'UUA': 'L', 
        'CUA': 'L', 'AUA': 'I', 'GUA': 'V', 'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V', 'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 
        'GCU': 'A', 'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A', 'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 'UCG': 'S', 
        'CCG': 'P', 'ACG': 'T', 'GCG': 'A', 'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D', 'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 
        'GAC': 'D', 'UAA': 'Stop', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E', 'UAG': 'Stop', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E', 
        'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G', 'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G', 'UGA': 'Stop', 'CGA': 'R', 
        'AGA': 'R', 'GGA': 'G', 'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}

        # Given a fasta file, returns a dictionary of labels as keys and sequences as values
        f = open(FASTA, "r")
        txt = f.read()
        fasta = txt.split(">")
        dict_FASTA = {}
        for i in fasta:
            seq =  i.split("\n")
            if len(seq) > 1:
                dict_FASTA[seq[0 ]] = ''.join(seq[1:])

        self.dict_FASTA = dict_FASTA

        dict_FASTA = self.dict_FASTA.values()
        # self.arr_FASTA = np.array([list(x) for x in dict_FASTA], dtype=object)
 
    def _gc_content(seq):
        """Finds % of CG content in a single sequence"""
        try:
            return (seq.count("G") + seq.count("C"))/len(seq)*100
        except(ZeroDivisionError):
            return "Division by zero occured when dividing by the total number of nucleotides.\nCheck that your file is in the correct format."

    def gc_content(self):
        """Finds the sequence containing the greatest GC content.
        Returns the name of sequence and the % of GC as a tuple"""

        n_GC = {}
        for i in self.dict_FASTA:
            gcl = Read_Genome._gc_content( self.dict_FASTA[i] )

            if type(gcl) == str: # check if an exception occured and return the error message
                return gcl
            else:                 # else add element to the list
                n_GC[i] = gcl

        name = max(n_GC, key=n_GC.get)
        percentage = max(n_GC.values())

        return '{}\n\n{:.2f}% of GC'.format(name, percentage)

    def count_nucleotides(self):
        """Takes a FASTA file and returns the count of each Nucleotide on every row
        in a list form"""
        
        arr_FASTA = np.array([list(x) for x in self.dict_FASTA.values()], dtype=object)
        cons_li = []
        for nuc in self.DNA_nucleotides:
            counts = np.count_nonzero(arr_FASTA==nuc, axis=0)
            cons_li.append(counts)

        return cons_li

    def common_ancestor(self):
        """Takes a FASTA file and returns count of each Nucleotide on every row in a dictionary form,
        and a possible common ancestor gene"""
        cons_li = Read_Genome.count_nucleotides(self)
        ancestor = []
        for row in np.array(cons_li, dtype=object).T:
            ancestor.append(np.argmax(row))
        map_idx = dict (zip( [0,1,2,3], self.DNA_nucleotides))
        ancestor = list(map( map_idx.get, ancestor))
        ancestor_str = ''.join(ancestor)

        return ancestor_str

    def _frequency(self, seq):
        """Find the frequency of each Nucleotide in a sequence"""
        freq = []
        for nuc in self.DNA_nucleotides:
            freq.append( seq.count(nuc) )
        return freq
    
    def frequency(self):
        freq = {}
        for key in self.dict_FASTA.keys():
            freq[key] = self._frequency(self.dict_FASTA[key])
        table = pd.DataFrame.from_dict(freq, orient='index', columns=self.DNA_nucleotides)
        return table
