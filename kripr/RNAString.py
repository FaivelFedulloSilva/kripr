from .BioString import BioString 


complementary_bases = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}

class RNAString(BioString):
    
    def complement(self):
        complementary_strand = ''
        for x in self.seq:
            complementary_strand += complementary_bases[x]
        return complementary_strand

    def reverse(self):
        return self.seq[::-1]

    def reverse_complement(self):
        return self.complement()[::-1]

    def count_base(self, base: chr) -> int:
        #TODO Cambiar lista por extension a definicion IUPAC para DNA
        if base not in ['A', 'U', 'G', 'C']:
            raise ValueError("Only DNA IUPAC codes are accepted")
        return self.seq.count(base)

    def GC_content(self):
        gc_count = 0
        for b in self.seq:
            if b in ['G', 'C']:
                gc_count += 1
        return gc_count/len(self.seq)

    def get_sub_sequence(self, start: int, end: int):
        """ Returns a new RNAString object containing the subsequence from start to end (both inclusive). """
        return RNAString(self.seq[start-1:end])

    #TODO Agregar RNA -> DNA
    #TODO gregar RNA -> AA