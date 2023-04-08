from .BioString import BioString

complementary_bases = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

class DNAString(BioString):

    def get_as_string(self):
        return self.seq

    def complement(self):
        complementary_strand = ''
        for x in self.seq:
            complementary_strand += complementary_bases[x.upper()]
        return DNAString(complementary_strand)

    def reverse(self):
        return DNAString(self.seq[::-1])

    def reverse_complement(self):
        return DNAString(self.complement().get_as_string()[::-1])

    def count_base(self, base: chr) -> int:
        #TODO Cambiar lista por extension a definicion IUPAC para DNA
        if base not in ['A', 'T', 'G', 'C']:
            raise ValueError("Only DNA IUPAC codes are accepted")
        return self.seq.count(base)

    def GC_content(self):
        gc_count = 0
        for b in self.seq:
            if b in ['G', 'C']:
                gc_count += 1
        return gc_count/len(self.seq)

    def get_sub_sequence(self, start: int, end: int):
        """ Returns a new DNAString object containing the subsequence from start to end (both inclusive). """
        return DNAString(self.seq[start-1:end])
    
    # TODO Agregar DNA -> RNA

# dna = DNAString('AAAATGGC')
# print(dna.complement())
# print(dna.reverse())
# print(dna.reverse_complement())
# print(dna.get_sub_sequence(0,4))
# print(dna.count_base('A'))