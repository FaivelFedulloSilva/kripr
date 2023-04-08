from .DNAString import DNAString
from .DNAStringSet import DNAStringSet
import polars as pl

class FastaHandler:
    def __init__(self, path = '', polars=False) -> None:
        def read_fasta(path: str):
            referenceDNA = DNAStringSet()
            current_reference = []
            current_id = None
            with open(path, 'r') as file:
                line = file.readline()
                while line != '':
                    if line[0] == '>':
                        if current_id != None:
                            referenceDNA.add_sequence(current_id, DNAString(''.join(current_reference)))
                        current_id = line[1:].split(' ')[0].strip()
                        current_reference= []
                    else:
                        current_reference.append(line.strip())
                    line = file.readline()
            return referenceDNA
        if polars:
            self.reference = pl.scan_csv(file='reference.pls', sep='\t', has_header=False)
        else:
            self.reference = read_fasta(path)

    def get_reference(self):
        return self.reference 

    def save_as_polars_file(self):
        with open('reference.pls', 'x') as file:
            for key,value in self.reference.sequences.items():
                file.write(f'{key}\t{value.get_as_string()}\n')


      
if __name__ == '__main__':
    faHandler = FastaHandler(r'./Data/reference/hg38.fa')
    faHandler.save_as_polars_file()