from dataclasses import dataclass
from .DNAString import DNAString


@dataclass
class DNAStringSet:

    sequences: dict[str, DNAString] = None

    def add_sequence(self, name: str, seq: DNAString):
        """ Adds a DNAString to the set with its name. """
        if self.sequences == None:
            self.sequences = {}
        if name not in self.sequences:
            self.sequences[name] = seq
    
    def get_sequence(self, name: str):
        """ Given a name, retrives the associated DNAString. """
        return self.sequences.get(name, None)
