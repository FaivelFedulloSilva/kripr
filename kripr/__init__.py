from .bam_handler import BAMhandler
from .BioString import BioString
from .DNAString import DNAString
from .DNAStringSet import DNAStringSet
from .fasta_handler import FastaHandler
from .gtf_handler import GTFhandler
from .gtf_object import GTFobject
from .RNAString import RNAString
from .main import   get_feature_sequence, \
                    get_feature_sequence_parallel, \
                    save_to_file, \
                    get_mapped_reads_from_bam,\
                    get_deseq2_matrix, \
                    normalize_base_with_deseq2, \
                    get_cds_sequence