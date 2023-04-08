from gtfparse import read_gtf 
import polars as pl
from .gtf_object import GTFobject

pl.Config.set_tbl_cols(16)
class GTFhandler:

    def __init__(self, path: str) -> None:
        self.__path: str = path
        self.__gtf_object = read_gtf(path)
    
    def get_gtf_object(self):
        return GTFobject(self.__gtf_object)

    def reload(self):
        self.__gtf_object = read_gtf(self.__path)


# gtf = GTFhandler('./Data/genesFiltrada.gtf')
# # print(gtf._gtf_object.describe())
# cds_df = gtf.filter_by_feature('CDS')     
# print(cds_df.sample(10))

if __name__ == '__main__':
    GTF_PATH = 'Data/genesFiltrada.gtf'
    gtf_file_handler = GTFhandler(GTF_PATH)
    gtf = gtf_file_handler.get_gtf_object()
    gtf_exon: GTFobject = gtf.filter_by_feature('exon')
    gtf_cds: GTFobject = gtf.filter_by_feature('CDS')
    print(gtf_exon.get_transcripts_ids())