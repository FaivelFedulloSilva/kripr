import polars as pl
from polars import DataFrame

class GTFobject:

    def __init__(self, gtf: DataFrame) -> None:
        self._gtf = gtf


    def filter_by_feature(self, feature: str):
        """Return a new GTFobject that contains only the rows of type feature"""

        filtered_gtf = self._gtf.filter(
            pl.col("feature") == feature
        ) 
        return GTFobject(filtered_gtf)

    
    def get_transcripts_ids(self) -> list[str]:
        """Return a list of all the transcripts_id"""
        print(type(self._gtf))
        query = (self._gtf
                .lazy()
                .select(pl.col('transcript_id')).unique()
            )
        
        return query.collect().to_dict(as_series=False)['transcript_id']


    def get_cds_data(self, transcript_id):
        df = self._gtf.filter(self._gtf['transcript_id'] == transcript_id)
        return (
            df.lazy()
            .groupby(['protein_id'])
            .agg(
                [
                    (pl.col('seqname').unique().take(0)),
                    (pl.col('start')),
                    (pl.col('end')),
                    (pl.col('strand').unique().take(0)),
                    (pl.col('transcript_id').unique().take(0)),
                    (pl.col('gene_id').unique().take(0)),
                ]
            )
        ).collect()

    def get_transcripts_data(self, transcripts_ids: list[str]) -> DataFrame:
        query = (self._gtf
                .lazy()
                .groupby(by='transcript_id')
                .agg(
                    [pl.col('seqname').unique().first(), pl.col('strand'), pl.col('start'),  pl.col('end'), pl.col('gene_id').unique().first()]
                    )
                .filter(pl.col('transcript_id').is_in(transcripts_ids))
            ).collect()

        return query
    
    def get_transcript_id_data(self, transcript_id: str) -> DataFrame:
        return self.get_transcripts_data([transcript_id])

    def get_full_transcripts_data(self, transcripts_ids: list[str]) -> DataFrame:
        query = (self._gtf
                .lazy()
                .filter(pl.col('transcript_id').is_in(transcripts_ids))
                .groupby(by='p_id')
                .agg(
                    [
                        pl.col('seqname').unique().first(),
                        pl.col('feature'),
                        pl.col('strand'),
                        pl.col('start'),  
                        pl.col('end')
                    ]
                )
            ).collect()

        return query

    def get_full_transcript_data(self, transcript_id: str) -> DataFrame:
        return self.get_full_transcripts_data([transcript_id])
    
    def generate_igv_file(self, gene_id_list: list[str]):
        query = (self._gtf
                .lazy()
                .filter(pl.col('gene_id').is_in(gene_id_list))
                .groupby(by='gene_id')
                .agg(
                    [
                        pl.col('seqname').unique().first(),
                        pl.col('transcript_id'),
                        pl.col('feature'),
                        pl.col('strand'),
                        pl.col('start'),  
                        pl.col('end')
                    ]
                )
            ).collect()

        return query

    def format_gene_igv(self, gene_id: str):
        dt = self.generate_igv_file([gene_id])
        iter_df = dt.iterrows(named=True)
        genes = {}
        for row in iter_df:
            transcripts = {}
            for i in range(len(row.transcript_id)):
                if row.transcript_id[i] not in transcripts:
                    transcripts[row.transcript_id[i]] = {'exon': [], 'CDS': [], 'start_codon': [], 'stop_codon': []}
                if row.strand[i] == '+':
                    transcripts[row.transcript_id[i]][row.feature[i]].append([row.start[i],row.end[i]])
                else:
                    transcripts[row.transcript_id[i]][row.feature[i]].append([row.end[i],row.start[i]])
            genes[row.gene_id] = transcripts
        return genes

    
