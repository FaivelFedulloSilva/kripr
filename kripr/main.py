
from collections import defaultdict
from multiprocessing import Pool
import os
import shutil
from kripr import DNAStringSet, BAMhandler, GTFobject, GTFhandler
import numpy as np
import polars as pl
from polars import DataFrame
import pysam

def get_second_elements(l: list[str,str])->list[int]:
    return [int(e[1]) for e in l]

def general_per_base_coverage(cov):
    """ Return the general per base coverage of the transcript"""
    return np.count_nonzero(cov)/np.size(cov)

def save_to_file(
    dt: DataFrame, 
    file_name: str,
    file_type: str = 'parquet'
):
    if file_type == 'parquet':
        dt.write_parquet(f'{file_name}.parquet')
    elif file_type == 'json':
        dt.write_json(f'{file_name}.json')
    elif file_type == 'csv':
        dt.write_csv(f'{file_name}.csv')

def read_from_file(
    file_path: str,
    file_type: str = 'parquet'
):
    dt = None
    if file_type == 'parquet':
        dt = pl.read_parquet(file_path)
    elif file_type == 'json':
        dt = pl.read_json(file_path)
    elif file_type == 'csv':
        dt = pl.read_csv(file_path)
    return dt

def task(row, ref_seq, bam):
    # print(f'Working from process {getpid()} - Chromosome {bam}')
    bam_handler = BAMhandler(bam)
    feature_sequence = {}
    feature_seq = []
    feature_depth = []
    for index in range(len(row[3])):
        seq = ref_seq.get_sub_sequence(row[3][index], row[4][index])
        if row[2][index] == '-':
            seq = seq.reverse_complement()
        cov = bam_handler.region_coverage(row[1], row[3][index], row[4][index])
        feature_depth += get_second_elements(cov)
        feature_seq.append(seq.get_as_string())
    feature_sequence['transcript_id'] = row[0]
    feature_sequence['gene_id'] = row[5]
    feature_sequence['coverage_percentage'] = general_per_base_coverage(feature_depth)
    feature_sequence['sequence'] = "".join(feature_seq)
    feature_sequence['coverage_per_base'] = (feature_depth)
    return feature_sequence

def get_feature_sequence_parallel(
    chunk_size, 
    transcript_ids: list[str], 
    gtf: GTFobject = None,
    dna: DNAStringSet = None, 
    bam: BAMhandler = None,
    bams_folder: str = None,
) -> DataFrame:        
        transcripts_data = gtf.get_transcripts_data(transcript_ids)
        iter_df = list(transcripts_data.iterrows())
        feature_sequences = {
            'transcript_id': [],
            'gene_id': [],
            'coverage_percentage': [],
            'sequence': [],
            'coverage_per_base': []
            }
        items = []
        for row in iter_df:
            items.append((list(row), dna.get_sequence(row[1]), f'{bams_folder}/{row[1]}.bam'))
        # create the process pool
        with Pool(processes=3) as pool:
            results = pool.starmap(task, items, chunksize=chunk_size)
            
        return pl.DataFrame(results)


#TODO: Que hace esta funcion? Borrarla?
def get_mapped_reads_from_bam(bam_path: str, option: int):
    if option == 0:
        with pysam.AlignmentFile(bam_path, "rb") as bamfile:
            mapped_reads = 0
            unmapped_reads = 0
            for read in bamfile.fetch(until_eof=True):
                if read.is_unmapped:
                    unmapped_reads += 1
                else:
                    mapped_reads += 1
        print(mapped_reads, unmapped_reads)
        return mapped_reads, unmapped_reads
    else:
        response = pysam.flagstat(bam_path)
        print(response)



def get_feature_sequence(
    transcript_ids: list[str], 
    gtf: GTFobject = None,
    dna: DNAStringSet = None, 
    bam: BAMhandler = None
) -> DataFrame:
    '''
    The row is a tuple of the following for:
        (
        'transcript_id': str,   [0]
        'seqname': str,         [1]
        'strand': list[str],    [2] 
        'start': list[int]',    [3]
        'end': list[int],       [4]
        'gene_id': str,         [5]
        )
    '''    
    transcripts_data = gtf.get_transcripts_data(transcript_ids)
    iter_df = transcripts_data.iterrows()
    feature_sequences = {
        'transcript_id': [],
        'gene_id': [],
        'coverage_percentage': [],
        'sequence': [],
        'coverage_per_base': []
        }
    sum_count = 0
    for row in iter_df:
        feature_seq = []
        feature_depth = []
        chromosome = f'chr{row[1]}' if row[1].isdigit() else row[1]
        ref_seq = dna.get_sequence(chromosome)
        for index in range(len(row[3])):
            seq = ref_seq.get_sub_sequence(row[3][index], row[4][index])
            if row[2][index] == '-':
                seq = seq.reverse_complement()
            cov = bam.region_coverage(chromosome, row[3][index], row[4][index])
            feature_depth += get_second_elements(cov)
            feature_seq.append(seq.get_as_string())
        feature_sequences['transcript_id'].append(row[0])
        feature_sequences['gene_id'].append(row[5])
        feature_sequences['coverage_percentage'].append(general_per_base_coverage(feature_depth))
        feature_sequences['sequence'].append("".join(feature_seq))
        sum_count += sum(feature_depth)
        feature_sequences['coverage_per_base'].append(feature_depth)
    print(sum_count)
    return pl.DataFrame(feature_sequences)
'''
{'transcript_id': 'NR_026823_1', 'seqname': 'chr1', 'strand': ['-', '-', '-'], 'start': [205129, 205793, 206237], 'end': [205690, 205995, 206597], 'gene_id': 'FAM138D'}
'''

def get_cds_sequence(
    transcript_id: str, 
    gtf: GTFobject = None,
    dna: DNAStringSet = None
) -> DataFrame:
    '''
    '''    
    transcripts_data = gtf.get_cds_data(transcript_id)
    iter_df = transcripts_data.iterrows()
    cds_sequences = {
        'protein_id': [],
        'transcript_id': [],
        'gene_id': [],
        'sequence': [],
        'seqname': [],
        'start': [],
        'end': [],
        'strand': [],
        }
    for row in iter_df:
        cds_seq = []
        cds_start = []
        cds_end = []
        ref_seq = dna.get_sequence(f'chr{row[1]}')
        for index in range(len(row[2])):
            cds_start.append(row[2][index])
            cds_end.append(row[3][index])
            seq = ref_seq.get_sub_sequence(row[2][index], row[3][index])
            if row[4] == '-':
                seq = seq.reverse_complement()
            cds_seq.append(seq.get_as_string())
        cds_sequences['protein_id'].append(row[0])
        cds_sequences['transcript_id'].append(row[5])
        cds_sequences['gene_id'].append(row[6])
        cds_sequences['sequence'].append("".join(cds_seq))
        cds_sequences['seqname'].append(row[1])
        cds_sequences['start'].append(cds_start)
        cds_sequences['end'].append(cds_end)
        cds_sequences['strand'].append(row[4])

    return pl.DataFrame(cds_sequences)

def get_deseq2_matrix(
        bam_files: list[BAMhandler], 
        transcripts_ids: list[str],
        gtf: GTFobject,
        output_path: str
    ):
        transcripts_data = gtf.get_transcripts_data(transcripts_ids)

        count_matrix = {
            'transcript_id': [],
            'pseudo_ref': [],
        }

        for bam in bam_files:
            count_matrix[bam.path] = []
            count_matrix[f'{bam.path}.ratio'] = []
        
        iter_df = transcripts_data.iter_rows(named=True)
        for row in iter_df:
            count_matrix['transcript_id'].append(row['transcript_id'])
            total_mult = 1
            for bam in bam_files:
                bam_count = 0
                for index in range(len(row['start'])):
                    bam_count += bam.read_count(row['seqname'], row['start'][index], row['end'][index])
                total_mult *= bam_count
                count_matrix[bam.path].append(bam_count)
            count_matrix['pseudo_ref'].append(total_mult ** (1/len(bam_files)))
            for bam in bam_files:
                try:
                    ratio = count_matrix[f"{bam.path}"][-1]/count_matrix["pseudo_ref"][-1]
                except ZeroDivisionError:
                    ratio = 0
                count_matrix[f'{bam.path}.ratio'].append(ratio)
        df = DataFrame(count_matrix)
        for bam in bam_files:
            df = df.with_columns(
                [
                    (pl.col(f'{bam.path}') / df.select(f'{bam.path}.ratio').to_series().median()).alias(f'{bam.path}.norm'),
                ]
            )
        for bam in bam_files:
            new_df = df.select([
                pl.col('transcript_id'), 
                pl.col(f'{bam.path}').alias('raw_count'),
                pl.col(f'{bam.path}.norm').alias('normalized_count')
            ])
            new_df.write_csv(f'{bam.path}.deseq2')
        # df.write_csv(output_path)


def normalize_base_with_deseq2(transcript_ids: list[str], base_count_file: str, deseq2_matrix_file: str):
    feature_count_per_base = read_from_file(base_count_file, 'parquet')
    deseq2_matrix = read_from_file(deseq2_matrix_file, 'csv')

    transcript_norm_factor = deseq2_matrix.filter(
            (
                pl.col("transcript_id") == transcript_ids[0]
            ) 
        ).select(['raw_count','normalized_count']).to_dict()
    norm_factor = transcript_norm_factor['normalized_count'].item()/transcript_norm_factor['raw_count'].item()
    normalized_base_count = []
    present_transcript_ids = []
    for transcript in transcript_ids:
        raw_base_count = feature_count_per_base.filter(
            (
                pl.col("transcript_id") == transcript
            ) 
        ).select('coverage_per_base').to_series().to_list()
        if len(raw_base_count) > 0:
            raw_base_count = raw_base_count[0]
            try:
                normalized_base_count.append([base_count*norm_factor for base_count in raw_base_count])
            except ZeroDivisionError:
                normalized_base_count.append([0]*len(raw_base_count))
            present_transcript_ids.append(transcript)
    new_df = DataFrame({'transcript_id': present_transcript_ids, 'normalized_count':normalized_base_count})
    new_df = feature_count_per_base.join(new_df, on='transcript_id', how='inner')
    new_df.write_parquet('new_rpf.parquet')

