import math
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from collections import defaultdict, OrderedDict
from multiprocessing import Pool
import os
from random import shuffle


import pysam
from kripr import get_cds_sequence, normalize_base_with_deseq2, get_deseq2_matrix, BioString, DNAString, DNAStringSet, GTFobject, BAMhandler, FastaHandler, GTFhandler, save_to_file, get_feature_sequence, get_feature_sequence_parallel
import numpy as np
import polars as pl
from polars import DataFrame
import time

class Timer:
    def __init__(self, name) -> None:
        self.name = name
        self.start = None

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self,exception_type, exception_value, exception_traceback):
        print(f'{self.name} has ended in {time.time() - self.start} seconds')



def use_case_get_transcript_data_to_exploratory_analisis():
    REF_PATH = '../RPF-Tesis/Data/reference/hg38.fa'
    GTF_PATH = '../RPF-Tesis/Data/genesFiltrada.gtf'
    RNA_BAMS = [
        "../RPF-Tesis/Data/BAMS/RNA/ZB/accepted_hits_06.bam",
        "../RPF-Tesis/Data/BAMS/RNA/ZB/accepted_hits_10.bam",
        "../RPF-Tesis/Data/BAMS/RNA/ZB/accepted_hits_11.bam",
    ]
    RPF_BAMS = [
        "../RPF-Tesis/Data/BAMS/RPF/ZB/accepted_hits_06.bam",
        "../RPF-Tesis/Data/BAMS/RPF/ZB/accepted_hits_10.bam",
        "../RPF-Tesis/Data/BAMS/RPF/ZB/accepted_hits_11.bam",
    ]

    ref = FastaHandler(REF_PATH).get_reference()

    gtf = GTFhandler(GTF_PATH).get_gtf_object()

    RPF_BAMS_HANDLERS = []
    for rpf in RPF_BAMS:
        RPF_BAMS_HANDLERS.append(BAMhandler(rpf))
        RPF_BAMS_HANDLERS[-1].parse_bam()
    RNA_BAMS_HANDLERS = []
    for rna in RNA_BAMS:
        RNA_BAMS_HANDLERS.append(BAMhandler(rna))
        RNA_BAMS_HANDLERS[-1].parse_bam()


    transcript_ids = gtf.get_transcripts_ids()
    rpf_bam: BAMhandler = RPF_BAMS_HANDLERS[0]
    rna_bam: BAMhandler = RNA_BAMS_HANDLERS[0]
    
    with Timer("split bam file"):
        rpf_bam.split_bam_by_chromosome('./data/rpf_splitted_bam')

    quantity = 100
    # shuffle(transcript_ids)
    start = time.time()
    dt_rna = get_feature_sequence(transcript_ids[:quantity],gtf, ref, bam=rpf_bam)
    end = time.time()
    print(end - start)
    print(dt_rna)

    start = time.time()
    dt_rna_parallel = get_feature_sequence_parallel(20, transcript_ids[:quantity],gtf, ref, bam=rpf_bam, bams_folder='./data/rpf_splitted_bam')
    end = time.time()
    print(end - start) 
    print(dt_rna)


def map_coordinates(arr: list[list[int]]):
    new_arr = []
    for l in arr:
        new_arr += list(range(l[0], l[1]+1))
    return new_arr

def complete_coordinates(coordinates: list[int], counts: list[int], frame: int, cds_len: int):
    index = frame
    i = 0
    new_coordinates = [coordinate*3 + frame for coordinate in range(0, (cds_len//3) + 1)]
    new_counts = []
    print(len(coordinates))
    for index in new_coordinates:
        if i >= len(coordinates):
            i += 1
            new_counts.append(0)
        else:
            if coordinates[i] == index:
                new_counts.append(counts[i])
                i += 1
            else:
                new_counts.append(0)
    return new_coordinates, new_counts

def get_cds_per_codon_counts(
    transcript_ids: list[str],
    shift: dict[str,int],
    bam: BAMhandler,
    reference: FastaHandler,
    gtf: GTFobject,
    min_expression: int =1,
    min_len: int = 20,
    max_len: int = 30,
):
    def shift_by_length(read_length, shift):
        shift = {
            '20':	10,
            '21':	10,
            '22':	10,
            '23':	10,
            '24':	11,
            '25':	11,
            '26':	11,
            '27':	11,
            '28':	12,
            '29':	12,
            '30':	14,
            '31':	14,
            'default':	14
        }
        read_length_str = str(read_length)
        if read_length > 31:
            read_length_str = 'default'
        return shift[read_length_str]
    bam = bam.bam_object
    cdss = {
        'protein_ids': [],
        'frame': [],
        'p_site_start': [],
        'read_count': []
    }
    for t_idx in range(len(transcript_ids)):
        seq_dt = get_cds_sequence(dna=reference.get_reference(), gtf=gtf, transcript_id=transcript_ids[t_idx])
        strand = '+'
        protein_id = ''
        cds_sections = []
        cds_sequence = ''
        iter_df = seq_dt.iter_rows()
        for row in iter_df:
            cds_sequence = row[3]
            protein_id = row[0]
            strand = row[7][0]
            for i in range(len(row[5])):
                cds_sections.append([row[4],row[5][i], row[6][i]])
        if len(cds_sections) == 0:
            continue

        cds_sections.sort(key = lambda x: x[1])

        full_mapped_cds = map_coordinates([[x[1],x[2]] for x in cds_sections])
        
        if strand == '-':
            full_mapped_cds.reverse()

        counts2 = {
            'read_name': [],
            'shift': [],
            'p_site_start': [],
            'frame': []
        } 

        per_frame = {
            20: [0,0,0],
            21: [0,0,0],
            22: [0,0,0],
            23: [0,0,0],
            24: [0,0,0],
            25: [0,0,0],
            26: [0,0,0],
            27: [0,0,0],
            28: [0,0,0],
            29: [0,0,0],
            30: [0,0,0],
            31: [0,0,0],
            32: [0,0,0],
            33: [0,0,0],
            34: [0,0,0],
            35: [0,0,0],
            36: [0,0,0],
            37: [0,0,0],
            38: [0,0,0],
            39: [0,0,0],
            40: [0,0,0],
        }
        if strand == '+':
            reads = bam.fetch(region=f'chr{cds_sections[0][0]}:{full_mapped_cds[0]}-{full_mapped_cds[-1]}')
        else:
            reads = bam.fetch(region=f'chr{cds_sections[0][0]}:{full_mapped_cds[-1]}-{full_mapped_cds[0]}')

        for r in reads:
            if r.query_length >= min_len and r.query_length <= max_len:
                read_name = r.query_name
                shift = shift_by_length(r.query_length, shift)  
                try:
                    p_site_start = full_mapped_cds.index(r.reference_start + 1 + shift)
                except:
                    #TODO: Los reads descartados pasan a una "bolsa" que se retorna en caso de que algun parametro asi lo indique
                    continue
                counts2['read_name'].append(read_name)
                counts2['shift'].append(shift)
                counts2['p_site_start'].append(p_site_start)
                counts2['frame'].append((p_site_start)%3)
                print(f'name {r.query_name} : readStart(0-base) {r.reference_start} : p-site-position(0-base) {p_site_start} : length {r.query_length} : shift {shift} : frame {p_site_start%3} : query_alignment_start {r.query_alignment_start}')
                per_frame[r.query_length] = per_frame.get(r.query_length, [0,0,0])
                per_frame[r.query_length][(p_site_start)%3] += 1
        # print(len(counts2['read_name']))
        input()
        if len(counts2['read_name']) >= min_expression:
            df = DataFrame(counts2)
            print(protein_id)
            print(strand)
            print('length\tframe 0\tframe 1\tframe 2')
            for l,f in per_frame.items():
                print(f'{l}\t{f[0]}\t{f[1]}\t{f[2]}')
            frame0_percent = sum([f[0] for f in per_frame.values()])/len(counts2['read_name'])
            frame1_percent = sum([f[1] for f in per_frame.values()])/len(counts2['read_name'])
            frame2_percent = sum([f[2] for f in per_frame.values()])/len(counts2['read_name'])

            # print(f' \t{frame0_percent:.2}%\t{frame1_percent:.2}%\t{frame2_percent:.2}%')
            # print(df)
            if len(df.select('read_name').to_series().to_list()) > 0:
                df_0 = df.filter(pl.col('frame') == 0)
                df_1 = df.filter(pl.col('frame') == 1)
                df_2 = df.filter(pl.col('frame') == 2)

                df_0 = (
                    df_0.lazy()
                    .groupby("p_site_start")
                    .agg(
                        [
                            pl.col('read_name').count()
                        ]
                    )
                    .sort("p_site_start")
                ).collect()

                df_0 = df_0.to_dict()
                cdss['protein_ids'].append(protein_id)
                cdss['frame'].append(0)
                new_starts, new_counts = complete_coordinates(df_0['p_site_start'], df_0['read_name'], 0, len(full_mapped_cds))
                cdss['p_site_start'].append(new_starts)
                cdss['read_count'].append(new_counts)


                df_1 = (
                    df_1.lazy()
                    .groupby("p_site_start")
                    .agg(
                        [
                            pl.col('read_name').count()
                        ]
                    )
                    .sort("p_site_start")
                ).collect()

                df_1 = df_1.to_dict()
                cdss['protein_ids'].append(protein_id)
                cdss['frame'].append(1)
                new_starts, new_counts = complete_coordinates(df_1['p_site_start'], df_1['read_name'], 1, len(full_mapped_cds))
                cdss['p_site_start'].append(new_starts)
                cdss['read_count'].append(new_counts)

                df_2 = (
                    df_2.lazy()
                    .groupby("p_site_start")
                    .agg(
                        [
                            pl.col('read_name').count()
                        ]
                    )
                    .sort("p_site_start")
                ).collect()

                df_2 = df_2.to_dict()
                cdss['protein_ids'].append(protein_id)
                cdss['frame'].append(2)
                new_starts, new_counts = complete_coordinates(df_2['p_site_start'], df_2['read_name'], 2, len(full_mapped_cds))
                cdss['p_site_start'].append(new_starts)
                cdss['read_count'].append(new_counts)


                # print(df_0)
                # print(f"|{' END ':=^80}|")
                # input()
                dt_cdss = DataFrame(cdss)
        # print('--------------END----------------')
    print(dt_cdss)
    dt_cdss.write_json('cds_by_frame.json')

if __name__ == '__main__':

    CORESPONDANCE_PATH = './coorespondance.csv'


    DATA_PATH = './data/'
    RNA_PATH = './data/rna.parquet'
    RPF_PATH = './data/rpf.parquet'
    
    REF_PATH = '../RPF-Tesis/Data/reference/hg38.fa'
    GTF_PATH = '../RPF-Tesis/Data/cds_only.gtf'
    RNA_BAMS = [
        "../RPF-Tesis/Data/BAMS/RNA/D_PLUS/accepted_hits_02.bam",
        "../RPF-Tesis/Data/BAMS/RNA/D_PLUS/accepted_hits_05.bam",
        "../RPF-Tesis/Data/BAMS/RNA/D_PLUS/accepted_hits_09.bam",
    ]
    RPF_BAMS = [
        "../RPF-Tesis/Data/BAMS/RPF/D_PLUS/accepted_hits_02.bam",
        "../RPF-Tesis/Data/BAMS/RPF/D_PLUS/accepted_hits_05.bam",
        "../RPF-Tesis/Data/BAMS/RPF/D_PLUS/accepted_hits_09.bam",
    ]
    
    



    RPF_BAMS_HANDLERS: list[BAMhandler] = []
    for rpf in RPF_BAMS:
        RPF_BAMS_HANDLERS.append(BAMhandler(rpf))
        RPF_BAMS_HANDLERS[-1].parse_bam()
    RNA_BAMS_HANDLERS: list[BAMhandler] = []
    for rna in RNA_BAMS:
        RNA_BAMS_HANDLERS.append(BAMhandler(rna))
        RNA_BAMS_HANDLERS[-1].parse_bam()

    ref = FastaHandler(REF_PATH)

    gtf: GTFobject = GTFhandler(GTF_PATH).get_gtf_object()
    print(gtf.get_gene_name('ENSG00000071889', 'gene_id', CORESPONDANCE_PATH))
    input()
    '''
    1- Obtener CDS pieces para un transcripto 
    2- Unir los CDS pieces y obtener la referencia
    '''
    shift = {
        '20':	10,
        '21':	10,
        '22':	10,
        '23':	10,
        '24':	11,
        '25':	11,
        '26':	11,
        '27':	11,
        '28':	12,
        '29':	12,
        '30':	14,
        '31':	14,
        'default':	14
    }
    # transcript_ids = ['ENST00000641399','ENST00000488573','ENST00000482101','ENST00000473497','ENST00000456483','ENST00000446550','ENST00000443427','ENST00000440973','ENST00000427531','ENST00000415488','ENST00000406599','ENST00000404742','ENST00000338799','ENST00000206249','ENST00000643001','ENST00000481743','ENST00000461472','ENST00000379328','ENST00000346208','ENST00000652686','ENST00000514699','ENST00000510170','ENST00000508760','ENST00000505058','ENST00000504572','ENST00000504336','ENST00000503701','ENST00000503201','ENST00000502892','ENST00000502500','ENST00000424646','ENST00000415690','ENST00000394466','ENST00000394464','ENST00000343796','ENST00000231509','ENST00000557538','ENST00000557446','ENST00000557206','ENST00000556827','ENST00000556237','ENST00000555014','ENST00000553999','ENST00000547430','ENST00000539097','ENST00000394997','ENST00000337138','ENST00000323441','ENST00000699988','ENST00000673719','ENST00000652095','ENST00000651694','ENST00000651595','ENST00000651158','ENST00000650658','ENST00000620651','ENST00000614754','ENST00000477468','ENST00000466278','ENST00000463769','ENST00000422851','ENST00000407029','ENST00000402630','ENST00000402042','ENST00000355630','ENST00000265354','ENST00000602216','ENST00000600238','ENST00000600147','ENST00000599898','ENST00000599870','ENST00000597648','ENST00000222247','ENST00000493317','ENST00000485756','ENST00000483765','ENST00000473558','ENST00000466550','ENST00000461690','ENST00000326092','ENST00000319826','ENST00000311549','ENST00000272274']
    transcript_ids = ['ENST00000222247']
    get_cds_per_codon_counts(transcript_ids, shift, RPF_BAMS_HANDLERS[2], ref,gtf)
    input('done')
    dt_rpf_parallel = get_feature_sequence(transcript_ids,gtf, ref.get_reference(), bam=RPF_BAMS_HANDLERS[2])
    input('done')
    dt_rna_parallel = get_feature_sequence(transcript_ids,gtf, ref.get_reference(), bam=RNA_BAMS_HANDLERS[2])
    input('done')
    save_to_file(dt_rpf_parallel, 'rpf_tam')
    save_to_file(dt_rna_parallel, 'rna_tam')


    # transcript_ids = gtf.get_transcripts_ids()#[104:107]
    
            

    # get_deseq2_matrix(RNA_BAMS_HANDLERS, transcript_ids, gtf, './output.csv')

    # normalize_base_with_deseq2(transcript_ids, RNA_PATH, '/home/faivel/Documents/Projects/Tesis/RPF-Tesis/Data/BAMS/RNA/ZB/accepted_hits_06.bam.deseq2')

    
    # transcripts_data = gtf.get_transcripts_data(transcript_ids)

    # '''
    # 100 transcripts
    #     stdout: 0.003507375717163086
    #     file: 0.0002105236053466797

    # 1000 transcripts
    #     stdout: 0.031742095947265625
    #     file: 0.002426624298095703 
    

    # '''


    

    # # Define a function to count the reads for a single transcript and BAM file
    # def count_reads(bam_file, transcripts):
    #     bam = pysam.AlignmentFile(bam_file, "rb")
    #     counts = {}
    #     for transcript in transcripts:
    #         count = 0
    #         for i in range(len(transcript['start'])):
    #             count += bam.count(transcript['seqname'], transcript['start'][i], transcript['end'][i])
    #         counts[transcript['transcript_id']] = count
    #     bam.close()
    #     return (bam_file, transcript, counts)

    # def get_matrix(bam_files, chunk, transcript_size):
    #     for bam in bam_files:
    #         split_bam_by_chromosome(bam, f'{bam}.folder', False)

    #     pre_list = [(f'{bam_file}.folder/{transcript["seqname"]}.bam', transcript) for bam_file in bam_files for transcript in transcripts[:transcript_size]]
    #     post_list = defaultdict(list)

    #     for row in pre_list:
    #         post_list[row[0]].append(row[1])
        
        
    #     # Create a pool of worker processes and map the count_reads function over the transcript and BAM file pairs
    #     pool = Pool(processes=os.cpu_count() - 1)
    #     if chunk == 0:
    #         results = pool.starmap(count_reads, list(post_list.items()))
    #     else:
    #         results = pool.starmap(count_reads, list(post_list.items()), chunksize=chunk)
    #     pool.close()
    #     pool.join()

    #     print(results)
    #     print()
    #     print()
    #     # print(results)

    #     # # Create a dictionary to store the read counts for each transcript and each BAM file
    #     # counts = {}
    #     # for bam_file in bam_files:
    #     #     counts[bam_file] = {}
    #     # for bam_file, transcript, count in results:
    #     #     counts[bam_file][transcript['transcript_id']] = count

    #     # # Print the count matrix
    #     # print("\t".join(["Transcript"] + bam_files))
    #     # for transcript in transcripts:
    #     #     row = [transcript['transcript_id']]
    #     #     for bam_file in bam_files:
    #     #         row.append(str(counts[bam_file][transcript['transcript_id']]))
    #     #     print("\t".join(row))

    # # Define the list of BAM files to process
    # bam_files = RNA_BAMS

    # # Define the list of transcripts to count reads for
    # transcripts = transcripts_data.to_dicts()
    # for chunk in [0,10]:
    #     for size in [100, 1000]:
    #         with Timer(f'{chunk} chunks and {size} transcripts'):
    #             get_matrix(bam_files, chunk, size)



    '''
    Funcion. El problema esta en que te devuelve el conteo normalizado por transcripto. La deteccion de pausas se realiza por base. Como deberia de utilizar esta informacion?
    opcion 1- Quedarme unicamente con el factor de normalizacion y utilizarlo para normalizar el conteo de base
    opcion 2 - No se realiza analisis de pausa normalizado con deseq2
    '''
    

    # def get_deseq2_matrix(
    #     bam_files: list[BAMhandler], 
    #     transcripts_ids: list[str],
    #     gtf: GTFobject,
    #     output_path: str
    # ):
    #     transcripts_data = gtf.get_transcripts_data(transcripts_ids)

    #     count_matrix = {
    #         'transcript_id': [],
    #         'pseudo_ref': [],
    #     }

    #     for bam in bam_files:
    #         count_matrix[bam.path] = []
    #         count_matrix[f'{bam.path}.ratio'] = []
        
    #     iter_df = transcripts_data.iter_rows(named=True)
    #     for row in iter_df:
    #         count_matrix['transcript_id'].append(row['transcript_id'])
    #         total_mult = 1
    #         for bam in bam_files:
    #             bam_count = 0
    #             for index in range(len(row['start'])):
    #                 bam_count += bam.read_count([row['seqname'], row['start'][index], row['end'][index]])
    #             total_mult *= bam_count
    #             count_matrix[bam.path].append(bam_count)
    #         count_matrix['pseudo_ref'].append(total_mult ** (1/len(bam_files)))
    #         for bam in bam_files:
    #             try:
    #                 ratio = count_matrix[f"{bam.path}"][-1]/count_matrix["pseudo_ref"][-1]
    #             except ZeroDivisionError:
    #                 ratio = 0
    #             count_matrix[f'{bam.path}.ratio'].append(ratio)
    #     df = DataFrame(count_matrix)
    #     for bam in bam_files:
    #         print(bam, df.select(f'{bam.path}.ratio').to_series().median())
    #         df = df.with_columns(
    #             [
    #                 (pl.col(f'{bam.path}') / df.select(f'{bam.path}.ratio').to_series().median()).alias(f'{bam.path}.norm'),
    #             ]
    #         )
    #     df.write_csv(output_path)
    

    # transcript_ids = gtf.get_transcripts_ids()
    # bam_files = RPF_BAMS_HANDLERS

    # get_deseq2_matrix(bam_files, transcript_ids, gtf, './output.csv')
        # count_matrix = {
        #     'transcript_id': [],
        #     'BAM1': [],
        #     'BAM2': [],
        #     'BAM3': [],
        #     'pseudo_ref': [],
        #     'BAM1_ratio': [],
        #     'BAM2_ratio': [],
        #     'BAM3_ratio': [],
        # }
        
        # iter_df = transcripts_data.iterrows(named=True)
        # for row in iter_df:
        #     count_matrix['transcript_id'].append(row['transcript_id'])
        #     total_mult = 1
        #     for i in range(len(RPF_BAMS_HANDLERS)):
        #         bam_count = 0
        #         for index in range(len(row['start'])):
        #             bam_count += RPF_BAMS_HANDLERS[i].read_count([row['seqname'], row['start'][index], row['end'][index]])
        #         total_mult *= bam_count
        #         count_matrix[f"BAM{i+1}"].append(bam_count)
        #     count_matrix['pseudo_ref'].append(total_mult ** (1/len(RPF_BAMS_HANDLERS)))
        #     for i in range(len(RPF_BAMS_HANDLERS)):
        #         try:
        #             ratio = count_matrix[f"BAM{i+1}"][-1]/count_matrix["pseudo_ref"][-1]
        #         except ZeroDivisionError:
        #             ratio = 0
        #         count_matrix[f"BAM{i+1}_ratio"].append(ratio)
        # df = DataFrame(count_matrix)
        # new_df = df.with_columns(
        #     [
        #         (pl.col('BAM1') / df.select('BAM1_ratio').to_series().median()).alias("new_BAM1"),
        #         (pl.col('BAM2') / df.select('BAM2_ratio').to_series().median()).alias("new_BAM2"),
        #         (pl.col('BAM3') / df.select('BAM3_ratio').to_series().median()).alias("new_BAM3"),

        #     ]
        # )
        # new_df.write_csv("martix.csv")
        # bam1_ratio = df.select('BAM1_ratio').to_series().to_list()
        # print(bam1_ratio)
        # print()
        # print(np.median(np.array(bam1_ratio)))
        # print()
        # print()

        # bam2_ratio = df.select('BAM2_ratio').to_series().to_list()
        # print(bam2_ratio)
        # print()
        # print(np.median(np.array(bam2_ratio)))
        # print()
        # print()
        
        # bam3_ratio = df.select('BAM3_ratio').to_series().to_list()
        # print(bam3_ratio)
        # print()
        # print(np.median(np.array(bam3_ratio)))


    # rpf_bam: BAMhandler = RPF_BAMS_HANDLERS[0]
    # rna_bam: BAMhandler = RNA_BAMS_HANDLERS[0]
    
    # with Timer("split bam file"):
    #     split_bam_by_chromosome(rpf_bam.path, './data/rpf_splitted_bam')

    # quantity = 100
    # # shuffle(transcript_ids)
    # start = time.time()
    # dt_rna = get_feature_sequence(transcript_ids[:quantity],gtf, ref, bam=rpf_bam)
    # end = time.time()
    # print(end - start)
    # print(dt_rna)

    # start = time.time()
    # dt_rna_parallel = get_feature_sequence_parallel(20, transcript_ids[:quantity],gtf, ref, bam=rpf_bam, bams_folder='./data/rpf_splitted_bam')
    # end = time.time()
    # print(end - start) 
    # print(dt_rna)



    '''
    if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

    BiocManager::install("DESeq2")
    '''
    '''
    library(DESeq2)
    '''


    '''
    dds.rpfs <- DESeqDataSetFromMatrix(countData = cts.rpfs[,c(1,3:6)],
                              colData = coldata[c(1,3:6),],
                              design = ~ condition)'''