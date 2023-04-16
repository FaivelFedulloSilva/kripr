from collections import defaultdict
import os
import shutil
import pysam
from pysam import AlignmentFile
from dataclasses import dataclass


@dataclass
class BAMhandler:

    path: str
    bam_object: AlignmentFile = None

    def parse_bam(self, options='b'):
        if options == 'b':
            self.bam_object = pysam.AlignmentFile(self.path, 'rb')
        elif options == 's':
            self.bam_object = pysam.AlignmentFile(self.path, 'rs')


    def region_coverage(self, contig: str, start: int, end: int):
        coverage = pysam.depth(self.path, '-a', '-r', f"{contig}:{start}-{end}")
        coverage = coverage.split('\n')
        coverage_return = []
        for line in coverage[:-1]:
            splitted_line = line.split('\t')
            coverage_return.append([splitted_line[1], splitted_line[2]])
        return coverage_return
    
    def read_count(
        self, 
        gene: tuple[str,int,int] = None,
        contig: str = None,
        start: int = None,
        end: int = None,
        region: str= None
    ):
        bam_iter = self.bam_object.count(gene[0], gene[1], gene[2])
            
        return bam_iter

            
    def get_read_length_histogram(self) -> tuple[dict[int,int], int, int]:
        """
        Returns a tuple containing a histogram of read lengths, and the number of mapped and unmapped reads.

        Returns:
        tuple: A tuple containing the following elements:
            - A dictionary representing the histogram of read lengths, where the keys are the lengths and the values are the number of reads with that length.
            - An integer representing the total number of mapped reads.
            - An integer representing the total number of unmapped reads.
        """
        histogram = defaultdict(int)
        mapped_reads = 0
        unmapped_reads = 0
        bam = None
        if self.bam_object:
            bam = self.bam_object
        else:
            bam = pysam.AlignmentFile(self.path, 'rb')
        for read in bam.fetch(until_eof=True):
            histogram[len(read.query_sequence)] += 1
            if read.is_unmapped:
                unmapped_reads += 1
            else:
                mapped_reads += 1
        
        if not self.bam_object:
            bam.close()

        return histogram, mapped_reads, unmapped_reads

    def split_bam_by_chromosome(self, output_folder: str, overwrite: bool = False) -> None:
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        else:
            if overwrite:
                shutil.rmtree(output_folder)
                os.mkdir(output_folder)
            else:
                return
        print('chau')
        bam_file = pysam.AlignmentFile(self.path, 'rb')

        for chromosome in bam_file.references:
            print(chromosome)
            output_bam = pysam.AlignmentFile(f"{output_folder}/{chromosome}.bam", "wb", template=bam_file)

            for read in bam_file.fetch(chromosome):
                # write the read to the output BAM file for the current chromosome
                output_bam.write(read)
            # close the output BAM file for the current chromosome
            output_bam.close()
            # create an index for the output BAM file for the current chromosome
            pysam.index(f"{output_folder}/{chromosome}.bam")

    def filter_bam(
        self,
        output_path: str,
        filter_by_length: bool,
        filter_by_mapping_quality: bool,
        min_length: int,
        max_length: int,
        min_mapping_quality: int,
    ):
        # Define a function to filter reads based on length
        def filter_read_length(read, min_length: int, max_length: int) -> bool:
            return len(read.query_sequence) >= min_length and len(read.query_sequence) <= max_length

        # Define a function to filter reads based on mapping quality
        def filter_mapping_quality(read, min_mapping_quality: int) -> bool:
            return read.mapping_quality >= min_mapping_quality
        
        # Open BAM file for reading
        if self.bam_object:
            bamfile = self.bam_object
        else:
            bamfile = pysam.AlignmentFile(self.path, "rb")  
        filtered_reads = [read for read in self.bamfile]
        
        if filter_by_length:
            filtered_reads = [read for read in filtered_reads if filter_read_length(read, min_length, max_length)]
        
        if filter_by_mapping_quality:
            filtered_reads = [read for read in filtered_reads if filter_mapping_quality(read, min_mapping_quality)]

        filtered_bamfile = pysam.AlignmentFile(output_path, "wb", template=bamfile)
        for read in filtered_reads:
            filtered_bamfile.write(read)

        filtered_bamfile.close()  # Close filtered BAM file
        if not self.bam_object:
            bamfile.close()


    def close(self):
        self.bam_object.close()


# RPF_PATH = './Data/RPF_1/accepted_hits_01.bam'
# TOTAL_PATH = './Data/totalRNA_01/accepted_hits_01.bam'

# contig = "chr15"
# start = 30903852
# end = 30903887
# regio = f"{contig}:{start}-{end}"
# bam = BAMhandler(RPF_PATH, 'b')
# print(bam.region_coverage(contig, start, end))
# counted = bam._bam_object.count(region=regio)
# print(counted)
# count_coverage = bam._bam_object.pileup(contig, start, end)
# for i in count_coverage:
#     print('->', i)

    

if __name__ == '__main__':
    RPF_PATH = r"./../../RPF-Tesis/Data/BAMS/RPF/ZB/accepted_hits_10.bam"
    bam = BAMhandler(RPF_PATH, 'b').parse_bam()
    bam.do()