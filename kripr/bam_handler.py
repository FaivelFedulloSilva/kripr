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
        """
        Parses a BAM file and returns a pysam.AlignmentFile object.

        Args:
            options: The mode to open the BAM file in. Can be 'b' for binary, 
            if a BAM file is provided or 's' for text, if a SAM file is provided.

        Returns:
            A pysam.AlignmentFile object.

        Raises:
            ValueError: If the options argument is not 'b' or 's'.

        """
        if options not in ('b', 's'):
            raise ValueError('Options must be "b" or "s"')

        self.bam_object = pysam.AlignmentFile(self.path, f'r{options}')


    def region_coverage(self, contig: str, start: int, end: int):
        """
        Returns the coverage of a region in a BAM file. The region is represented 
        as 1-based close intervals

        Args:
            contig: The name of the contig.
            start: The start position of the region.
            end: The end position of the region.

        Returns:
            A list of list, where each internal list contains the position and 
            coverage of a base in the region.

        """

        coverage = pysam.depth(self.path, '-a', '-r', f"{contig}:{start}-{end}")
        coverage = coverage.split('\n')
        coverage_return = []
        for line in coverage[:-1]:
            splitted_line = line.split('\t')
            coverage_return.append([splitted_line[1], splitted_line[2]])
        return coverage_return
    

    def read_count(
        self, 
        contig: str = None,
        start: int = None,
        end: int = None,
        region: str= None
    ):
        """
        Returns the number of reads in the BAM file that overlap the given region. 
        If region is not provided, then the contig, start and end shall be. 
        Region is considered as 1-based close interval, while contig, start, end 
        is considered as 0-based half open interval. 

        Args:
            contig: The name of the contig.
            start: The start position of the region.
            end: The end position of the region.
            region: A Samtools compatible string representing a genome region 

        Returns:
            A list of list, where each internal list contains the position and coverage of a base in the region.

        """
        if region:
            bam_iter = self.bam_object.count(region)
        else:
            bam_iter = self.bam_object.count(contig, start, end)
            
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
        """
        Splits a BAM file into separate BAM files, one for each chromosome.

        Args
        ----------
        output_folder : str
            The output folder to write the split BAM files to.
        overwrite : bool, optional
            If True, the output folder will be overwritten if it already exists. 
            If False, an error will be raised if the output folder already exists.
        """
        
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        else:
            if overwrite:
                shutil.rmtree(output_folder)
                os.mkdir(output_folder)
            else:
                return
        bam_file = pysam.AlignmentFile(self.path, 'rb')

        for chromosome in bam_file.references:
            output_bam = pysam.AlignmentFile(f"{output_folder}/{chromosome}.bam", "wb", template=bam_file)

            for read in bam_file.fetch(chromosome):
                output_bam.write(read)
            output_bam.close()
            pysam.index(f"{output_folder}/{chromosome}.bam")

    # TODO: Si se setea filter by length pero no se da un minimo y un maximo explota. Acomodar 
    def filter_bam(
        self,
        output_path: str,
        filter_by_length: bool,
        filter_by_mapping_quality: bool,
        min_length: int,
        max_length: int,
        min_mapping_quality: int,
    ):
        """
        Filters a BAM file by length and/or mapping quality.

        Parameters
        ----------
        output_path : str
            The path to the output BAM file.
        filter_by_length : bool, optional
            If True, reads will be filtered by length.
        filter_by_mapping_quality : bool, optional
            If True, reads will be filtered by mapping quality.
        min_length : int, optional
            The minimum length of a read to be included in the output BAM file.
        max_length : int, optional
            The maximum length of a read to be included in the output BAM file.
        min_mapping_quality : int, optional
            The minimum mapping quality of a read to be included in the output BAM file.

        Returns
        -------
        None

        """

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