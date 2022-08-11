import re
import sys
from itertools import chain

sam_file = sys.argv[2] + "/" + sys.argv[1] + ".final.sam"
insertion_site_file = sys.argv[2] + "/" + sys.argv[1] + ".insertion_site.txt"

strand = []
with open(sam_file, 'r') as sam_file, open(insertion_site_file, 'w') as insertion_site_file:
    for line in sam_file:
        samflag = line.split()[1]
        if samflag == "99":
            strand = "+"
        elif samflag == "83":
            strand = "-"
        chr_number = line.split()[2]
        position = int(line.split()[3])
        CIGAR = line.split()[5]
        mate_info = line.split()[11]
        TLEN = line.split()[8]
        read_ID= line.split()[0]
        insertion_site_start = position
        insertion_site_end = position

        digits = list(map(int, re.findall(r'\d+', CIGAR)))
        read_length = sum(digits)

        if samflag == "83":
            insertion_site_start = insertion_site_start + read_length -2
            insertion_site_end = insertion_site_end + read_length -2
        elif samflag == "99":
            insertion_site_start = insertion_site_start
            insertion_site_end = insertion_site_end

        result = [chr_number, str(insertion_site_start), str(insertion_site_end), str(CIGAR), str(mate_info), str(TLEN), str(read_ID), str(strand)]
        insertion_site_file.write('\t'.join(result)+'\n')
