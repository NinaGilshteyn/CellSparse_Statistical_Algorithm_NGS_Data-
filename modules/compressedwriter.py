
# CellSparse: A Statistical Algorithm for Detecting Upregulated Genes in scRNA-seq >

#**Author:** Nina Gilshteyn, M.S.
#**Affiliation:** Deeds Lab, UCLA — Department of Integrative Biology & Physiology
#**Status:** Manuscript in preparation — *"The Biology is in The Tails: Skewness and>
#**Profile:** [ninagilshteyn.github.io](https://ninagilshteyn.github.io) · [github.c>



import sys
import os
import gzip
import csv

def compressed_write(listt, output_file):
#write the transformed data to a compressed csv file using gzip
    with gzip.open(output_file, 'wt', newline = '') as f: #OUTPUt file is the path and 'wt' is writing in text mode (instead of reading binary mode) 
        csv_writer = csv.writer(f)
        csv_writer.writerows(listt.values.tolist())
