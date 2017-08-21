# Function to call vcf2maf for mitochondrial variants

import os, sys, pdb, numpy as np, scipy as sp, pandas as pd,time, argparse, pysam

# Parse necessary arguments
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("-f", "--datafile",type=str, help="full path to BAM fle")
parser.add_argument("-o","--outdir",type=str,help="directory for MAF files")
parser.add_argument("-h","--help",action='help', default=argparse.SUPPRESS,
                    help='Count reads for estimating mtDNA copy number from exome sequencing data.')

args = parser.parse_args()

# Read in the arguments
f = args.datafile
outdir = args.outdir

# Make sure the output directories are created
if not os.path.exists(outdir):
    os.makedirs(outdir)
    
MTcall = ' '.join(['samtools','view','-c','-q','30','-f','2','-F','1536',f,'MT','>',outdir+f.split('/')[-1] + '.MTcounts'])
nuclearcall = ' '.join(['samtools','view','-c','-q','30','-f','2','-F','1536', f,'>',outdir+f.split('/')[-1] + '.nuclearcounts'])

fullcall = MTcall + ';' + nuclearcall
print(fullcall)
os.system(fullcall)

    
    