# Script to split a single cell RNA BAM that is aggregated across all cells into distinct files per cell. Data is stored in /ifs/e63data/leslielab/pelossof/public/sc/

import pysam, os, sys, numpy as np, pandas as pd, pdb
from collections import Counter

b = str(sys.argv[1])
bamfile = pysam.AlignmentFile(b,'rb')

datadir = '/ifs/work/schultz/reznik/colorectalsc/'
vcfdir = '/ifs/e63data/schultzlab/reznik/mtimpact/scratch/'
outdir = '/ifs/e63data/schultzlab/reznik/mtimpact/mafs_colorectalsplit/'

# first, iterate through all the mitochondrial reads and get the distinct cell barcodes
print('Getting unique cells...')
cellids = list()
for read in bamfile.fetch('MT',1,16569):
    
    # for now, keep only high quality read
    if read.mapping_quality != 255:
        continue
    cellids.append( read.query_name.split(':')[1] ) 
    
# identify ids that have sufficiently sequencing depth and make unique bam files for them
cellcounts = Counter(cellids)
highcells = [item[0] for item in cellcounts.items() if cellcounts[item[0]] > 1000]

# Write out counts
pd.DataFrame(cellcounts.values(),index = cellcounts.keys()).to_csv(datadir + b.split('/')[-1] + 'CellCounts_mtDNA.csv')
#pdb.set_trace()
#highcells = ['CAAGTTGTCGATAGAA'] # for testing

print('Going through each cell with sufficient reads')
for uqcell in highcells:
    umidict = dict()
    dropctr = 0
    print(uqcell)
    
    splitbamname =  b.split('/')[-1].split('_')[0] + '_' + uqcell + '.bam'
    bamout = pysam.AlignmentFile(datadir + splitbamname, 'wb', template=bamfile)
    
    for read in bamfile.fetch('MT',1,16569):
    
        # for now, keep only high quality read
        if read.mapping_quality != 255 or read.query_name.split(':')[1] != uqcell:
            continue
        else:
        
            # Save the UMI
            umi = read.query_name.split(':')[2]
            startpos = read.reference_start
            if umi not in umidict.keys():
                bamout.write(read)
                if startpos is not None:
                    umidict[umi] = [startpos]
            else:
                #print('Checking umi ' + str(umi))
                
                # Get prior entries
                priorumi = umidict[umi]
                
                # Check if the reads sharing UMIs are near each other
                closestval = [item - startpos for item in priorumi]
                if np.min( closestval ) < 10:
                    #print('Found a duplicate UMI that needs to be collapsed:' + umi)
                    dropctr += 1
                    continue
                else:
                    bamout.write(read)
                    if startpos is not None:
                        umidict[umi] = umidict[umi].append(startpos)
                    
    print('This many reads were dropped from this sample: ' + str(dropctr))
    bamout.close()
    
    # Index
    os.system('samtools index ' + datadir + splitbamname)
    
    # Now call the variant-calling pipeline
    vcall = ' '.join(['python /home/reznik/work/mtimpact/analysis/MTvariantpipeline.py', '-d',datadir, '-v',vcfdir, '-o',outdir, '-b',splitbamname, '-g', 'GRCh37','-q', '0', '-Q', '0','--strand','1'])
    
    #os.system(vcall)