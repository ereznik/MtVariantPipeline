# Function to call vcf2maf for mitochondrial variants

import os, sys, pdb, numpy as np, scipy as sp, pandas as pd,time, argparse, pysam

# Parse necessary arguments
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("-d", "--datadir",type=str, help="directory for BAM files")
parser.add_argument("-v", "--vcfdir", type=str, help="directory for intermediate VCF files")
parser.add_argument("-o","--outdir",type=str,help="directory for MAF files")
parser.add_argument("-b","--bamfiles",type=str,help="path to tab-delimited, no header file of tumor/normal bams. BAM names should not include the full path, just the name of the files in datadir to call variants on")
parser.add_argument("-q","--mapq",type=int,help="minimum mapping quality, default = 10",default = 10)
parser.add_argument("-Q","--baseq",type=int,help="minimum base quality, default = 10",default = 10)
parser.add_argument("-h","--help",action='help', default=argparse.SUPPRESS,
                    help='A simple variant calling and annotation pipeline for mitochondrial DNA variants. Accepts both individual BAM files and paired tumor/normal BAM files. Because of the hairiness of multiallelic calling, we keep all positions in the VCF and then filter to remove any positions without at least 10 reads supporting the putative variant. FIX TO DO THIS FOR BOTH NORMAL AND TUMOR SAMPLES! To send a call to bsub, try bsub -R "rusage[mem=16]" -M 32 -We 120 -W 4800 -e $HOME/work/mtimpact/scratch/ -o /home/reznik/work/mtimpact/scratch/ python MTvariantpipeline.py ... with the suitable options specified. Note that we use the CMO version of b37 (which seems to use rCRS), but is named HG19.')

args = parser.parse_args()

# Read in the arguments
datadir = args.datadir
vcfdir = args.vcfdir
outdir = args.outdir

# Set key parameters
minmapq = args.mapq
minbq = args.baseq

# Make sure the output directories are created
if not os.path.exists(vcfdir):
    os.makedirs(vcfdir)
if not os.path.exists(outdir):
    os.makedirs(outdir)
    
# Read in the annotaitons
trna = pd.read_csv('../data/MitoTIP_March2017.txt',header = 0,sep = '\t')
mitimpact = pd.read_csv('../data/MitImpact_db_2.7.txt',header = 0,sep = '\t',decimal = ',') # note that the decimal point here is indicated as a comma, mitimpact is funny

# Make the indices searchable for annotation for tRNA data
trna.index = [trna.at[item,'rCRS base'] + str(trna.at[item,'Position']) + trna.at[item,'Change'] if trna.at[item,'Change'] != 'del' else trna.at[item,'rCRS base'] + str(trna.at[item,'Position']) + trna.at[item,'Change'] + trna.at[item,'rCRS base'] for item in trna.index]

# For Mitimpact, consider only snps, and make data searchable
mitimpact = mitimpact[ mitimpact['Start'] == mitimpact['End'] ]
mitimpact.index = [mitimpact.at[item,'Ref'] + str(mitimpact.at[item,'Start']) + mitimpact.at[item,'Alt'] for item in mitimpact.index]

# Make sure there are no duplicate indices for the annotation
trna = trna[~trna.index.duplicated(keep='first')]
mitimpact = mitimpact[~mitimpact.index.duplicated(keep='first')]

# Indicate the columns to keep for trna and mitimpact when annotating
trna_cols = ['Predictive score']
mitimpact_cols = ['APOGEE_boost_mean_prob', 'Mitomap_Dec2016_Status', 'Mitomap_Dec2016_Disease']

# These are the "bad" mutations we automatically call pathogenic
badmuts = ['Nonsense_Mutation','Nonstop_Mutation','Frame_Shift_Del','Frame_Shift_Ins']

# Read in the paired BAM files
bamfiles = pd.read_csv(args.bamfiles,header = None,sep = '\t')
fs = bamfiles.ix[:,0]

# Also create a dataframe that stores the depth of MT coverage for each sample
mtcounts = pd.DataFrame( columns = ['MTCounts'] )

for ii in range(bamfiles.shape[0]):
    
    f = bamfiles.at[ii,0].strip() # remove leading and trailing whitespace
    
    print('Working on ' + f + '...')
    
    # Try to get the readcounts for the file
    try:
        mt = pysam.view('-c',datadir + f,'-q 10','-F 1536','MT')
        mtcounts.at[f,'MTCounts'] = int(mt)
    except:
        print('Error in getting read counts for ' + f + ', moving on...')
        continue
        
    # Check if we have a normal file
    if pd.isnull(bamfiles.at[ii,1]):
        normalflag = False
    else:
        normalbam = bamfiles.at[ii,1].strip() # remove leading and trailing whitespace
        normalflag = True
    
    # Part 1: Variant calling. Some samples have normal bam files, others do not. Do not confuse the two.
    
    if normalflag:
        
        # We have a normal bam
        countcall = ' '.join(["/opt/common/CentOS_6-dev/bin/current/samtools mpileup --region MT --count-orphans --no-BAQ --min-MQ ",str(minmapq), "--min-BQ", str(minbq), "--ignore-RG --excl-flags UNMAP,SECONDARY,QCFAIL,DUP --BCF --output-tags DP,AD,ADF,ADR --gap-frac 0.005 --tandem-qual 80 --fasta-ref /ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta", datadir  + f, datadir + normalbam + "| bcftools call --multiallelic-caller --ploidy GRCh37 --keep-alts ", ">", vcfdir + f + "_temp.vcf"])
    
        mafcall = "cmo_vcf2maf --version develop --input-vcf " + vcfdir + f + ".vcf " + "--output-maf " + outdir + f + ".maf --vcf-tumor-id " + datadir + f + " --vcf-normal-id " + datadir + normalbam + " --tumor-id " + f + " --normal-id " + normalbam
    
    else:
        # We don't have a normal bam
        print('We do not have a normal bam file for ' + f)
        
        countcall = ' '.join(["/opt/common/CentOS_6-dev/bin/current/samtools mpileup --region MT --count-orphans --no-BAQ --min-MQ ",str(minmapq), "--min-BQ", str(minbq), "--ignore-RG --excl-flags UNMAP,SECONDARY,QCFAIL,DUP --BCF --output-tags DP,AD,ADF,ADR --gap-frac 0.005 --tandem-qual 80 --fasta-ref /ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta", datadir + f, "| bcftools call --multiallelic-caller --ploidy GRCh37 --keep-alts ", ">", vcfdir + f + "_temp.vcf"])
    
        mafcall = "cmo_vcf2maf --version develop --input-vcf " + vcfdir + f + ".vcf " + "--output-maf " + outdir + f + ".maf --vcf-tumor-id " + datadir + f + " --tumor-id " + f
    
    
    # Make the VCF file
    #print(countcall)
    os.system(countcall)
    
    # Write comment lines to final VCF
    os.system( "cat " + vcfdir + f + "_temp.vcf | grep -F '#' > " + vcfdir + f + ".vcf")
    
    # Read in the VCF file and only keep rows with at least 3 reads supporting the variant
    
    # But first, count the number of comment lines
    skip_rows = 0
    with open(vcfdir + f + '_temp.vcf', 'r') as ftemp:
        for line in ftemp:
            if line.startswith('##'):
                skip_rows += 1
            else:
                break
    
    vcf = pd.read_csv( vcfdir + f + '_temp.vcf', header = 0,sep = '\t',skiprows = skip_rows )
    vcf = vcf[vcf['ALT'] != '.'] # drop any sites with no alterative variant
    #vcf['AltDepths_Tumor'] = [vcf.ix[item,9].split(':')[-1].split(',')[1:] for item in vcf.index]
    #vcf['MaxDepth_Tumor'] = [ np.max([int(item2) for item2 in item]) for item in vcf['AltDepths_Tumor'] ]
    #vcf = vcf[vcf['MaxDepth_Tumor'] >= 10]
    
    # Write out VCF, making sure to remove the extra allele depth columns
    #vcf = vcf.drop( ['AltDepths_Tumor','MaxDepth_Tumor'], axis = 1 )
    vcf.to_csv( vcfdir + f + ".vcf", sep = '\t', header = None, index = None, mode = 'a' )
    
    #print(mafcall)
    os.system(mafcall)
    
    #pdb.set_trace()
    
    ####################################################################################
    ####################################################################################
    
    # Part 2: Variant annotation. The general workflow is to annotate SNPs with tRNA and mitimpact data, and to assume that all frameshifts/nonsense are potentially pathogenic.
    
    # Read in the MAF file
    maf = pd.read_csv(outdir + f + '.maf',header = 0,sep = '\t',comment = '#')
    
    # Make a short name for each variant. 
    maf['ShortVariantID'] = maf['Reference_Allele'] + maf['Start_Position'].map(str) + maf['Tumor_Seq_Allele2']
    maf.ix[maf['Variant_Type']!='SNP','ShortVariantID'] = 'NA'
    
    # Add the annotations
    for col in trna_cols + mitimpact_cols:
        colorder = maf.columns.tolist()
        maf = pd.concat( [maf,pd.DataFrame( columns = [col] )] )
        
        # Keep column corder
        maf = maf[colorder + [col]]
        
    maf[trna_cols] = trna.ix[ maf[ 'ShortVariantID' ], trna_cols].reset_index()[trna_cols]
    maf[mitimpact_cols] = mitimpact.ix[ maf['ShortVariantID'],:].reset_index().ix[:, mitimpact_cols]
    
    # Add tumor and normal VAF
    maf['TumorVAF'] = maf['t_alt_count']/maf['t_depth']
    maf['NormalVAF'] = maf['n_alt_count']/maf['n_depth']
    
    ####################################################################################
    ####################################################################################
    #pdb.set_trace()
    # Part 3: Assign conditions for pathogenicity
    maf['Pathogenic_Reason'] = ''
    maf['Pathogenic_mtDNA_Variant'] = False
    
    # Anything with APOGEE score greater than 0.9
    maf.ix[maf['APOGEE_boost_mean_prob'] > 0.9,'Pathogenic_mtDNA_Variant'] = True
    maf.ix[maf['APOGEE_boost_mean_prob'] > 0.9,'Pathogenic_Reason'] = 'APOGEE'
    
    # Anything confirmed in MITOMAP
    maf.ix[ maf['Mitomap_Dec2016_Status'].isin(['Confirmed','Cfrm']), 'Pathogenic_mtDNA_Variant'] = True
    maf.ix[ maf['Mitomap_Dec2016_Status'].isin(['Confirmed','Cfrm']), 'Pathogenic_Reason'] = 'MITOMAP'
    
    # Anything with tRNA score greater than 18
    maf.ix[maf['Predictive score'] > 16.2,'Pathogenic_mtDNA_Variant'] = True
    maf.ix[maf['Predictive score'] > 16.2,'Pathogenic_Reason'] = 'tRNA Predictive Score'
    
    # Any frameshift/nonsense variants
    maf.ix[maf['Variant_Classification'].isin(badmuts),'Pathogenic_mtDNA_Variant'] = True
    maf.ix[maf['Variant_Classification'].isin(badmuts),'Pathogenic_Reason'] = 'Frameshift/Nonsense'
    
    ####################################################################################
    ####################################################################################
    
    # Part 3: Write out file
    maf.to_csv(outdir + f + '.maf',index = None,sep = '\t')

# Write out the counts file as well
mtcounts.to_csv(outdir + args.bamfiles.split('/')[-1] + 'Counts.txt')
    