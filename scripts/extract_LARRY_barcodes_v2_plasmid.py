# grep a LARRY-specific sequence
# Ask Alejo whether they have a pipeline to call lineage barcodes from FASTQ files, otherwise DIY
import os, pickle, gzip
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Extract the LARRY barcode.')
parser.add_argument('--fastq', help='The fastq with the barcodes')
parser.add_argument('--barcode', help='The cell barcode file')
parser.add_argument('--output', help='The output directory')
parser.add_argument('--name', help='The sample name')

args = parser.parse_args()

successes = []
no_fastq = []
no_abundant_bcs = []
if not os.path.exists(args.output + '/LARRY_tmp'):
    os.mkdir(args.output + '/LARRY_tmp')
    
# check that there is at least 1 *.fastq.sorted.fastq.gz file
sorted_fastqs = [args.fastq]

# check that there is an abundant bcs file
abundant_bcs = pd.read_csv(args.barcode, sep = '\t')
abundant_bcs = abundant_bcs.index.values.tolist()
    
# first, grep out lines with the barcode sequence...
gfp_seq = 'gctaggagagaccatatgggatccgat'.upper()

for i,ff in enumerate(sorted_fastqs):
    print('Retrieving barcode reads from '+ff.split('/')[-1])
    os.system('zcat ' + ff + ' | grep ' + gfp_seq + ' -B 1 >> '+ args.output + '/LARRY_tmp/BC_'+repr(i)+'.txt')


if len(no_abundant_bcs) > 0:
    print('\nThe following libraries had no abundant_barcodes.pickle file')
    for f in no_abundant_bcs: print(f)

if len(no_fastq) > 0:
    print('\nThe following libraries had no fastq file')
    for f in no_fastq: print(f)
print('\n')

# Filtering and combining barcode reads

def is_valid(bc):
    return bc[10:12]=='TG' and bc[16:18]=='CA' and bc[22:24]=='GT' and bc[28:30]=='AG'

out = gzip.open(args.output + '/' + args.name + '_LARRY_sorted_and_filtered_barcodes.fastq.gz','wb')
out_larry_id = open(args.output + '/' + args.name + '_LARRY_cell_barcode_fluorophore_map.txt','wb')

for f in os.listdir(args.output + '/LARRY_tmp'):
    print('Filtering reads from ',f)
    #ab_bc = pickle.load(open(lib+'/abundant_barcodes.pickle','rb'))
    cell_barcodes = pd.read_csv(args.barcode, sep='\t')
    ab_bc = cell_barcodes.index
    
    for xx in open(args.output + '/LARRY_tmp/'+f).read().strip('\n').split('--'):
        xx = xx.strip('\n').split('\n')
        gfp_bc = xx[1].split(gfp_seq)[1]
        cell_bc = xx[0].split('RG:Z:')[1]
        if len(gfp_bc) >= 34:
            gfp_bc = gfp_bc[:34]
            if is_valid(gfp_bc):
                if cell_bc in ab_bc:
                    tag = '>'+cell_bc
                    out.write((tag+'\n').encode('utf-8'))
                    out.write((gfp_bc+'\n').encode('utf-8'))
                    txt = cell_bc + ',' + gfp_bc[0:6]
                    out_larry_id.write((txt+'\n').encode('utf-8'))
                    
out.write('\n'.encode('utf-8'))
out.close()
out_larry_id.close()

