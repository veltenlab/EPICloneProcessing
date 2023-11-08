import numpy as np
import networkx as nx
import pickle, json, gzip
import argparse
import csv

parser = argparse.ArgumentParser(description='Transform LARRY barcodes to clonal annotation.')
parser.add_argument('--fastq', help='The fastq with the cell and LARRY barcodes')
parser.add_argument('--barcode', help='The cell barcode file')
parser.add_argument('--output', help='The output directory')
parser.add_argument('--nreads', type=int, help='The minimum number of reads', default=10)
parser.add_argument('--nhamming', help='The hamming distance used for merging', default=3)

args = parser.parse_args()

N_READS = int(args.nreads)
N_HAMMING = int(args.nhamming)
cell_bcs = open(args.barcode)

counts = {}
f = gzip.open(args.fastq)
l = f.readline().decode("utf-8").strip('\n')
current_tag = []
i = 0
print('Reading in all barcodes')
while not (l == '' and len(current_tag)==0):
    i += 1
    if i % (3*10**6)==0: print('Processed '+repr(int(i/3))+' reads')
    if l == '':
        current_tag = []
    elif l[0] == '>':
        current_tag = l[1:].split(',')
    elif l != '' and len(current_tag)==1:
        current_tag.append(l)
        current_tag = tuple(current_tag)
        if not current_tag in counts: counts[current_tag] = 0
        counts[current_tag] += 1
        
    l = f.readline().decode("utf-8").strip('\n')

counts_filtered = {k:v for k,v in counts.items() if v >= N_READS}
print('Retaining '+repr(len(counts_filtered))+ ' out of '+repr(len(counts))+' (Cell-BC,GFP-BC) combinations')

def hamming(bc1,bc2): return np.sum([x1 != x2 for x1,x2 in zip(bc1,bc2)])

all_gfp_bcs = sorted(set([k[1] for k in counts_filtered]))
good_gfp_bcs = []
bc_map = {}
for i,bc1 in enumerate(all_gfp_bcs):
    if i > 0 and i % 500 == 0: print('Mapped '+repr(i)+' out of '+repr(len(all_gfp_bcs))+' barcodes')
    mapped = False
    for bc2 in good_gfp_bcs:
        if hamming(bc1,bc2) <= N_HAMMING:
            mapped = True
            bc_map[bc1] = bc2
            break
    if not mapped:
        good_gfp_bcs.append(bc1)

counts_filtered_copy = counts_filtered.copy()
for k in counts_filtered_copy.keys():
    if k[1] in bc_map.keys():
        updated_barcode = tuple([k[0], bc_map[k[1]]])
        old_item = counts_filtered[k]
        del counts_filtered[k]
        counts_filtered[updated_barcode] = old_item

print('\nCollapsed '+repr(len(bc_map))+' barcodes')
for bc in good_gfp_bcs: bc_map[bc] = bc

BC_set = sorted(set([bc for bc in good_gfp_bcs if bc != '']))
clone_mat = np.zeros((len(good_gfp_bcs),len(BC_set)))
for i,bc in enumerate(good_gfp_bcs):
    if bc != '':
        j = BC_set.index(bc)
        clone_mat[i,j] = 1
clone_mat = np.array(clone_mat,dtype=int)

with open(args.output + '_filtered_clone_combinations.csv', 'w') as f:
    for key in counts_filtered.keys():
        f.write("%s,%s\n"%(key,counts_filtered[key]))
 
np.savetxt(args.output + '_clone_mat.csv',clone_mat,delimiter=',',fmt='%i');
np.save(args.output + '_clone_mat.npy',clone_mat);
open(args.output + '_barcode_list.txt','w').write('\n'.join(BC_set));


