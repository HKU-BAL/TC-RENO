import os
import collections
import sys

gene_reads = set()

path = sys.argv[1]#'/autofs/bal19/qhli/RNA/review/process/no_replicate/cluster/'
f = open(path + '/intersect_reads_gene.bed','r')
rep_cluster = collections.defaultdict(set)
for line in f.readlines():
    line = line.split('\n')[0]
    line = line.split('\t')
    if line[23] != '.':
       name = line[3]
       rep = line[23] + '|' + line[5]#.split(',')[0]
       rep_cluster[rep].add(name)

f.close()


f = open(sys.argv[2],'r')
read_info = {}
for line in f.readlines():
    info = line
    line = line.split('\n')[0]
    line = line.split('\t')
    info = line[0]
    for n in range(1,12):
        info += '\t' + line[n]
    read_info[line[3]] = info + '\n'
f.close()

try:
    os.mkdir(path + '/annotation_round1/')
except OSError:
    pass

for rep in rep_cluster:
#    print rep
    f_write = open(path + '/annotation_round1/'+rep,'w')
    for information in rep_cluster[rep]:
        f_write.write(read_info[information])
    f_write.close()

