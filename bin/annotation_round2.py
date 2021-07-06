import os
import collections
import sys

rep_cluster = collections.defaultdict(set)
path = sys.argv[1]#'/autofs/nas4/qhli/RNA/NAR_col_cluster/'
f = open(path + '/merge_round1.bed','r')
for line in f.readlines():
    line = line.strip().split('\t')
    rep = line[0] + ':' + line[1] + '-' + line[2] + '|' + line[3]
    names = line[4].split(',')
    if len(names)>=1:
    
       rep_cluster[rep] = []
       for name in names:
           if len(name) > 5:
              rep_cluster[rep].append(name)

f.close()

#path = '/autofs/bal19/qhli/RNA/review/process/validate_depth/NAR_col0/'
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
    os.mkdir(path + '/annotation_round2/')
except OSError:
    pass

for rep in rep_cluster:
    f_write = open(path + '/annotation_round2/'+rep,'w')
    for information in rep_cluster[rep]:
        f_write.write(read_info[information])
    f_write.close()







