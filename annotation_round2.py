import sys
import os
import collections

sub = collections.defaultdict(set)
path = sys.argv[1]
f = open(path + 'intersect_novel_rep_round1.bed','r')
for line in f.readlines():
    line = line.strip().split('\t')
    if float(line[12])/(int(line[2]) - int(line[1])) > 0.8:
       sub[line[3]].add(line[9])
f.close()

f = open(path + 'merge_novel.bed','r')
rep_cluster = {}
for line in f.readlines():
    line = line.split('\n')[0]
    line = line.split('\t')
    names = line[4].split(',')
    rep = names[0]
    rep_cluster[rep] = []
    for name in names:
        name = name
        rep_cluster[rep].append(name)
f.close()

f = open(sys.argv[2],'r')
read_info = {}
for line in f.readlines():
    info = line
    line = line.split('\n')[0]
    line = line.split('\t')
    read_info[line[3]] = info
f.close()

try:
    os.mkdir(sys.argv[1] + '/annotation_round2/')
except OSError:
    pass

for rep in rep_cluster:
    f_write = open(sys.argv[1]  + 'annotation_round2/' + rep,'w')
    existing = set()
    for read in rep_cluster[rep]:
        f_write.write(read_info[read])
        for subread in sub[line[3]]:
            if subread not in existing:
               f_write.write(read_info[subread])
               existing.add(subread)

    f_write.close()




