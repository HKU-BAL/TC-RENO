import os
import sys

path = sys.argv[1]
f = open(path + 'merge_rep_round2.bed','r')
rep_cluster = {}
for line in f.readlines():
    line = line.split('\n')[0]
    line = line.split('\t')
    names = line[4].split(',')
    
    if len(names)>0:

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
f = open(sys.argv[3],'r')
for line in f.readlines():
    info = line
    line = line.split('\n')[0]
    line = line.split('\t')
    read_info[line[3]] = info
f.close()

try:
    os.mkdir(sys.argv[1] + '/annotation_round3/')
except OSError:
    pass

for rep in rep_cluster:
    
    f_write = open(path + '/annotation_round3/'+rep,'w')
    for information in rep_cluster[rep]:
        f_write.write(read_info[information])
    f_write.close()


