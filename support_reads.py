import os
import sys
from collections import Counter
import collections

round_subreads_rep = collections.defaultdict(set)
round_rep_subreads = collections.defaultdict(set)

exp = {}
remove = set()
path = sys.argv[1]
f = open(path + '/rep_round2.support.bed','r')
for line in f.readlines():
    score = line.strip().split('|')[1]
    reads = line.strip().split('|')[0].split('\t')
    line = line.strip().split('\t')
    exp[line[0]] = float(score)
    for read in reads:
        round_rep_subreads[line[0]].add(read)
f.close()

remove_rep = collections.defaultdict(set)
cutoff = float(sys.argv[4])
files = os.listdir(path + '/cluster_round3/')

for file in files:
    f = open(path + '/cluster_round3/' + file, 'r')
    for line in f.readlines():
        line = line.strip().split('\t')
        names = line[12].split('|')[0].split(',')
        for name in names:
            if len(name)>0 and '.' not in name and exp[name] <cutoff:
               remove.add(name)

               remove_rep[name].add(line[3])
    f.close()


for isoform in remove:
    for rep in remove_rep[isoform]:
        for read in round_rep_subreads[isoform]:
            round_rep_subreads[rep].add(read)

reads = []
for rep in round_rep_subreads:
    if rep not in remove:
       for i in round_rep_subreads[rep]:
           reads.append(i)
read_num = Counter(reads)
print len(remove)
f_write = open(sys.argv[2],'w')
f_write1 = open(sys.argv[3],'w')
f_write1.write('isoform\tsupport_reads\tread_count\n')
f = open(path + 'rep_round2.bed','r')
for line in f.readlines():
    info = line.strip()
    line = line.strip().split('\t')
    if line[3] not in remove:
       f_write.write(info+'\n')      
       num = 0
       reads = ''
       exp = 0
       for read in round_rep_subreads[line[3]]:
           if '.' not in read:
              num += 1.0/read_num[read]
           if read != line[3]:
              reads += read +','
     
       f_write1.write(str(line[3]) + '\t' + reads+ '\t' + str(num) + '\n')
f_write.close()
f_write1.close()
f.close()




