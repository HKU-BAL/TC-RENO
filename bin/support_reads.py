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
#    print reads
    
    exp[reads[0]] = float(score)
    for read in reads:
        if len(read) > 5:
           round_rep_subreads[reads[0]].add(read)
f.close()

remove_rep = collections.defaultdict(set)
cutoff = float(sys.argv[5])
files = os.listdir(path + '/cluster_round3_TC_RENO/')

for file in files:
    f = open(path + '/cluster_round3_TC_RENO/' + file, 'r')
    for line in f.readlines():
        line = line.strip().split('\t')
        names = line[12].split('|')[0].split(',')
        for name in names:
            
            if len(name)>3 and exp[name] <cutoff:
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

read_pos = {}
f = open(sys.argv[6],'r')
for line in f.readlines():
    info = line
    line = line.strip().split('\t')
    read_pos[line[3]] = info
f.close()


cutoff1 = float(sys.argv[4])
f_write = open(sys.argv[2],'w')
f_write1 = open(sys.argv[3],'w')



f_write2 = open(path + '/read_percent','w')
for r in read_num:
    f_write2.write(r + '\t' + str(read_num[r]) + '\n')
f_write2.close()


f_write1.write('isoform\tsupport_reads\tread_count\n')
f = open(path + '/rep_round2.support.bed','r')
for line in f.readlines():
  #  info = line.strip()
    line = line.strip().split('|')[0].split('\t')
    if line[0] not in remove:
  
       num = 0
       reads = ''
      
      
       for read in round_rep_subreads[line[0]]:
    
           num += 1.0/read_num[read]
          
           if read != line[0]:
              reads += read +','
       if exp[line[0]] >=cutoff1:
              f_write.write(read_pos[line[0]])
              f_write1.write(str(line[0]) + '\t' + reads+ '\t' + str(num) + '\n')
f_write.close()
f_write1.close()
f.close()




