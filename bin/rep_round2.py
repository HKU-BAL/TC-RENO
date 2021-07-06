import os
import collections
from collections import Counter
import sys

read_info = {}
round_subreads_rep = collections.defaultdict(set)
round_rep_subreads = collections.defaultdict(set)

round_rep_reads =  collections.defaultdict(set)
f = open(sys.argv[1]+'/rep_round1.support.bed','r')
for line in f.readlines():
    line = line.strip().split('\t')
    for n in range(1, len(line)):
        if len(line[n]) > 5:
           round_rep_subreads[line[0]].add(line[n])
           round_subreads_rep[line[n]].add(line[0])
    round_rep_subreads[line[0]].add(line[0])
f.close()
'''
f = open(sys.argv[1]+'/rep_round1.bed','r')
for line in f.readlines():
    info = line.strip()
    line = line.strip().split('\t')
    read_info[line[3]] = info + '\n'
f.close()
'''
path = sys.argv[1]+'/cluster_round2_TC_RENO/'
path1 = sys.argv[1]+ '/annotation_round2/'
path3 = sys.argv[1]+ '/annotation_round3/'

try:
    os.mkdir(path3)
except OSError:
    pass


files = os.listdir(path)
for file in files:
    f = open(path1 + file, 'r')
    for line in f.readlines():
        info = line
        line = line.strip().split('\t')
        read_info[line[3]] = info 
    f.close()

    f = open(path + file ,'r')
    f_write = open(path3 + file,'w')
    for line in f.readlines():
        line = line.strip().split('\t')
        rep = line[3]
        names = line[12].split('|')[0].split(',')
        f_write.write(read_info[line[3]] )
        for name in names:
            if len(name)>3:
               round_rep_subreads[rep].add(name)
               round_subreads_rep[name].add(rep)
        round_rep_subreads[rep].add(rep)
    f.close()
    f_write.close()



remove = set()
for rep in round_rep_subreads:

    if rep in round_subreads_rep:
       for rep1 in round_subreads_rep[rep]:
           for read in round_rep_subreads[rep]:
               round_rep_subreads[rep1].add(read)
       remove.add(rep)


reads = []
for rep in round_rep_subreads:
    if rep not in remove:
       for i in round_rep_subreads[rep]:
           reads.append(i)
read_num = Counter(reads)  

path2 = sys.argv[1]
#f_write = open(path2 + 'rep_round2.bed','w')
f_write1 = open(path2 + '/rep_round2.support.bed','w')
#f_write_short = open(path2 + 'rep_round2.short.bed','w')
for rep in round_rep_subreads:
    if rep not in remove:
       number = 0
       for i in round_rep_subreads[rep]:
            number += 1.0/read_num[i]
#       print rep, number
      # if number >= cutoff:
 #      f_write.write(read_info[rep])
      # information = read_info[rep].split('\t')
    #      f_write_short.write(information[0] + '\t' + information[1] + '\t' + information[2] + '\t' + information[3] + '\t' + information[4] + '\t' + information[5] +'\n')
       info = ''
       for i in round_rep_subreads[rep]:
           if i != rep:
              info += '\t' + i
       f_write1.write(rep + info + "|" + str(number) + '\n' )
#f_write.close()
f_write1.close()
#f_write_short.close()




