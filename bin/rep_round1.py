import os
import collections
import sys

subreads_rep = collections.defaultdict(set)
rep_subreads = collections.defaultdict(set)
#path = '/autofs/nas4/qhli/RNA/NAR_col_cluster/'
path = sys.argv[2]
files = os.listdir(path + '/cluster_round1_TC_RENO/')
for file in files:
    f = open(path + '/cluster_round1_TC_RENO/' + file ,'r')
    for line in f.readlines():
        line = line.strip().split('\t')
        rep = line[3]
        names = line[12].split('|')[0].split(',')
        for name in names:      
#            if len(name)> 5:
                rep_subreads[rep].add(name)
                subreads_rep[name].add(rep)
        rep_subreads[rep].add(line[3])
    f.close()

remove = set()
for rep in rep_subreads:
    if rep in subreads_rep:
       for rep1 in subreads_rep[rep]:
           for read in rep_subreads[rep]:
               rep_subreads[rep1].add(read)
       remove.add(rep)

read_info = {}
f = open(sys.argv[1],'r')
for line in f.readlines():
    info = line.strip()
    line = line.strip().split('\t')
    read_info[line[3]] = line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t0\t' + line[5] + '\n'
f.close()

reads = set()
f_write = open(path + '/rep_round1.bed','w')
f_write1 = open(path + '/rep_round1.support.bed','w')
for rep in rep_subreads:
    if rep not in remove:
       f_write.write(read_info[rep])
       info = ''
       for i in rep_subreads[rep]:
           if i != rep:
              info += '\t' + i #+ '\t'
              reads.add(i)
       reads.add(rep)
       f_write1.write(rep+ info +'\n' )
#f_write.close()
#f_write1.close()

for read in read_info:
    if read not in reads:
       f_write.write(read_info[read])
       f_write1.write(read+'\n' )
f_write.close()
f_write1.close()




       


