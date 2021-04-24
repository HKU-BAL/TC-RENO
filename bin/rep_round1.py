import os
import sys
import collections

read_info = {}
length = {}
f = open(sys.argv[2],'r')
for line in f.readlines():
    info = line.strip()
    line = line.strip().split('\t')
    exon  = line[10].split(',')
    number = 0 
    for e in exon:
        if len(e)>0:
           number += int(e)
    length[line[3]] = number
    read_info[line[3]] = info + '\n'
f.close()


subreads_rep = collections.defaultdict(set)
rep_subreads = collections.defaultdict(set)
path = sys.argv[1]+  '/cluster_round1/'
path1 = sys.argv[1]+ '/annotation_round1/'
files = os.listdir(path)
for file in files:
    f = open(path + file ,'r')
    for line in f.readlines():
        line = line.strip().split('\t')
        rep = line[3]
        exon_len = 0
        names = line[12].split('|')[0].split(',')
        if '.' not in rep:
          for name in names:
            if '.' in name:
                if exon_len == 0:
                   exon_len = length[name]
                   rep = name
                elif length[name] > exon_len:
                   exon_len = length[name]
                   rep = name
        for name in names:
            if len(name)>0 and name != rep:
               rep_subreads[rep].add(name)
               subreads_rep[name].add(rep) 
        rep_subreads[rep].add(line[3]) 
        if rep != line[3]:
           subreads_rep[line[3]].add(rep)               
    f.close()

all_reads = set()
f = open(sys.argv[3],'r')
for line in f.readlines():
    info = line.strip()
    line = line.strip().split('\t')
    all_reads.add(line[3])
    read_info[line[3]] = info + '\n'
f.close()


remove = set()
for rep in rep_subreads:
    if rep in subreads_rep:
       for rep1 in subreads_rep[rep]:
           for read in rep_subreads[rep]:
               rep_subreads[rep1].add(read)
       remove.add(rep)

path2 = sys.argv[1] 
f_write = open(path2 + 'rep_round1.bed','w')
f_write_short  = open(path2 + 'rep_round1.short.bed','w')
f_novel = open(path2 + 'novel.bed','w')
f_write1 = open(path2 + 'rep_round1.support.bed','w')
gene_reads = set()
for rep in rep_subreads:
    if rep not in remove:
       f_write.write(read_info[rep] )
       information = read_info[rep].split('\t')
       f_write_short.write(information[0] + '\t' + information[1] + '\t' + information[2] + '\t' + information[3] + '\t' + information[4] + '\t' + information[5] +'\n')
       info = ''
       for i in rep_subreads[rep]:
           gene_reads.add(i)
           if i != rep:
              info += '\t' + i #+ '\t' 
       f_write1.write(rep+ info +'\n' )
f_write.close()
f_write1.close()
f_write_short.close()

remain = all_reads - gene_reads

for r in remain:
    information = read_info[r].split('\t')
    f_novel.write(information[0] + '\t' + information[1] + '\t' + information[2] + '\t' + information[3] + '\t' + information[4] + '\t' + information[5] +'\n')
f_novel.close()



