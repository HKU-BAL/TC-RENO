import sys, argparse
import sys
import copy
#sys.path.insert(0, "/mnt/bal19/qhli/RNA/TC_pipeline/bin/")

from tracklist import read_bigg, write_bigg, bigglist_to_bedfile
import collections
from collections import OrderedDict
import os

def group_bigg_by_gene(bigglist):
    gene_bigg=OrderedDict()
    for bigg in bigglist:
        try:
            gene_bigg[bigg.geneName].append(bigg)
        except KeyError:
            gene_bigg[bigg.geneName]=[]
            gene_bigg[bigg.geneName].append(bigg)
    return gene_bigg
read_n={}

try:
    os.mkdir(sys.argv[4] + 'cluster/')
except OSError:
    pass

#print sys.argv
read_gene={}
f=open(sys.argv[1],'r') #args.intersect)
for line in f.readlines():
#    print line
    line_l=line.split("\t")
    name=line_l[3]+"$"+line_l[1]+"$"+line_l[2]
    gene=line_l[-3]
    if name not in read_gene :
       read_gene[name]  = set()
       read_gene[name].add(gene)
    else:
       read_gene[name].add(gene)
f.close()

nano_bigg=read_bigg(sys.argv[2])#arg.i)
nano_bigg_new=[]
for bigg in nano_bigg:
    try:
        genes=read_gene[bigg.name+'$'+str(bigg.chromStart)+'$'+str(bigg.chromEnd)]
        for gene in genes:
           bigg_copy = copy.deepcopy(bigg)
           bigg_copy.geneName=gene
           nano_bigg_new.append(bigg_copy)
    except KeyError:
        pass
gene_nano=group_bigg_by_gene(nano_bigg_new)
anno_bigg=read_bigg(sys.argv[3])#args.r)
gene_anno=group_bigg_by_gene(anno_bigg)

try:
    os.mkdir(sys.argv[4] + '/annotation_round1/')
except OSError:
    pass
for gene, nano_bigg in gene_nano.iteritems():
    anno_bigg=gene_anno[gene]
    try:    
        os.mkdir(sys.argv[4] + "/annotation_round1/"+gene)
    except OSError:
        pass

    anno_out= sys.argv[4] + "/annotation_round1/{gene}/{gene}_gff.bed".format(gene=gene)
    nano_out= sys.argv[4] + "/annotation_round1/{gene}/{gene}_nano.bed".format(gene=gene)
    write_bigg(anno_bigg, anno_out)
    write_bigg(nano_bigg, nano_out)


