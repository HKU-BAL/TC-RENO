import sys
import copy
from pre_round1 import read_bigg, bigg_count_write,write_bigg, bigglist_to_bedfile
import collections
from collections import OrderedDict
import os
from pre_cluster_round1 import *
from utils import parmap

#from trackcluster.batch import *


def group_bigg_by_gene(bigglist):
    gene_bigg=OrderedDict()
    for bigg in bigglist:
        try:
            gene_bigg[bigg.geneName].append(bigg)
        except KeyError:
            gene_bigg[bigg.geneName]=[]
            gene_bigg[bigg.geneName].append(bigg)
    return gene_bigg
#sys.argv
subk=[]
files = os.listdir(sys.argv[1])

try:
    os.mkdir(sys.argv[2])
except OSError:
    pass

for file in files:
    subk.append(file)
def cluster_isoforms(key,batchsize=1000, intronweight=float(sys.argv[3]), full=False):
    gff_file = sys.argv[1] + key + "/" + key + "_gff.bed"
    nano_file = sys.argv[1] + key + "/" + key + "_nano.bed"
    biggout = sys.argv[2] + key +".bed"
    if full is False:
        if os.stat(nano_file).st_size == 0: 
            return 0
        if os.path.isfile(biggout):
            return 0
    bigg_gff = read_bigg(gff_file)
    bigg_nano_raw = read_bigg(nano_file)
    try:
        bigg_nano = prefilter_smallexon(bigg_nano_raw, bigg_gff, key, cutoff=50,)
        n_count = 100
        n = 0
        if bigg_nano is None:
            return 0

        try:
            while n < n_count and len(bigg_nano) > batchsize:

                bigg_1 = bigg_nano[:batchsize]
                bigg_2 = bigg_nano[batchsize:]
                _, bigg_list_by1 = flow_cluster(bigg_1, bigg_gff, key, intronweight=intronweight)
                bigg_nano = add_subread_bigg(bigg_list_by1 + bigg_2)
                n += 1
            D, bigg_nano_new = flow_cluster(bigg_nano, bigg_gff,key,  intronweight=intronweight)
            bigg_nano_new = add_subread_bigg(bigg_nano_new)
            for bigg in bigg_nano_new:
                bigg.write_subread()

            bigg_count_write(bigg_nano_new, out=biggout)
        except Exception as e:
            print 'c',e
    except Exception as e:
        print 'd',e
    return 1
#from utils import parmap
parmap(cluster_isoforms, subk,int(sys.argv[4]))

