from track import bigGenePred
from collections import defaultdict
import scipy
import operator
from pre_round3 import summary,pandas_summary1,wrapper_bedtools_intersect3,wrapper_bedtools_intersect2, bigglist_to_bedfile, pandas_summary, add_subread_bigg, get_readall_bigg
from utils import del_files


def flow_cluster(bigg_nano, intronweight=2.3):
    cutoff= 0.005
    bigg_nano.sort(key=operator.attrgetter("chromStart"))
    bigg_list=add_subread_bigg(bigg_nano)
    D2, bigg_list_by2 =cal_distance(bigg_list, intronweight=intronweight)
    D_remain, bigg_l3=filter_D(D2, bigg_list_by2,cutoff=cutoff)
    return D_remain, bigg_l3

def getij(bigg_list):
    ij_list=[]
    for i in range(len(bigg_list)):
        for j in range(len(bigg_list)):
            ij_list.append((i,j))
    return ij_list


def get_pos_dic(bigg_list):
    pos_dic={}
    for n, bigg in enumerate(bigg_list):
        pos_dic[bigg.name]=n
    return pos_dic


def select_list(bigg_list, keep):
    bigg_list_new=[]
    for i in keep:
        bigg_list_new.append(bigg_list[i])
    return bigg_list_new

def select_D(D, keep):
    D=D[keep,:]
    D=D[:,keep]

    return D


def cal_distance(bigg_list, intronweight=2.3):
    bigg_list.sort(key=operator.attrgetter("chromStart"))

    for i in bigg_list:
        i.get_exon()
        i.to_bedstr()
    length=len(bigg_list)
    D_exon=scipy.ones([length, length])
    D_intron=scipy.ones([length, length])
    D =scipy.ones([length, length])
    pos_dic=get_pos_dic(bigg_list)
    file_exon, file_intron = bigglist_to_bedfile(bigg_list)
    pos_dic=get_pos_dic(bigg_list)
    D_intersect = scipy.zeros([length, length])
    exon_out=wrapper_bedtools_intersect3(file_exon, file_exon)
    exon_i=summary(exon_out, pos_dic,D_intersect )
    D_intersect = scipy.zeros([length, length])
    intron_out=wrapper_bedtools_intersect3(file_intron, file_intron)
    intron_i=summary(intron_out,  pos_dic,D_intersect )
    ij_list=getij(D)
    for i,j in ij_list:
        if exon_i[i,j]!= 0:
        
           D_exon[i, j] =1- float(exon_i[i,j])/bigg_list[i].exonlen
           if intron_i[i,j]!=0:
               D_intron[i, j] = (1- float(intron_i[i,j])/bigg_list[i].intronlen)
           elif  bigg_list[i].intronlen ==0:
               D_intron[i, j] = 0
    D=(D_exon+intronweight*D_intron)/float(1+intronweight)
    del_files([exon_out, intron_out, file_exon, file_intron])
    return D, bigg_list

def filter_D(D, bigg_list, cutoff=0.005):
    fullset=set(range(len(D)))
    drop=set()
    ij_list=getij(D)
    for n in range(5):
        for i,j in ij_list:
            if D[i,j]<cutoff and "." not in bigg_list[i].name:
                if i==j:
                    pass
                elif bigg_list[i].exonlen<bigg_list[j].exonlen:
                    drop.add(i)
                    bigg_list[j].subread.add(bigg_list[i].name)
                    bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
        for i,j in reversed(ij_list):
            if D[i,j]<cutoff and "." not in bigg_list[i].name:
                if i==j:
                    pass
                elif bigg_list[i].exonlen<bigg_list[j].exonlen:
                    drop.add(i)
                    bigg_list[j].subread.add(bigg_list[i].name)
                    bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
    keep=fullset-drop
    for n, bigg in enumerate(bigg_list):
        if bigg.ttype=="isoform_anno":
            keep.add(n)
    keepl = sorted(list(keep))
    bigg_list_new = select_list(bigg_list, keepl)
    D = select_D(D, keepl)

    return D, bigg_list_new





