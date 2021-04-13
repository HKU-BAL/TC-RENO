#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/9/2018 10:01 AM
# @Author  : Runsheng     
# @File    : cluster.py
"""
Generating similarity matrix using differ between gene models
make dendrogram for further plotting
The input is a list of
"""
# self import
from track import bigGenePred
from collections import defaultdict
# third part import
import scipy
import operator
from tracklist import summary,pandas_summary1,wrapper_bedtools_intersect3,wrapper_bedtools_intersect2, bigglist_to_bedfile, pandas_summary, add_subread_bigg, get_readall_bigg
from utils import del_files


def flow_cluster(bigg_nano,  by="ratio_all", cutoff="auto", intronweight=2.3):
    bigg_nano.sort(key=operator.attrgetter("chromStart"))

    if by=="ratio_all":
        by1="ratio"
        by2="ratio_short"
    else:
        by1=by
        by2=by

    if cutoff=="auto":
        cutoff1= 0.05
        cutoff2= 0.001
    else: # expect cutoff as a tuple (0.05, 0.01)
        cutoff1, cutoff2= cutoff
#    for i in bigg_nano:
 #       print 'ippy',i.name,i.chromStart
    # hard code first filter of overalpping of 50 bp
    #bigg_l1=prefilter_smallexon(bigg_nano, bigg_gff, cutoff=50) # using default cutoff 0.95
    bigg_list=add_subread_bigg(bigg_nano)
   # print bigg_list
  #  for i in bigg_list:
   #    print 'happy',i.name,i.chromStart
    # can be change filters
    #print 'end'
    D1, bigg_list_by1=cal_distance(bigg_list, intronweight=intronweight, by=by1)
#    for i in bigg_list_by1:
 #      print 'ippy',i.name,i.chrom
#    print 'finish compute'
    D_remain, bigg_l3 =filter_D(D1, bigg_list_by1, by=by1, cutoff=cutoff1)
    #for i in bigg_l2:
     #  print 'ippy',i.name,i.chrom
#    D2, bigg_list_by2=cal_distance(bigg_l2, intronweight=intronweight, by=by2)
 #   D_remain, bigg_l3=filter_D(D2, bigg_list_by2, by=by2, cutoff=cutoff2)
#    for i in bigg_l3:
 #       print 'ippy',i.name,i.chrom

    # add sanity check
    # the bigg_l3 subreads number together with read number+ bigg_l3=bigg_ll

    #missed_2=get_readall_bigg(bigg_list_by1)-get_readall_bigg(bigg_l2)
    #missed_3=get_readall_bigg(bigg_list_by2)-get_readall_bigg(bigg_l3)

    print "flow cluster",len(bigg_list),  len(bigg_l3)#, len(bigg_l3)

    return D_remain, bigg_l3


def getij(bigg_list):
    ij_list=[]
    for i in range(len(bigg_list)):
        for j in range(1, len(bigg_list)):
            ij_list.append((i,j))
    return ij_list
def getij1(bigg_list):
    ij_list=[]
    for i in range(len(bigg_list)):
        for j in range(i+1, len(bigg_list)):
            ij_list.append((i,j))
    return ij_list

def get_pos_dic(bigg_list):
    pos_dic={}
    for n, bigg in enumerate(bigg_list):
        pos_dic[bigg.name]=n
    return pos_dic


def select_list(bigg_list, keep):
    # re_order D and bigg_list
    bigg_list_new=[]
    for i in keep:
        bigg_list_new.append(bigg_list[i])
    return bigg_list_new

def select_D(D, keep):
    D=D[keep,:]
    D=D[:,keep]

    return D

def prefilter_smallexon(bigg_list,bigg_list_gff, cutoff=50):
    """
    remove two kind of reads:

    1. not same strand as in current annotation
    2. with less than cutoff intersection with current annotation

    :param bigg_list:
    :param cutoff: at least 50bp intersection with current annotation
    :return: retained list
    """
    if len(bigg_list_gff)==0:
        return bigg_list
    #print 'begin',len(bigg_list)
    # filter 1
    #print 'begin',len(bigg_list) 
    strand=bigg_list_gff[0].strand
  #  print strand
  #  for x in bigg_list :
        # 'x',x.strand 
    bigg_list_strand=[x for x in bigg_list if 1==1]
#    print bigg_list_strand
#    for x in bigg_list_strand :
 #       print  'x',x#,x.strand 
    if len(bigg_list_strand)==0:
        return None
  #  for i in bigg_list:
   #     print i,i.chrom,i.name
    # filter 2
   # print 'first',len(bigg_list_strand)
    nano_exon, nano_intron=bigglist_to_bedfile(bigg_list_strand)
    gff_exon, gff_intron=bigglist_to_bedfile(bigg_list_gff)
#    print 'aaaa'
    exonfile=wrapper_bedtools_intersect2(nano_exon, gff_exon)
#    print exonfile
    out_d=pandas_summary1(exonfile)
  #  print 'ccc'
    keep_name=set()
    for k, intersection in out_d.items():
        nano_name, gff_name=k
#        print k, intersection  
       # print nano_name,intersection 
        if intersection > cutoff:
            keep_name.add(nano_name)
   # print len(keep_name)
    bigg_list_new=[]
    for bigg in bigg_list:
        if bigg.name in keep_name:
            bigg_list_new.append(bigg)
   # print 'filter',len(bigg_list_new)
    ### clean up
    del_files([exonfile, nano_intron, gff_intron])
  #  for i in bigg_list_new:
     #  print 'filter',i,i.chrom,i.name
    return bigg_list_new


def cal_distance(bigg_list, intronweight=2.3, by="ratio"):
    """
    :param bigg_list:
    :param intronweight: if 0, do not cal the intron to save time
    :param by: used to cal the distance between two bigg object, can be "ratio", "ratio_short", "length", "length_short"
    :return: D: distance matrix
    """
    #wkdir=set_tmp()
    #os.chdir(wkdir)

    bigg_list.sort(key=operator.attrgetter("chromStart"))

    for i in bigg_list:
        i.get_exon()
        i.to_bedstr()
#        print 'info',i.name,i.chrom,
    length=len(bigg_list)
    D_exon=scipy.ones([length, length])
    D_intron=scipy.ones([length, length])
    D =scipy.ones([length, length])
#    for i in bigg_list:
 #       print i.name,i.chrom
    # get an pos combination and the name of bigg for each i
    # ij_list=getij(bigg_list)
    pos_dic=get_pos_dic(bigg_list)

    # flow begin
    file_exon, file_intron = bigglist_to_bedfile(bigg_list)
 #   f = open(file_exon,'r')
 #   print file_exon
  #  for line in f.readlines():
   #      print line
   # f.close()
    pos_dic=get_pos_dic(bigg_list)
    D_intersect = scipy.zeros([length, length])
    exon_out=wrapper_bedtools_intersect3(file_exon, file_exon)
 #   print 'begin read'
    exon_i=summary(exon_out, pos_dic,D_intersect )
  #  print exon_out
  #  print 'exon_out',exon_out
    D_intersect = scipy.zeros([length, length])
    intron_out=wrapper_bedtools_intersect3(file_intron, file_intron)
 #   print D_intersect

    intron_i=summary(intron_out,  pos_dic,D_intersect )
   # print intron_out
  #  print len(intron_i)
  #  print 'intron_out',intron_i
    ij_list=getij1(D)
    for i,j in ij_list:
        if len(bigg_list[j].blockSizes)==len(bigg_list[i].blockSizes) and exon_i[i,j]!= 0:
           D_exon[i, j] =1- float(exon_i[i,j])/ (bigg_list[i].exonlen + bigg_list[j].exonlen - exon_i[i,j])

           #if bigg_list[i].intronlen ==0 and bigg_list[j].intronlen ==0 :
            #  D_intron[i, j] = 0
           if intron_i[i,j]!=0:
                D_intron[i, j] = (1- float(intron_i[i,j])/ (bigg_list[i].intronlen + bigg_list[j].intronlen - intron_i[i,j]))
           elif bigg_list[i].intronlen ==0 and bigg_list[j].intronlen ==0 :
                D_intron[i, j] = 0
#           print D_exon[i, j], D_intron[i, j], bigg_list[i].exonlen,bigg_list[j].exonlen, exon_i[i,j], bigg_list[i].name,bigg_list[j].name
    D=(D_exon+intronweight*D_intron)/float(1+intronweight)
    '''
    for k, intersection in exon_i.items():
        name1, name2=k
        i=pos_dic[name1]
        j=pos_dic[name2]

        min_length = min(bigg_list[i].exonlen, bigg_list[j].exonlen)
        union = bigg_list[i].exonlen + bigg_list[j].exonlen - intersection
       # print 'exon',k,intersection,min_length,union
        # debug insanity
        if union <=0:
            print "exon", name1, name2, bigg_list[i].exonlen,  bigg_list[j].exonlen, union, intersection
        # debug over

        if by == "ratio":
            # exon could be 0?
            if min_length == 0:
                D_exon[i, j] = 1
                if bigg_list[i].intronlen ==0 and bigg_list[j].intronlen ==0 :
                    D_intron[i, j] = 0
            else:
                similar = float(intersection) / union
                D_exon[i, j] = 1 - similar
                if bigg_list[i].intronlen ==0 and bigg_list[j].intronlen ==0 :
                    D_intron[i, j] = 0
#            print "ratio", name1,name2, D_intron[i, j],D_exon[i, j]
        elif by == "ratio_short":
            # intron could be 0
            if min_length == 0:
                D_exon[i, j] = 1
                if bigg_list[i].intronlen ==0 and bigg_list[j].intronlen ==0 :
                   D_intron[i, j] = 0
            else:
                D_exon[i, j] = 1 - float(intersection) / min_length
                if bigg_list[i].intronlen ==0 and  bigg_list[j].intronlen ==0 :
                   D_intron[i, j] = 0
 #           print "ratio_short", name1,name2, D_intron[i, j],D_exon[i, j]
    
    for k, intersection in intron_i.items():
      #  print k,intersection
        name1, name2 = k
        i = pos_dic[name1]
        j = pos_dic[name2]
       # print 'i j',i,j
        min_length = min(bigg_list[i].intronlen, bigg_list[j].intronlen)
        #print min_length  
        union = bigg_list[i].intronlen + bigg_list[j].intronlen - intersection
       # print 'intron',k,intersection,min_length,union
        #### debug
        if union <=0:
            print "intron",name1, name2, bigg_list[i].intronlen,  bigg_list[j].intronlen, union, intersection
        #### debug over

      #  if by == "ratio":
            # intron could be 0
        if min_length == 0:
                D_intron[i, j] = intronweight
        else:
                #print union
                similar = float(intersection) / union
                D_intron[i, j] =intronweight*( 1 - similar)

#        elif by == "ratio_short":
 #           # intron could be 0
  #          if min_length == 0:
   #             D_intron[i, j] = intronweight
    #        else:
     #           D_intron[i, j] = intronweight*(1 - float(intersection) / min_length)
    for k, intersection in exon_i.items():
        name1, name2=k
        i=pos_dic[name1]
        j=pos_dic[name2]

        min_length = min(bigg_list[i].exonlen, bigg_list[j].exonlen)
        union = bigg_list[i].exonlen + bigg_list[j].exonlen - intersection
       # print 'exon',k,intersection,min_length,union
        # debug insanity
        if union <=0:
            print "exon", name1, name2, bigg_list[i].exonlen,  bigg_list[j].exonlen, union, intersection
        # debug over

        #if by == "ratio":
            # exon could be 0?
        if min_length == 0:
            D_exon[i, j] = 1
            if bigg_list[i].intronlen ==0 and bigg_list[j].intronlen ==0 :
                    D_intron[i, j] = 0
                    D[i,j] = D_exon[i, j]
            else:
                  D[i,j] = (D_exon[i, j]+D_intron[i, j])/float(1+intronweight)
        else:
             similar = float(intersection) / union
             D_exon[i, j] = 1 - similar
             if bigg_list[i].intronlen ==0 and bigg_list[j].intronlen ==0 :
                 D_intron[i, j] = 0
                 D[i,j] = D_exon[i, j]
             else:
                 D[i,j] = (D_exon[i, j]+D_intron[i, j])/float(1+intronweight)
    #    print intersection,union,bigg_list[i].exonlen , bigg_list[j].exonlen , D_exon[i, j],D_intron[i, j], bigg_list[i].name,bigg_list[j].name
#            print "ratio", name1,name2, D_intron[i, j],D_exon[i, j]
       
        elif by == "ratio_short":
            # intron could be 0
            if min_length == 0:
                D_exon[i, j] = 1
                if bigg_list[i].intronlen ==0 and bigg_list[j].intronlen ==0 :
                   D_intron[i, j] = 0
            else:
                D_exon[i, j] = 1 - float(intersection) / min_length
                if bigg_list[i].intronlen ==0 and  bigg_list[j].intronlen ==0 :
                   D_intron[i, j] = 0
        '''
   # D=(D_exon+intronweight*D_intron)/float(1+intronweight)



    # cleanup
    del_files([exon_out, intron_out, file_exon, file_intron])

    # debug:
    #print("D_exon",D_exon)
    #print("D_intron", D_intron)
    #print("D",D)

    #cleanup(remove_all=True)

    return D, bigg_list


def write_D(D, bigg_list_new, outfile="./test/d.csv"):
    bigg_name=[x.name for x in bigg_list_new]

    if outfile is None:
        pass
    else:
        with open(outfile, "w") as fw:
            fw.write(",".join(bigg_name))
            fw.write("\n")
            for i in D:
                str_l=[str(x) for x in i]
                fw.write(",".join(str_l))
                fw.write("\n")


def filter_D(D, bigg_list, by="ratio", cutoff="auto", add_miss=False):

    """
    cutoff selection:
    learn from unc52, <0.025 in ratio_short
    return: index of the matrix that can be retained
    """
    if cutoff=="auto":
        if by=="ratio":
            cutoff= 0.05
        elif by=="ratio_short":
            cutoff= 0.001 # may need to add to 0.01
        elif by=="length" or by=="length_short":
            cutoff=100

    else: # expect two numbers for the cutoff
        cutoff=cutoff
#    print cutoff
    # hard code a cutoff for sw score of SL
    sw_score=11

    fullset=set(range(len(D)))
    drop=set()

    # first
    #print "first"
    #print len(get_readall_bigg(bigg_list))

    #for bigg in bigg_list:
    #    bigg.get_exon()


    # same list
    ij_list=getij(D)

    # add a filter function for D, if the distance between exon is 0, merge the small one with
    # unless need to parer the intronD and exonD separately, or else the filter should be outer function
    for i,j in ij_list:
        if D[i,j]<cutoff  and len(bigg_list[j].blockSizes)==len(bigg_list[i].blockSizes):
            if i==j:
                pass
            else:
                if by=="ratio":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:
                        drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)

                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
                            drop.add(i)
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                            drop.add(j)
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        drop.add(j)
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                        # to add a subread add here

    for i,j in ij_list:
        if D[i,j]<cutoff  and len(bigg_list[j].blockSizes)==len(bigg_list[i].blockSizes):
            if i==j:
                pass
            else:
                if bigg_list[i].exonlen<bigg_list[j].exonlen:
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                elif bigg_list[i].exonlen==bigg_list[j].exonlen:          
                        if i<j:
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

    for i,j in reversed(ij_list):
        if D[i,j]<cutoff  and len(bigg_list[j].blockSizes)==len(bigg_list[i].blockSizes):
            if i==j:
                pass
            else:
                if bigg_list[i].exonlen<bigg_list[j].exonlen:
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j:
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

    for i,j in ij_list:
        if D[i,j]<cutoff and len(bigg_list[j].blockSizes)==len(bigg_list[i].blockSizes):
            if i==j:
                pass
            else:
                if bigg_list[i].exonlen<bigg_list[j].exonlen:
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j:
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
    for i,j in reversed(ij_list):
        if D[i,j]<cutoff and len(bigg_list[j].blockSizes)==len(bigg_list[i].blockSizes):
            if i==j:
                pass
            else:
                if bigg_list[i].exonlen<bigg_list[j].exonlen:
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j:
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

    for i,j in ij_list:
        if D[i,j]<cutoff and len(bigg_list[j].blockSizes)==len(bigg_list[i].blockSizes):
            if i==j:
                pass
            else:
                if bigg_list[i].exonlen<bigg_list[j].exonlen:
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j:
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
    for i,j in reversed(ij_list):
        if D[i,j]<cutoff and len(bigg_list[j].blockSizes)==len(bigg_list[i].blockSizes):
            if i==j:
                pass
            else:
                if bigg_list[i].exonlen<bigg_list[j].exonlen:
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j:
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
    for i,j in ij_list:
        if D[i,j]<cutoff and len(bigg_list[j].blockSizes)==len(bigg_list[i].blockSizes):
            if i==j:
                pass
            else:
                if bigg_list[i].exonlen<bigg_list[j].exonlen:
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j:
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
    for i,j in reversed(ij_list):
        if D[i,j]<cutoff and len(bigg_list[j].blockSizes)==len(bigg_list[i].blockSizes):
            if i==j:
                pass
            else:
                if bigg_list[i].exonlen<bigg_list[j].exonlen:
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j:
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
    for i,j in ij_list:
        if D[i,j]<cutoff and len(bigg_list[j].blockSizes)==len(bigg_list[i].blockSizes):
            if i==j:
                pass
            else:
                if bigg_list[i].exonlen<bigg_list[j].exonlen:
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j:
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
    for i,j in reversed(ij_list):
        if D[i,j]<cutoff and len(bigg_list[j].blockSizes)==len(bigg_list[i].blockSizes):
            if i==j:
                pass
            else:
                if bigg_list[i].exonlen<bigg_list[j].exonlen:
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j:
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
    for i,j in ij_list:
        if D[i,j]<cutoff and len(bigg_list[j].blockSizes)==len(bigg_list[i].blockSizes):
            if i==j:
                pass
            else:
                if bigg_list[i].exonlen<bigg_list[j].exonlen:
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j:
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
    keep=fullset-drop
    # change the default score of gene, no need to add
    for n, bigg in enumerate(bigg_list):
        if bigg.ttype=="isoform_anno":
            keep.add(n)


    # re_order D and bigg_list
    keepl = sorted(list(keep))
    bigg_list_new = select_list(bigg_list, keepl)


    #----------------------------------------#
    #### sanity check for missed ones
    ## collect the missed ones
    pos_dic=get_pos_dic(bigg_list)
    missed_name=get_readall_bigg(bigg_list)-get_readall_bigg(bigg_list_new)

    #print missed_name
    #print "inside"
    #print len(get_readall_bigg(bigg_list_new))
    if add_miss:
        if len(missed_name)>0:
            #print "{} missing bigg found, added back but may affect the isoforms".format(len(missed_name))
            missed_num=set()
            for k in missed_name:
                missed_num.add(pos_dic[k])

            keep=keep.union(missed_num)

            keepl=sorted(list(keep))
            bigg_list_new=select_list(bigg_list, keepl)

            #### end of sanity check


    D = select_D(D, keepl)

    return D, bigg_list_new

def filter_D_no(D, bigg_list, by="ratio", cutoff="auto", add_miss=False):

    """
    cutoff selection:
    learn from unc52, <0.025 in ratio_short
    return: index of the matrix that  in bigg_l2:
     #  print 'ippy',i.name,i.chromi an be retained
    """
    if cutoff=="auto":
        if by=="ratio":
            cutoff=0.99#0.1
        elif by=="ratio_short":
            cutoff=0.99#0.001 # may need to add to 0.01
        elif by=="length" or by=="length_short":
            cutoff=100

    else: # expect two numbers for the cutoff
        cutoff=cutoff

    # hard code a cutoff for sw score of SL
   # sw_score=0

    fullset=set(range(len(D)))
    rep=set()
    drop = set()
    rep_label= set()
    # first
    #print "first"
    #print len(get_readall_bigg(bigg_list))

    #for bigg in bigg_list:
    #    bigg.get_exon()


    # same list
    ij_list=getij1(D)
    sub_reads_info = {}
    cluster_reads =set()
    # add a filter function for D, if the distance between exon is 0, merge the small one with
    # unless need to parer the intronD and exonD separately, or else the filter should be outer function
    dict_read_num = {}
    
    for i,j in ij_list: 
 #       print i,j, bigg_list[i].name,bigg_list[j].name, D[i,j]
#        print bigg_list[i].name,bigg_list[j].name, D[i,j]
        dict_read_num[bigg_list[i].name] = i
#        print 'ength',len(bigg_list[j].blockSizes)
 #       print 'ccc',i,bigg_list[i].name,dict_read_num[bigg_list[i].name]
     #   print bigg_list[i].name,bigg_list[j].name,D[i,j]
        if D[i,j]<cutoff and len(bigg_list[j].blockSizes)==len(bigg_list[i].blockSizes):
#            if i!=j:
 #              cluster_reads.add(i)
  #             cluster_reads.add(j)
    #        if i==j:
     #           pass
           # if i!=j:
             
                cluster_reads.add(i)
                cluster_reads.add(j)
#                print 'intersect',bigg_list[i].name,bigg_list[j].name
                if by=="ratio":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:
                        rep_label.add(j)
                        if j not in sub_reads_info:
                           sub_reads_info[j] = set()
                           sub_reads_info[j].add(i)
                        else:
                           sub_reads_info[j].add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
      #                  print 'intersect',bigg_list[i].name,bigg_list[j].name

#                        if bigg_list[i].subread in bigg_list[j].subread:
 #                          drop.add(i)
#                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
#                        print '1',i,j,bigg_list[i].name,bigg_list[j].name, bigg_list[j].subread
                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
                            rep_label.add(j)
                            if j not in sub_reads_info:
                               sub_reads_info[j] = set()
                               sub_reads_info[j].add(i)
                            else:
                               sub_reads_info[j].add(i)
                          #  drop.add(i)
                          #  print '2',bigg_list[i].name,bigg_list[j].name
       #                     print 'intersect',bigg_list[i].name,bigg_list[j].name

                            bigg_list[j].subread.add(bigg_list[i].name)
        #                    if bigg_list[i].subread in bigg_list[j].subread:
         #                        drop.add(i)
 #                           bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
 #                           print '2',bigg_list[i].name,bigg_list[j].name,bigg_list[j]
                    
                        else:
                            rep_label.add(i)
                            if i not in sub_reads_info:
                               sub_reads_info[i] = set()
                               sub_reads_info[i].add(j)
                            else:
                               sub_reads_info[i].add(j)
                           # drop.add(j)
                            bigg_list[i].subread.add(bigg_list[j].name)
      #                      if bigg_list[j].subread in bigg_list[i].subread:
       #                          drop.add(j)
        #                    print 'intersect',bigg_list[j].name,bigg_list[i].name 
            #               print 3,bigg_list[j].name,bigg_list[i].name
#                            bigg_list[i].subread.add(bigg_list[j].name)
  #                          bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        rep_label.add(i)
                       # rep.add(j)
                        if i not in sub_reads_info:
                            sub_reads_info[i] = set()
                            sub_reads_info[i].add(j)
                        else:
                            sub_reads_info[i].add(j)
   #                     print 4,bigg_list[j].name,bigg_list[i].name
                        bigg_list[i].subread.add(bigg_list[j].name)
     #                   if bigg_list[j].subread in bigg_list[i].subread:
    #                       drop.add(j)i
         #               print 'intersect',bigg_list[j].name,bigg_list[i].name
   #                     bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
#                        print 4,j,i,bigg_list[j].name,bigg_list[i].name,bigg_list[i].subread
  #                      # to add a subread add here
                '''  
                if by=="ratio_short":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:# and bigg_list[i].score<sw_score:
                        drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
       #                 print '1a',i,j,bigg_list[i].name,bigg_list[j].name, bigg_list[j].subread

                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
                            drop.add(i)
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
        #                    print '2a',i,j,bigg_list[i].name,bigg_list[j].name, bigg_list[j].subread

                        else:
                            drop.add(j)
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
         #                   print '3a',j,i,bigg_list[j].name,bigg_list[i].name,bigg_list[i].subread

                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:# and bigg_list[j].score<sw_score:
                        drop.add(j)
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
          #              print '4a',j,i,bigg_list[j].name,bigg_list[i].name,bigg_list[i].subread
                '''  
    for i in rep_label:
      #  sub = []
#        if bigg_list[i].name =="AT1G62830.1":
 #          print "AT1G62830.1",bigg_list[i].subread       
        for j in sub_reads_info[i]:
#            print 'info',i,j,bigg_list[i].name,bigg_list[j].name
            if bigg_list[i].subread> bigg_list[j].subread:

 #               print 'drop_1111111',bigg_list[i].name,bigg_list[j].name
   #             print bigg_list[i].subread
  #              print bigg_list[j].subread
                drop.add(j)
         #   elif bigg_list[i].subread == bigg_list[j].subread:
      #          print 'drop',bigg_list[i].name, bigg_list[j].name     
       #         print bigg_list[i].subread
        #        print bigg_list[j].subread                          
        
    '''
    for i,j in ij_list:
  #      print bigg_list[i].name,bigg_list[j].name, D[i,j]
        dict_read_num[bigg_list[i].name] = i
 #       print 'ccc',i,bigg_list[i].name,dict_read_num[bigg_list[i].name]
        if D[i,j]<cutoff and len(bigg_list[j].blockSizes)==len(bigg_list[i].blockSizes):
    #        for info in bigg_list[i].subread:
#                print 'info',info
            if i==j:
                pass
            else:
                if by=="ratio":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:
  #                      print 'a',bigg_list[j].name,bigg_list[i].name
  #                      rep.add(j)
 #                       bigg_list[j].subread.add(bigg_list[i].name)
                        sub = []
                        rep = []
                        for infos in bigg_list[i].subread:
 #                           print infos
                            if len(infos)>1:
                               sub.append(infos)
                        for infos in bigg_list[j].subread:
                            if len(infos)>1:
                               rep.append(infos)
 #                       print sub
#                        print rep
                        name = 0 
#                        print 'sub',bigg_list[i].subread
 #                       print 'rep',bigg_list[j].subread
                        for infos in sub:
                             if infos not in rep:
                                 name =1
                        if name ==0 and len(sub)>0:                            # else:#if sub in rep or len(sub)==0:

   #                             print '3a',j,i,bigg_list[j].name,bigg_list[i].name#,bigg_list[i].subread

                                drop.add(i)
                        if len(sub)==0:
                                drop.add(i)
#                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
#                        print '1',i,j,bigg_list[i].name,bigg_list[j].name, bigg_list[j].subread
                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
                          #  rep.add(j)
                          #  drop.add(i)
                            sub = []
      #                      print 'a',bigg_list[j].name,bigg_list[i].name
     #                       print sub
    #                        print rep
                            rep = []
                            for infos in bigg_list[i].subread:
                                if len(infos)>1:
                                  sub.append(infos)
                            for infos in bigg_list[j].subread:
                                if len(infos)>1:
                                   rep.append(infos)
                            name = 0
                            for infos in sub:
                                 if infos not in rep:
                                    name = 1
                            if name ==0 and len(sub)>0:#if sub in rep or len(sub)==0:
       #                             print '3a',j,i,bigg_list[j].name,bigg_list[i].name#,bigg_list[i].subread

                                    drop.add(i)
                            if len(sub)==0:
                                    drop.add(i)
                            #for infos in bigg_list[j].subread:
                             #   if len(infos)>1:   
        #                      #     rep.append(infos)
         #                   print sub,rep
                          #  print '2',bigg_list[i].name,bigg_list[j].name
                            #bigg_list[j].subread.add(bigg_list[i].name)
                            #if sub in rep or len(sub)==0:
                             #    drop.add(i)
                              #   print '3a',j,i,bigg_list[j].name,bigg_list[i].name#,bigg_list[i].subread
 #                           bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
 #                           print '2',bigg_list[i].name,bigg_list[j].name,bigg_list[j]
                        else:
                          #  rep.add(i)
                           # drop.add(j)
                            sub = []
          #                  print 'a',bigg_list[i].name,bigg_list[j].name
                            rep = []
                            for infos in bigg_list[j].subread:
                               if len(infos)>1:   
                                  sub.append(infos)
                            for infos in bigg_list[i].subread:
                                if len(infos)>1:   
                                   rep.append(infos)
                          #  print '2',bigg_list[i].name,bigg_list[j].name
                            #bigg_list[j].subread.add(bigg_list[i].name)
           #                 print sub
            #                print rep
                            name =0 
                            for infos in sub:
                                 if infos not in rep:
                                    name =1
                            if name ==0 and len(sub)>0:#if sub in rep or len(sub)==0:
             #                       print '3a',j,i,bigg_list[j].name,bigg_list[i].name#,bigg_list[i].subread

                                    drop.add(j)
                            if len(sub)==0:
                                    drop.add(j)
                        #    for infos in bigg_list[j].subread:
                         #       if len(infos)>1:
                          #         rep.append(infos)
                           # print sub,rep
                           # if sub in rep or len(sub)==0:
                         #   bigg_list[i].subread.add(bigg_list[j].name)
                           # if bigg_list[j].subread in bigg_list[i].subread:
                            #     drop.add(j)
                             #    print '3a',i,j,bigg_list[i].name,bigg_list[j].name#,bigg_list[j].subread
             #               print 3,bigg_list[j].name,bigg_list[i].name
                           # bigg_list[i].subread.add(bigg_list[j].name)
  #                          bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        #rep.add(i)
                       # print 'sub',bigg_list[j].subread
                       # print 'rep',bigg_list[i].subread
                        sub = []
                        rep = []
              #          print 'a',bigg_list[i].name,bigg_list[j].name

                        for infos in bigg_list[j].subread:
                            if len(infos)>1:   
                               sub.append(infos)

                        for infos in bigg_list[i].subread:
                            if len(infos)>1:   
                                rep.append(infos)
               #         print sub,rep
                     #   print 'sub',bigg_list[i].subread
                      #  print 'rep',bigg_list[j].subread
                     #   if sub in rep or len(sub)==0:                       ## rep.add(j)
                           #print '3a',i,j,bigg_list[i].name,bigg_list[j].name
                        #print 4,bigg_list[j].name,bigg_list[i].name
                       # bigg_list[i].subread.add(bigg_list[j].name)
                  #      if bigg_list[j].subread in bigg_list[i].subread:
                       #    print '3a',i,j,bigg_list[i].name,bigg_list[j].name
                      #     drop.add(j)
                        name = 0
                        for infos in sub:
                             if infos not in rep:
                                name = 1
                        if name ==0 and len(sub)>0:#if sub in rep or len(sub)==0:
                #                print '3a',j,i,bigg_list[j].name,bigg_list[i].name#,bigg_list[i].subread

                                drop.add(j)
                        if len(sub)==0:
                                drop.add(j)
    
 #   print 'iandj',type(ij_list),ij_list
  #  print 'reverse',ij_list    
    for i,j in reversed(ij_list):
      #  dict_read_num[bigg_list[i].name] = i
     #   print 'i j',i,j
        if D[i,j]<cutoff:
            if i==j:
                pass
            else:
                if by=="ratio":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:
           #             drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
      #                  print '1',i,j,bigg_list[i].name,bigg_list[j].name, bigg_list[j].subread
                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
            #                drop.add(i)
                          #  print '2',bigg_list[i].name,bigg_list[j].name
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
       #                     print '2',bigg_list[i].name,bigg_list[j].name,bigg_list[j]
                        else:
             #               drop.add(j)
        #                    print 3,bigg_list[j].name,bigg_list[i].name
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:
              #          drop.add(j)
                     #   print 4,bigg_list[j].name,bigg_list[i].name
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
         #               print 4,j,i,bigg_list[j].name,bigg_list[i].name,bigg_list[i].subread
                if by=="ratio_short":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:# and bigg_list[i].score<sw_score:
               #         drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)

                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
                #            drop.add(i)
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                 #           drop.add(j)
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:# and bigg_list[j].score<sw_score:
                  #      drop.add(j)
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)     
    for i,j in (ij_list):  
        if D[i,j]<cutoff:
            if i==j:
                pass
            else:
                if by=="ratio":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:
   #                     drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
      #                  print '1',i,j,bigg_list[i].name,bigg_list[j].name, bigg_list[j].subread
                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
    #                        drop.add(i)
                          #  print '2',bigg_list[i].name,bigg_list[j].name
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
       #                     print '2',bigg_list[i].name,bigg_list[j].name,bigg_list[j]
                        else:
     #                       drop.add(j)
        #                    print 3,bigg_list[j].name,bigg_list[i].name
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:
      #                  drop.add(j)
                     #   print 4,bigg_list[j].name,bigg_list[i].name
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
         #               print 4,j,i,bigg_list[j].name,bigg_list[i].name,bigg_list[i].subread
                if by=="ratio_short":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:# and bigg_list[i].score<sw_score:
       #                 drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)

                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
        #                    drop.add(i)
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
         #                   drop.add(j)
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:# and bigg_list[j].score<sw_score:
          #              drop.add(j)
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread) 
   
    for i,j in reversed(ij_list):
      #  dict_read_num[bigg_list[i].name] = i
     #   print 'i j',i,j
        if D[i,j]<cutoff:
            if i==j:
                pass
            else:
                if by=="ratio":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:
           #             drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
      #                  print '1',i,j,bigg_list[i].name,bigg_list[j].name, bigg_list[j].subread
                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
            #                drop.add(i)
                          #  print '2',bigg_list[i].name,bigg_list[j].name
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
       #                     print '2',bigg_list[i].name,bigg_list[j].name,bigg_list[j]
                        else:
             #               drop.add(j)
        #                    print 3,bigg_list[j].name,bigg_list[i].name
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:
              #          drop.add(j)
                     #   print 4,bigg_list[j].name,bigg_list[i].name
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
         #               print 4,j,i,bigg_list[j].name,bigg_list[i].name,bigg_list[i].subread
                if by=="ratio_short":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:# and bigg_list[i].score<sw_score:
               #         drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)

                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
                #            drop.add(i)
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                 #           drop.add(j)
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:# and bigg_list[j].score<sw_score:
                  #      drop.add(j)
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
  
    
    for i,j in (ij_list):
        if D[i,j]<cutoff:
            if i==j:
                pass
            else:
                if by=="ratio":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:
   #                     drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
      #                  print '1',i,j,bigg_list[i].name,bigg_list[j].name, bigg_list[j].subread
                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
    #                        drop.add(i)
                          #  print '2',bigg_list[i].name,bigg_list[j].name
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
       #                     print '2',bigg_list[i].name,bigg_list[j].name,bigg_list[j]
                        else:
     #                       drop.add(j)
        #                    print 3,bigg_list[j].name,bigg_list[i].name
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:
      #                  drop.add(j)
                     #   print 4,bigg_list[j].name,bigg_list[i].name
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
         #               print 4,j,i,bigg_list[j].name,bigg_list[i].name,bigg_list[i].subread
                if by=="ratio_short":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:# and bigg_list[i].score<sw_score:
       #                 drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)

                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
        #                    drop.add(i)
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
         #                   drop.add(j)
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:# and bigg_list[j].score<sw_score:
          #              drop.add(j)
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)


    for i,j in reversed(ij_list):
      #  dict_read_num[bigg_list[i].name] = i
     #   print 'i j',i,j
        if D[i,j]<cutoff:
            if i==j:
                pass
            else:
                if by=="ratio":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:
           #             drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
      #                  print '1',i,j,bigg_list[i].name,bigg_list[j].name, bigg_list[j].subread
                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
            #                drop.add(i)
                          #  print '2',bigg_list[i].name,bigg_list[j].name
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
       #                     print '2',bigg_list[i].name,bigg_list[j].name,bigg_list[j]
                        else:
             #               drop.add(j)
        #                    print 3,bigg_list[j].name,bigg_list[i].name
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:
              #          drop.add(j)
                     #   print 4,bigg_list[j].name,bigg_list[i].name
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
         #               print 4,j,i,bigg_list[j].name,bigg_list[i].name,bigg_list[i].subread
                if by=="ratio_short":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:# and bigg_list[i].score<sw_score:
               #         drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)

                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
                #            drop.add(i)
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                 #           drop.add(j)
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:# and bigg_list[j].score<sw_score:
                  #      drop.add(j)
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

    for i,j in (ij_list):
      #  dict_read_num[bigg_list[i].name] = i
     #   print 'i j',i,j
        if D[i,j]<cutoff:
            if i==j:
                pass
            else:
                if by=="ratio":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:
           #             drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
      #                  print '1',i,j,bigg_list[i].name,bigg_list[j].name, bigg_list[j].subread
                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
            #                drop.add(i)
                          #  print '2',bigg_list[i].name,bigg_list[j].name
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
       #                     print '2',bigg_list[i].name,bigg_list[j].name,bigg_list[j]
                        else:
             #               drop.add(j)
        #                    print 3,bigg_list[j].name,bigg_list[i].name
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:
              #          drop.add(j)
                     #   print 4,bigg_list[j].name,bigg_list[i].name
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
         #               print 4,j,i,bigg_list[j].name,bigg_list[i].name,bigg_list[i].subread
                if by=="ratio_short":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:# and bigg_list[i].score<sw_score:
               #         drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)

                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
                #            drop.add(i)
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                 #           drop.add(j)
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:# and bigg_list[j].score<sw_score:
                  #      drop.add(j)
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
    for i,j in reversed(ij_list):
      #  dict_read_num[bigg_list[i].name] = i
     #   print 'i j',i,j
        if D[i,j]<cutoff:
            if i==j:
                pass
            else:
                if by=="ratio":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:
           #             drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
      #                  print '1',i,j,bigg_list[i].name,bigg_list[j].name, bigg_list[j].subread
                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
            #                drop.add(i)
                          #  print '2',bigg_list[i].name,bigg_list[j].name
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
       #                     print '2',bigg_list[i].name,bigg_list[j].name,bigg_list[j]
                        else:
             #               drop.add(j)
        #                    print 3,bigg_list[j].name,bigg_list[i].name
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:
              #          drop.add(j)
                     #   print 4,bigg_list[j].name,bigg_list[i].name
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
         #               print 4,j,i,bigg_list[j].name,bigg_list[i].name,bigg_list[i].subread
                if by=="ratio_short":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:# and bigg_list[i].score<sw_score:
               #         drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)

                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
                #            drop.add(i)
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                 #           drop.add(j)
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:# and bigg_list[j].score<sw_score:
                  #      drop.add(j)
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

    for i,j in (ij_list):
      #  dict_read_num[bigg_list[i].name] = i
     #   print 'i j',i,j
        if D[i,j]<cutoff:
            if i==j:
                pass
            else:
                if by=="ratio":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:
           #             drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
      #                  print '1',i,j,bigg_list[i].name,bigg_list[j].name, bigg_list[j].subread
                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
            #                drop.add(i)
                          #  print '2',bigg_list[i].name,bigg_list[j].name
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
       #                     print '2',bigg_list[i].name,bigg_list[j].name,bigg_list[j]
                        else:
             #               drop.add(j)
        #                    print 3,bigg_list[j].name,bigg_list[i].name
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:
              #          drop.add(j)
                     #   print 4,bigg_list[j].name,bigg_list[i].name
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
         #               print 4,j,i,bigg_list[j].name,bigg_list[i].name,bigg_list[i].subread
                if by=="ratio_short":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:# and bigg_list[i].score<sw_score:
               #         drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)

                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
                #            drop.add(i)
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                 #           drop.add(j)
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:# and bigg_list[j].score<sw_score:
                  #      drop.add(j)
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
    for i,j in reversed(ij_list):
      #  dict_read_num[bigg_list[i].name] = i
     #   print 'i j',i,j
        if D[i,j]<cutoff:
            if i==j:
                pass
            else:
                if by=="ratio":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:
           #             drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
      #                  print '1',i,j,bigg_list[i].name,bigg_list[j].name, bigg_list[j].subread
                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
            #                drop.add(i)
                          #  print '2',bigg_list[i].name,bigg_list[j].name
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
       #                     print '2',bigg_list[i].name,bigg_list[j].name,bigg_list[j]
                        else:
             #               drop.add(j)
        #                    print 3,bigg_list[j].name,bigg_list[i].name
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:
              #          drop.add(j)
                     #   print 4,bigg_list[j].name,bigg_list[i].name
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
         #               print 4,j,i,bigg_list[j].name,bigg_list[i].name,bigg_list[i].subread
                if by=="ratio_short":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:# and bigg_list[i].score<sw_score:
               #         drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)

                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
                #            drop.add(i)
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                 #           drop.add(j)
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:# and bigg_list[j].score<sw_score:
                  #      drop.add(j)
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)


    '''
    keep = set()
  #  print 'hhappy' 
   # keep = fullset-drop-rep
    for infos in rep_label:
       if infos not in drop:

          keep.add(infos)
   #       print bigg_list[i].name
    kepp1 = fullset-cluster_reads
   # print cluster_reads
    for  infos in kepp1:
       # print 'keep',bigg_list[infos].name
        keep.add(infos)
    # change the default score of gene, no need to add
    for n, bigg in enumerate(bigg_list):
        if bigg.ttype=="isoform_anno":
            keep.add(n)
  #  print rep_label, drop
 #   print 'aaa',bigg_list[32].subread    
    # re_order D and bigg_list
    keepl = sorted(list(keep))
    bigg_list_new = select_list(bigg_list, keepl)
#    print 'keep',keepl
        
  #  for i in range(len(ij_list)):
#        print i,j, 'all' 
 #       print 'name',i,bigg_list[i].name 
#    print 'happy',dict_read_num
 #   print dict_read_num['4f3facf9-8dd9-48e5-8460-895ba248be39']
    '''
    for info  in keepl:
  #      print info, 
   #     print bigg_list[info].subread

        if len(bigg_list[info].subread)>1:
            for names in bigg_list[info].subread:
    #           print 'aaa', names, dict_read_num[names]
              if names in dict_read_num:  
                 label = dict_read_num[names] 
 
                 if label in bigg_list:
                     if len(bigg_list[label].subread)>1:
                        bigg_list[info].subread.add(bigg_list[label].name)
                        bigg_list[info].subread=bigg_list[info].subread.union(bigg_list[label].subread)
    '''
     #   print 'b', bigg_list[info].name
      #  print 'aloha',bigg_list[info].name,bigg_list[info].subread
    #----------------------------------------#
    #### sanity check for missed ones
    ## collect the missed ones
    pos_dic=get_pos_dic(bigg_list)
    missed_name=get_readall_bigg(bigg_list)-get_readall_bigg(bigg_list_new)

    #print missed_name
    #print "inside"
    #print len(get_readall_bigg(bigg_list_new))
    if add_miss:
        if len(missed_name)>0:
            #print "{} missing bigg found, added back but may affect the isoforms".format(len(missed_name))
            missed_num=set()
            for k in missed_name:
                missed_num.add(pos_dic[k])

            keep=keep.union(missed_num)

            keepl=sorted(list(keep))
            bigg_list_new=select_list(bigg_list, keepl)

            #### end of sanity check


    D = select_D(D, keepl)

    return D, bigg_list_new


