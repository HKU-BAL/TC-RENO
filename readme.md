# TC-RENO

# Pipeline for isoforms identification and quantification <br>
TC-RENO is optimized using the basics of TrackCluster (https://github.com/Runsheng/trackcluster) with modification in the clustering method and addition of novel functions. TC-RENO is used for nanopore DRS data. <br>

## Updates of the clustering method
1. The similarity matrix initialize value is changed to prevent two reads that do not intersect with each other from clustering together.  <br>
2. For reads supporting multiple isoforms, its contribution would be shared evenly among all supported isoforms across genes, i.e 1/num of supported isoforms.<br>
3. Each identified isoform has a fixed number of supporting reads, i.e, alteration of the input read order would not influence the clustering result.<br>
4. Denominator of Score 2 is updated from the minimum length of two intersecting isoforms to the length of each isoform, which can avoid merging some real isoforms together, such that isoforms with shorter exons length are identified as fragments and removed.<br>
Example:<br>
 ![Workflow](http://www.bio8.cs.hku.hk/RNA/tmp/isoform_github.png)<br> 


## Function extension
1. TC-RENO introduces ade novo isoforms discovery module, such that  not  only isoforms with existing genes intersect could be identified.<br>
2. Not all short isoforms which are included in long isoforms would be discarded. This can toretain some short novel isoforms with more supported reads. TC-RENO classifies isoforms as “standard” and “high-confident” based on the number of support reads. Short isoforms with enough supporting reads (high-confident) would be saved, while other short reads without enough support would be regarded as fragments and discarded.<br>


# Quick start
## 1. Preprocess
Use `bam2bigg.py` of TrackCluster to convert bam files to bed file (reads.bed).<br>
Use `gff2bigg.py` of TrackCluster to convert gff files to bed file (ref.bed).<br>
Noted: When identifying isoforms, the reference and reads are merged. So we use "." to distinguish reads and reference. The forth column of in ref.bed should include ".";  "." should not exist in the fourth column of the reads.bed file.

## 2. Isoform identification and quantification
``` 
python TC_RENO.py -i reads.bed   -r ref.bed -o isoforms.bed -q isoforms_exp.txt -tmp /dev/shm/tmp/
``` 
The identified isoforms are saved in cisoforms.bed` and the corresponding support reads are saved in `isoforms_exp.txt`.<br>

