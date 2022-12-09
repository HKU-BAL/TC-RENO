# TC-RENO

# Overview <br>
TC-RENO is optimized using the basics of TrackCluster (https://github.com/Runsheng/trackcluster) with modification in the clustering method and addition of novel functions. TC-RENO is used for nanopore DRS data. <br>
Workflow:<br>
 ![example](http://www.bio8.cs.hku.hk/novel/AF-NS_workflow.png)<br> 

## Updates of the clustering method
1. The similarity matrix initialize value is changed to prevent two reads that do not intersect with each other from clustering together.  <br>
2. For reads supporting multiple isoforms, its contribution would be shared evenly among all supported isoforms across genes, i.e 1/num of supported isoforms.<br>
3. Each identified isoform has a fixed number of supporting reads, i.e, alteration of the input read order would not influence the clustering result.<br>
4. Denominator of Score 2 is updated from the minimum length of two intersecting isoforms to the length of each isoform, which can avoid merging some real isoforms together, such that isoforms with shorter exons length are identified as fragments and removed.<br>
Example:<br>
 ![example](http://www.bio8.cs.hku.hk/novel/isoform_exp.png)<br> 
   TC-RENO identifies isoform A and isoform B as two independent isoforms rather than consider isoform A as a fragment of isoform B.

## Function extension
1. TC-RENO introduces ade novo isoforms discovery module, such that  not  only isoforms with existing genes intersect could be identified.<br>
2. Not all short isoforms which are included in long isoforms would be discarded. This can toretain some short novel isoforms with more supported reads. TC-RENO classifies isoforms as “standard” and “high-confident” based on the number of support reads. Short isoforms with enough supporting reads (high-confident) would be saved, while other short isoforms would be regarded as fragments and discarded.<br>


# Quick start
## 1. Preprocess
Use `bam2bigg.py` of TrackCluster to convert bam files to bed file (reads.bed).<br>
Use `gff2bigg.py` of TrackCluster to convert reference gff files to bed file (ref.bed).<br>
Then, the overlaping reference regions will be merged together
``` 
bedtools sort -i ref.bed |bedtools merge -s -c 4 -o distinct -i - | awk '{OFS="\t"}{print $1, $2,$3,$1":"$2"-"$3, 0, $4}' > merged_ref.bed
```

## 2. Isoform identification and quantification
``` 
python TC_RENO.py -i reads.bed   -r merged_ref.bed -o isoforms.bed -q isoforms_exp.txt -tmp /tmp/ -t 2 
``` 
The identified isoforms are saved in `isoforms.bed` and the corresponding supporting reads are saved in `isoforms_exp.txt`.<br>

