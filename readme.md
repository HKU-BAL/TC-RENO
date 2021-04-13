# TC-RENO

# Pipeline for isoforms identification and quantification <br>
TC-RENO modifies the clustering method of TrackCluster and expands the existing functions. <br>

## Update the clustering method
1. Existing isoforms reported in Araport11 would be clustered together even if the similarity score between the two isoforms is below the cut-off.
2. The clustering result varied by the ordering of reads. Alteration of the read order results in a different number of supported reads (i.e. subreads) in the final identified isoform.
3. Wrong method of computing isoforms’ supporting reads
4. TrackCluster mistakenly identifies some real isoforms as fragments and discarded them.


## Function extension:
1. The defaut TrackCluster only identifies novel isoforms that intersect with the existing gene annotation. For those remaining reads that did not align to the existing isoforms, we developed a new function, de novo identification, to further cluster them together.
2. We would save short isoforms with enough supporting reads rather than discard them directly.
