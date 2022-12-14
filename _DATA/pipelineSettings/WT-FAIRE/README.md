# WT FAIRE Timecourse Details
### Pipeline

Based on McKay lab FAIRE pipeline...
```
git clone git@github.com:mckaylabunc/faire-pipeline.git
cd faire-pipeline
git checkout 911509a375a30712abba9723187843dd000f9bac
```

With modifications:

- dm6 blacklist from [Boyle Lab](https://github.com/Boyle-Lab)
  - [dm6-blacklist.v2.bed.gz](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/dm6-blacklist.v2.bed.gz)
  - Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z

- `bowtie2` run with `--very-sensitve` 

