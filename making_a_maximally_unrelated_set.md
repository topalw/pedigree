Finding founders - Making a maximum independent set
================
Alexandros Topaloudis

# Making an unrelated set

Making an unrelated set of individuals in our pedigree structure can be
done through the beta kinship matrix in a graph approach. In this case
we set individuals as nodes in a graph and create edges between
individuals with beta kinship \> x where x a cutoff we define as kinship
allowed (ex. 0.0335). We then prune the graph to retain a set of
independent individuals. Simply this can be done by choosing one
individual at random, adding him to the independent set and then
removing all individuals related with them. Continuing through the
nodes, this would a give a maximal independent set but not the maximum
independent set (biggest possible maximal independent set). Multiple
iterations of the process each starting and progressing through the
nodes in different order would give multiple sets and after an
exhaustive search, the biggest set would be the maximum independent set.
Even with \~350 nodes this process is painstakingly slow. The problem is
solved in
![O(n^22^n)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;O%28n%5E22%5En%29 "O(n^22^n)")
by a brute force algorithm described above.

## Case 1 - manual modification and igraph

Lets start with the graph approach.

``` r
library(kinship2)
library(igraph)
library(knitr)
library(corrplot)
library(SNPRelate)
dp <- read.delim('depth_stats_all.idepth')
ped <- read.csv('genetically_fixed_pedigree.csv')
pedk <- kinship(ped$id, ped$dadid,ped$momid)
metadata <- read.csv('refpanel_metadata.csv')
b <- read.table('swiss_beta_matrix.table')
b <- b[-which(rownames(b)=='877825F2'),-which(rownames(b)=='877825F2')]
rownames(b) <- metadata$NEWname[match(rownames(b),metadata$rawVCFname)]
kept <- read.csv('individuals_to_keep.csv')
b <- b[match(kept$NEWname,rownames(b)),match(kept$NEWname,rownames(b))]
```

We can make a graph using an adjacency matrix we basically translate a
kinship matrix with values \> cutoff –\> 1 and values less than a cutoff
to 0. Pairs of individuals with 1 in the matrix represent nodes connect
with an edge while pairs with 0 are unconnected.

``` r
kin2adj <- function(bmat, cutoff){
  edges <- data.frame(which(bmat > cutoff,arr.ind = T))
  tmp <- matrix(data=rep(0,nrow(bmat)^2), nrow=nrow(bmat), ncol=nrow(bmat))
  for(i in 1:nrow(edges)){tmp[edges$row[i],edges$col[i]]<- 1}
  diag(tmp) <- 0
  return(tmp)
}
g1 <- graph_from_adjacency_matrix(kin2adj(b,0.2),mode = 'undirected')
plot(g1)
```

![](making_a_maximally_unrelated_set_files/figure-gfm/make%20a%20graph-1.png)<!-- -->

``` r
#mis <- largest_ivs(g1) # prohibitingly slow !!!
```

This function is as expected super slow so lets use the greedy algorithm
of kinship2 and then verify with the beta matrix.

## Case 2 - kinship2 function pedigree.unrelated() + extra prunning

Kinship2 has a function that finds unrelated individuals in a pedigree
using a greedy algorithm. However, this only gives a list of individuals
with kinship == 0 using a pedigree input. This means that we need to
provide a pedigree to the function and this pedigree should capture all
our links. By construction the best we can do is the improved pedigree
with genetic kinship. Lets see what this gives us…

``` r
unrel1  <- pedigree.unrelated(pedigree(ped$id,ped$dadid,ped$momid,ped$sex),avail=ped$seq)
```

    ## Warning in pedigree(ped$id, ped$dadid, ped$momid, ped$sex): More than 25% of the
    ## gender values are 'unknown'

``` r
bu <- as.matrix(b[rownames(b) %in% unrel1,rownames(b) %in% unrel1])
image(bu)
```

![](making_a_maximally_unrelated_set_files/figure-gfm/pedigree%20unrelated-1.png)<!-- -->

``` r
diag(bu) <- 0
hist(as.numeric(as.matrix(bu)),breaks=100) 
```

![](making_a_maximally_unrelated_set_files/figure-gfm/pedigree%20unrelated-2.png)<!-- -->

We still see that some individuals have links \> 0 which is due to the
incompleteness of our pedigree we fed to the kinship2 method. Now we can
use a heuristic algorithm that removes individuals based on the \# of \>
x links they have in the genomic matrix. At the same time we can take
into account another layer of information ex. Depth and remove
individuals sorted on \# of links and depth to make sure we keep the
highest depth individuals when possible.

``` r
find.names <- function(beta.mat,cutoff){
  df <- data.frame('id1'=NULL,
                   'id2'=NULL,
                   'beta'=NULL)
  end <- nrow(beta.mat) - 1
  for(i in 1:end){
    start <- i + 1
    for(j in start:nrow(beta.mat)){
        if(beta.mat[i,j] >= cutoff){
          tmp.df <- data.frame('id1'=rownames(beta.mat)[i],
                               'id2'=rownames(beta.mat)[j],
                               'beta'=beta.mat[i,j])
          df <- rbind(df,tmp.df)
        }
      }
  }
  return(df)
}
make.ind.df <- function(pairs){
  tmp <- data.frame('id'=unique(c(pairs$id1,pairs$id2)))  # make individual df
  # make negative to match decreasing order !!!
  tmp$dp <- -dp$MEAN_DEPTH[match(metadata$rawVCFname[match(tmp$id,metadata$NEWname)],dp$INDV)] # add dp
  tmp$links <- 0 #  add link count
  for(i in 1:nrow(tmp)){
  tmp$links[i] <- sum(c(pairs$id1,pairs$id2)==tmp$id[i])
  }
  tmp <- tmp[order(tmp[,3],tmp[,2], decreasing = T),]
  return(tmp)
}
bu <- as.matrix(b[rownames(b) %in% unrel1,rownames(b) %in% unrel1])
pairs <- find.names(bu,0.1) 
#kable(pairs,caption = 'pairs of unrelated that are related!') 
related <- make.ind.df(pairs)
kable(related[1:10,],caption='related individuals with # of links' )
```

|     | id      |        dp | links |
|:----|:--------|----------:|------:|
| 11  | M032051 |  -9.63268 |     3 |
| 13  | M033660 |  -9.44783 |     2 |
| 16  | M032002 | -29.14100 |     2 |
| 17  | M026459 |  -8.29615 |     1 |
| 15  | M017312 |  -8.73039 |     1 |
| 20  | M035431 |  -9.13105 |     1 |
| 6   | M006413 |  -9.44743 |     1 |
| 10  | M031581 |  -9.69714 |     1 |
| 4   | M004500 | -10.05470 |     1 |
| 7   | M011816 | -10.25360 |     1 |

related individuals with \# of links

Now we can prune this list by removing the individuals with the biggest
\# of pairs and then recalculating this list and then removing again
until list is empty.

``` r
bu <- as.matrix(b[rownames(b) %in% unrel1,rownames(b) %in% unrel1])
find.unrelated <- function(b,cutoff){
  bu <- b
  # initialize pairs   
  pairs <- find.names(b,cutoff) # get pairs 
  while(nrow(pairs)!=0){
    # find and remove worst candidate
    tmp <- make.ind.df(pairs)
    rm.ind <- tmp[1,1]
    bu <- bu[-match(rm.ind,rownames(bu)),-match(rm.ind,rownames(bu))]
    # recalculate pairs
    pairs <- find.names(bu,cutoff)
  }
  return(bu)
}
bu <- find.unrelated(bu,0.01)
dim(bu)
```

    ## [1] 41 41

``` r
hist(as.numeric(bu),breaks=140,xlim=c(-.1,.1),main='histogram of genomic kinship of unrelated samples',
     xlab='beta')
abline(v=0.0335,lty=2)
diag(bu) <- 0
corrplot(bu, tl.pos='n',diag = F,method = 'shade',is.corr = F)
#unrelated1 <- metadata[match(rownames(bu),metadata$NEWname),] # ~75 at 0.0335
#write.csv(unrelated1, 'unrelated_001_kinship_42_040722.csv',quote=F,row.names = F)
```

![](making_a_maximally_unrelated_set_files/figure-gfm/sequential%20removal-1.png)![](making_a_maximally_unrelated_set_files/figure-gfm/sequential%20removal-2.png)

This gives us a stringent set of unrelated individuals. Lowering the
limit in the manual part in the end can make the cutoff more stringent.
We choose 0.03125 (- 0.0335 a bit more to allow for variance in beta)
which represents 4th degree distant relatedness and prunes relatedness
well enough without reducing the sample size significantly. The code
works and is very fast (not our code - this one is slow but thankfully
it only has to prune a few individuals - since kinship2 does most of the
work)…

## Europe.

And for the Europeans we can follow the same procedure, but here a
strict cutoff on kinship will not be representative of true kinship
(kinship was normalized on average of whole dataset but family structure
and overrepresentation of CH has an impact). Specifically, the mean
relatedness is a bit higher than expected from an unrelated sample and
this means that we expect more negative values than normal. This should
not impact our cutoff.

``` r
geno <- snpgdsOpen('ref_panel_snps_f1_masked_maf05_miss05_LDpruned.gds')
beu <- snpgdsIndivBeta(geno,autosome.only = F)
```

    ## Individual Inbreeding and Relatedness (beta estimator):
    ## Excluding 0 SNP (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
    ##     # of samples: 508
    ##     # of SNPs: 85,246
    ##     using 1 thread
    ## Individual Beta:    the sum of all selected genotypes (0,1,2) = 51312622
    ## CPU capabilities: Double-Precision SSE2
    ## Mon Jul 11 09:54:08 2022    (internal increment: 65536)
    ## [..................................................]  0%, ETC: ---        [==================================================] 100%, completed, 1s
    ## Mon Jul 11 09:54:09 2022    Done.

``` r
beur <- snpgdsIndivBetaRel(beu,mean(beu$beta))
indx1 <- match(metadata$rawVCFname[metadata$pop !='CH'], beur$sample.id)
eub <- beur$beta[indx1, indx1]
rownames(eub) <- beur$sample.id[indx1]
nrow(eub)
```

    ## [1] 156

``` r
hist(eub,breaks=100)
abline(v=0.2)
abline(v=0.1)
abline(v=0.0335)
nrow(find.unrelated(eub,0.2))
```

    ## [1] 150

``` r
nrow(find.unrelated(eub,0.12))
```

    ## [1] 122

``` r
nrow(find.unrelated(eub,0.0335))
```

    ## [1] 50

![](making_a_maximally_unrelated_set_files/figure-gfm/europeans%20unrelated-1.png)<!-- -->
