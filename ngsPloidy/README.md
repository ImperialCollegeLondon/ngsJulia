# ngsPloidy
Inference of ploidy from short-read sequencing data.

In this application, we show how `ngsJulia` can be used to implement a statistical test of ploidy level assignment from sequencing data.
`ngsPloidy` takes in input gzipped mpileup files and returns several metrics of per-sample ploidy assignment. It also tests for aneuploidy within the sample. 
'ngsPloidy` supports data filtering based on quality and depth values.


## Model

With the following notation:
* P: probability
* D: observed sequencing data
* O: ploidy
* A: ancestral state

`ngsPloidy` calculates the following log-likelihood:
```
log(P(D, G, A|O_1=y_1, O_2=y_2, ..., O_n=y_n)) = sum_{samples} sum_{sites} P(D|G,A) P(G|O,A) P(A)
```
where:
* P(D|G,A): genotype likelihood
* P(G|O,A): genotype probability (derived from the expected population allele frequency from the site frequency spectrum)
* P(A): probability of correct assignment of ancestral state

The model assumes that we observed at most two alleles. 
It can also calculate the per-sample log-likelihood by omitting ```sum_{samples}```.
`ngsPloidy` will output the most likely array of marginal ploidies, as well as the likelihood of all samples having the same ploidy as a test against aneuploidy.

## Tutorial

### Initialisation

Let's initialise paths to Julia language and `ngsJulia`.
```
JULIA=~/Software/julia-1.6.1/bin/julia
NGSJULIA=~/Software/ngsJulia
```

Unless we wish to use uniform probabilities (as we will see later), first we need a create a file containing probabilities of genotypes and of the major allele being ancestral. These probability files can be generated using the following R script:
```
Rscript $NGSJULIA/ngsPloidy/writePars.R --help
```
where:
* '-k' denotes the shape of the site frequency spectrum:
        - k=1 : constant population size
        - k>1 : population bottleneck
        - k<1 : population growth
* '-n' is the effective population size
* '-s' if a flag and specifies that SNPs have beeen called; this makes sense only if only one samples is analysed with called SNPs
* '-p' specifies how to define the probability that the ancestral state is the major allele, if 0.5 this means you assume folded data, if -1 it will compute it using the site frequency spectrum
* '-h' prints a help message.

Therefore, we can create our file as:
```
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 1 > test.pars
# without ancestral allele (probability of being correct of 0.90)
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 0.90 > test.unk.pars
# with automatic calculation of probability of misassignment of ancestral allele
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p -1 > test.auto.pars
# with folded spectrum
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 0.5 > test.fold.pars
# with SNPs only
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 0.9 -n 100000 -p 1 -s > test.snp.pars
```

As an example,
```
cat test.pars
```
outputs in the first line the probabiolity of major allele being ancestral and derived, while in the remaining lines it gives the prior genotype probabilities for diallelic sites up to a maximum ploidy of 8.


### Simulations

We can simulate sequencing data in mpileup format by specifying the ploidy of each individual and other parameters of the sequencing experiment and species.
We can use the script provided:
```
Rscript $NGSJULIA/simulMpileup.R --help
```

### Case studies

Let's investigate several case studies to understand how `ngsPloidy` works.
We can run
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --help
```
to see all options available.

#### Case A: 2 haploids, 2 diploids, 2 triploids, 2 tetraploids, 2 pentaploids, sequenced a 20X at the haploid level

We simulate 100 sites all polymorphic in the population.
```
Rscript $NGSJULIA/simulMpileup.R --out test.A.txt --copy 1x2,2x2,3x2,4x2,5x2 --sites 100 --depth 20 | gzip > test.A.mpileup.gz

# true data
less -S test.A.txt

# observed sequencing data
less -S test.A.mpileup.gz
```

If we assume that we know the ancestral state and it's equivalent to the reference state, we can run:
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.A.mpileup.gz --fpars test.pars --nSamples 10 --keepRef 1
```
The option `--keepRef` forces the reference allele to be one of the two considered alleles and it is mandatory with `--fpars`.

Results are printed on the screen, and show:
* nr of analysed sites: vector of sites that passed filtering for each sample
* log-likelihoods of per-sample ploidies: a matrix of (nr_sample x nr_ploidies) with the log-likelihood of each sample having a certain ploidy (rows are separated by ;)
* #MLE vector of ploidies: the vector if individual maximum likelihood estimates of ploidy for each sample
* #log-likelikehood of MLE vector of ploidies: the log-likelihood of the above vector of estimated ploidies
* #LRT of aneuploidy: difference between MLE  vector of ploidies and the log-likliehood of all samples having the same ploidy, calculated for all tested ploidies

----------------------------------------

If `--fout` is given (optional), the program will print some statistics for each site, including the estimated allele frequency.
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.A.mpileup.gz --fpars test.pars --fout test.A.out.gz --nSamples 10 --keepRef 1

less -S test.A.out.gz
```
Specifically, this file reports:
* chrom: chromosome
* pos: position
* ref: reference allele
* depth: sequencing depth
* ref/anc: reference or ancestral allele (inferred)
* alt/der: alternate or derived allele (inferred)
* lrtSnp: LRT for the site being a SNP
* lrtBia: LRT for the site being biallelic
* lrtTria: LRT for the site being triallelic
* aaf: estimated alternate or ancestral allele frequency

---------------------------------------------

If `-fglikes` is given (optional), the program will ouput the per-site genotype likelihoods for each tested ploidy.
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.A.mpileup.gz --fpars test.pars --fglikes test.A.glikes.gz --nSamples 10 --keepRef 1

less -S test.A.glikes.gz
```

--------------------------

We can also change the genotype probabilities in input. For instance, we can impose an automatic setting of the probability of major allele being the ancestral state.
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.A.mpileup.gz --fpars test.auto.pars --fout test.A.out.gz --nSamples 10 --keepRef 1

less -S test.A.out.gz
```

---------------------------

One additional possibility is also to do not impose any genotype probability based on the site frequency spectrum P(G|O,A) (with `--fpars`) and use a uniform probability distribution instead (`--unif 1`).
With this option you are required that all sites where all the samples have data are processed (i.e. `--minSamples` should be equal to `--nSamples`), otherwise the program will throw an error.

```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.A.mpileup.gz --unif 1 --nSamples 10 --minSamples 10
```

---------------------------

Finally, we can even call genotypes and consider the likelihoods of called genotypes only (`-callGeno`), although this is not recommended for low-coverage sequencing data.

```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.A.mpileup.gz --unif 1 --callGeno 1 --nSamples 10 --minSamples 10
```


#### Case B: 10 triploids and no output file, simulating an error rate of assigning the ancestral state of 0.10

We can simulate such scenario with the following line:
```
Rscript $NGSJULIA/simulMpileup.R --out test.B.txt --copy 3x10 --sites 100 --depth 20 --panc 0.1 | gzip > test.B.mpileup.gz
```

As there is uncertainty in the polarisation, we can incorporate such error in the genotype probability file calcolated at the beginning.
As further illustration, we will filter out bases with a uqlaity score lower than 15 in Phred score.
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl	--fin test.B.mpileup.gz --fpars test.unk.pars --nSamples 10 --keepRef 1 --minQ 15
```
As you can see, the test for aneuploidy is rejected as the most likely vector of ploidies supports all samples being triploid. 

### Case C: 1 diploid, 8 triploids, 1 tetraploid with and folded allele frequencies

This scenario can be simulated with:
```
Rscript $NGSJULIA/simulMpileup.R --out test.C.txt --copy 2x1,3x8,4x1 --sites 100 --depth 20 | gzip > test.C.mpileup.gz
```

As the data is folded (we have no information on which allele is ancestral or derived), we can use the appropriate genotype probability file calculated above.
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.C.mpileup.gz --fpars test.fold.pars --nSamples 10 --keepRef 1
```

### Case D: 1 tetraploid with Ne=1e6 and experiencing population growth with SNP calling

In this case we have only one sample and we will attempt to infer its ploidy for called SNPs.
Let's assume we have 1000 sites polymorphic in the population. As we have one sample, we assume it is sequenced at higher depth.
This scenario can be simulated with (note that, as an illustration, we are not generating an output file for the real data):
```
Rscript $NGSJULIA/simulMpileup.R --copy 4x1 --sites 1000 --depth 30 --ksfs 0.90 --ne 100000 | gzip > test.D.mpileup.gz
```

Given that we wish to call SNPs, we need to indicate the appropriate file for genotype probabilities and a threshold for SNP calling (`--thSnp` in X^2 score value, here 6.64 equivalent to a p-value of 0.01).

```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.D.mpileup.gz --fpars test.snp.pars --keepRef 1 --nSamples 1 --thSnp 6.64
```
Please note that when `test.snp.pars` is used haploid likelihoods cannot be calculated.
Also note that in this case the number of analysed sites (SNPs) is less that the number of simulated sites.

## Further options

All options available can be retrieved by:
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --help
```
Several filtering options on base quality, depth and proportion of minor reads are available. Likewise, options to include only bialleic sites (i.e. exclude triallelic and multiallelic sites) can be used.
Results may vary depending on the filtering options and users are encourage to consider how their inferences are robust to the data processing pipeline.




