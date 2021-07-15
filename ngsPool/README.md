# ngsPool
Estimation of allele frequencies and other metrics from pooled-sequencing data.

In this application, we show how `ngsJulia` can be used to estimate allele frequencies with an implementation suitable for low-coverage pooled-sequencing data.
When the sample size is given in addition to a gzipped mpileup file in input, `ngsPool` calculates allele frequency likelihoods which can be used for further applications. 
As an illustration, we provide scripts for estimating the site frequency and conduct a 2x2 association test.
'ngsPool` supports data filtering based on quality and depth values.

## Tutorial

Let's initialise paths to Julia language and `ngsJulia`.
```
JULIA=~/Software/julia-1.6.1/bin/julia
NGSJULIA=~/Software/ngsJulia
```

### Simulations 

We can simulate pooled-sequencing data in mpileup format by specifying several experimental parameters.
We can use the script provided:
```
Rscript $NGSJULIA/simulMpileup.R --help
```
using the `--pool 1` option

As an illustration, let's simulate 10 diploids of 100 base pairs with an average depth of 20 and base quality of 20 in Phred score, from a population of 50,000 effective size under constant-size evolution.
```
Rscript $NGSJULIA/simulMpileup.R --out test.txt --copy 2x10 --sites 100 --depth 20 --qual 20 --ksfs 1 --ne 50000 --pool | gzip > test.mpileup.gz

# true data
less -S test.txt

# observed sequencing data
less -S test.mpileup.gz
```

#### Case A: SNP calling and estimation of allele frequencies with a frequentist approach

We can look at all available options in `ngsPool`:
```
$JULIA $NGSJULIA/ngsPool/ngsPool.jl --help
```

If no information on the number of sampled chromosomes is provided, `ngsPool` returns a simple frequentist estimation of allele frequencies.
If a threshold on the likelihood ratio test for being a SNP is provided (`--lrtSnp`), all inferred monomorphic sites are filtered out.
These options can be achieved with the following command:
```
$JULIA $NGSJULIA/ngsPool/ngsPool.jl --fin test.mpileup.gz --fout test.out.gz --lrtSnp 6.64 
```
with `6.64` equivalent to a p-value of 0.01.

The program produces a file
```
less -S test.out.gz
```
with the following information for each site
* chrom: chromosome
* position: position
* reference: reference allele
* nonreference: inferred alternate allele
* major: inferred major allele
* minor: inferred minor allele
* lrtSnp: LRT for the site being a SNP
* lrtBia: LRT for the site being biallelic
* lrtTria: LRT for the site being triallelic
* maf: estimated minor allele frequency
* freqMax: maximum likelihood estimate of allele frequency (disabled if `--nChroms 0`)
* freqE: expected value of allele frequency (disabled if `--nChroms 0`)

#### Case B: estimation of allele frequencies without SNP calling using a likelihood approach

When the sample size is given with `--nChroms`, `ngsPool` returns additional estimates of allele frequency.
`nChroms` must be set equal to the product between the number of individuals and their ploidy.
With the `--fsaf`` option, the program will also produce a file with the per-site sample allele frequency likelihoods which can be exploited for further applications.
This option can be achieved with
```
$JULIA $NGSJULIA/ngsPool/ngsPool.jl --fin test.mpileup.gz --fout test.out.gz --nChroms 20 --fsaf test.saf.gz --minDepth 390

less -S test.out.gz
```
As an illustration, we also remove sites with a depth lower than 390.
Please note that columns `freqMax` and `freqE` are now enabled in the output file.

With this option, a new file is created containing the sample allele frequency likelihoods for each site:
```
less -S test.saf.gz
```
These likelihoods can be used for custom applications suitable for low-coverage pooled-sequencing data, as illustrated below.

#### Custom applications

#### Case C: site frequency spectrum

As an illustration, we provide a script to estimate the site frequency spectrum (SFS) from sample allele frequency likelihoods.
```
Rscript $NGSJULIA/ngsPool/scripts/poolSFS.R test.saf.gz > sfs.txt

less -S sfs.txt
```
The script returns three different estimates, as specified by the name on the first column:
* count: calculate the SFS by counting over the most likely allele frequency
* fit_count: fit an exponential distribution from the counts of the most likely allele frequency
* fit_saf: fit an exponential distribution from allele frequency likelihoods

The last estimate should be more suitable for low-coverage data. The distribution is fitted by minimising the KL divergence. 
Please note that fitted SFS do have have values for fixed allele frequencies.

#### Case D: association test

Assume we have two mpileup files for just one SNP, one file for cases and one for controls.
The SNP is associated to a binary trait and the allele frequencies are different between cases and controls.
```
# cases with a population allele frequency of 0.09
Rscript $NGSJULIA/ngsPool/scripts/simulMpileup_qq.R --copy 2x50 --sites 1 --depth 5 --qq 0.09 --pool | gzip > cases.mpileup.gz

$JULIA $NGSJULIA/ngsPool/ngsPool.jl --fin cases.mpileup.gz --nChroms 100 --fsaf cases.saf.gz --fout /dev/stderr

# controls with a population allele frequency of 0.04
Rscript $NGSJULIA/ngsPool/scripts/simulMpileup_qq.R --copy 2x50 --sites 1 --depth 5 --qq 0.04 --pool | gzip > controls.mpileup.gz

$JULIA $NGSJULIA/ngsPool/ngsPool.jl --fin controls.mpileup.gz --nChroms 100 --fsaf controls.saf.gz --fout /dev/stderr

```

The script provided will calculate a likelihood ratio test of the two allele frequencies being different. Note that this test is done on allele frequency likelihoods and not counts.
```
Rscript $NGSJULIA/ngsPool/scripts/poolAssoc.R cases.saf.gz controls.saf.gz
```

## Further options

All options available can be retrieved by:
```
$JULIA $NGSJULIA/ngsPool/ngsPool.jl --help
```
Several filtering options on base quality and depth. 
Likewise, options to include only bialleic sites (i.e. exclude triallelic and multiallelic sites) can be used.
Results may vary depending on the filtering options and users are encourage to consider how their inferences are robust to the data processing pipeline.







