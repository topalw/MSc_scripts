# Simple calculation of basic popgen stats


I created these simple functions that are in no way runtime-optimized because I was using SNPRelate for exploratory data analysis
and wanted to calculate these metrics without going through the painful step of file conversion.


## Content


There are 2 scripts in this folder. Both are R functions that REQUIRE the SNPRelate package to load the Genotype matrix from a GDS file. The functions require a START and FINISH index position that signify the population in the genotype matrix. For example if rows 12:24 are a  population, then start = 12 and fin = 24. The third argument is the PATH to the gds file in order to read the genotype matrix.
**care** if your individuals of the population are not serial in the gds you have to modify a couple of lines to accomodate intervals since the functions are made for ordered gds files. 

* One is called Heterozygosity_calc.r and it returns the average Observed and expected heterozygosity in a population from a genotype matrix, as well as the Standard errors of the mean. The number of observed heterozygosity is the number of heterozygotes over the number of individuals. The number of expected heterozygosity is the familiar HWE, 2pq with the allele frequencies calculated from the population data.

* The other script is called Private_alleles.r and returns the number of population specific alleles. It iterates over the columns and checks the sum of the rest of the rows (population excluded). If the sum is either 0 or 2*individuals it considers that the other allele is only found in the population of interest. The calculation takes into account NAs. 

### Misc

A cool way to check that the population-specific alleles function is correct, is to exclude your population of interest from a SNPRelate analysis (ex. Fst calculation), then the SNPRelate calculation excludes non-variant sites and says how many of them there are in your sample, thus how many private alleles the excluded population carried. 

Another important fact for begginers (like me) in SNPRelate is that the sample_id order is the same as the sample order in the vcf file that was used to produce the gds. To check that just 
```shell 
grep '^#CHROM' your.vcf | head -1
``` 

Best
