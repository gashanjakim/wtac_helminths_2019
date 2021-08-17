# Genetic admixture notes


## F-statistics

### f3 (admixture)
- notes adapted from:
     - https://compvar-workshop.readthedocs.io/en/latest/contents/03_f3stats/f3stats.html


- F3 statistics are a useful analytical tool to understand population relationships.
- F3 statistics, just as F4 and F2 statistics measure allele frequency correlations between populations and were introduced by Nick Patterson in his Patterson 2012

- two main tests
     - test whether a target population (C) is admixed between two source populations (A and B)
     - measure shared drift between two test populations (A and B) from an outgroup (C)
- F3 statistics are in both cases defined as the product of allele frequency differences between population C to A and B, respectively:

```shell
F3(A,B;C)=⟨(c−a)(c−b)⟩
```
-  Here, ⟨⋅⟩ denotes the average over all genotyped sites, and a,b and c denote the allele frequency for a given site in the three populations A,B and C.

- interpretation
     - negative F3
          - allele frequency of the target population (C) is on average intermediate between the allele frequencies of A and B
          - if the entire statistics is negative, it suggests that in many positions, the allele frequency c is indeed intermediate, suggesting admixture between the two sources
     - positive F3
          - implies shared or correlated genetic drift between populations A and B

- additional notes:
     - if an F3 statistics is not negative, it does not proof that there is no admixture!
     - zscore
          - As general rule, a Z score of -3 or more suggests a significant rejection of the Null hypothesis that the statistic is not negative.




### f4 (also called D-statistics)







### f3 (outgroup)
- same basic concept to the admixture f3 statistics, however, instead of a target C and two source populations A and B, one now gives an outgroup C and two test populations A and B.
- the statistic F3(A, B; C) measures the branch length from C to the common ancestor of A and B, coloured red.
     - this statistic is simply a measure of how closely two population A and B are related with each other, as measured from a distant outgroup. It is thus a similarity measure: The higher the statistic, the more genetically similar A and B are to one another.
