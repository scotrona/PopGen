---
title: "HW1.R"
output: 
  html_document: 
    keep_md: yes
name: "Sierra Cotrona Gentry"
---


Chapter 4

Population genetics refers to expectations, probabilities, and the observed frequencies in relation to those expectations. R is a useful tool in this. 


```r
p <- 0.5 
p^2 
```

```
## [1] 0.25
```

```r
f <- 1/1000 
sqrt (f) 
```

```
## [1] 0.03162278
```

```r
1-(sqrt(f)) 
```

```
## [1] 0.9683772
```
We can set an allele frequency, define the frequency of homozygotes, introduce a mutant allele, and understand the probability of inheriting that allele.

```r
curve (x^2, 0, 1, 
       xlab = "Allele Frequencies", ylab = "Genotype Frequencies", 
       col = "green", lwd = 2)
text (0.6, 0.2, "Homozygotes", col = "green")
curve (2*x*(1-x), 0, 1, add = TRUE, 
       xlab = "Allele Frequencies", ylab = "Genotype Frequencies", 
       col = "blue", lwd = 2) 
text (0.25, 0.5, "Heterozygotes", col = "blue")
p <- 2/3
points (p, p^2, lwd=2, cex=2)
```

![](Project-1_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

Graphically represented is the expected homozygote and heterozygote genotype frequencies as a function of allele frequencies. There is a defined point on the graph where homozygotes and heterozygotes are equal in frequency. 


```r
allele <- c("A", "A", "a", "a", "a", "a", "a", "a", "a", "a")
allele <- c(rep ("A", 2), rep ("a", 8))
popsize <- 100
pop <- matrix (nrow = popsize, ncol = 2) 
for (i in 1: popsize) {
  pop [i, 1] <- sample (allele, 1)
  pop [i, 2] <- sample (allele, 1)
} #for loop uses variable i to keep track of position from one to popsize
head (pop) #shows heterozygotes/homozygotes for each individual
```

```
##      [,1] [,2]
## [1,] "A"  "a" 
## [2,] "a"  "A" 
## [3,] "a"  "a" 
## [4,] "a"  "a" 
## [5,] "a"  "a" 
## [6,] "A"  "a"
```
Simulating data is useful for understanding how data can be manipulated in R. We have simulated a list of alleles to create a matrix with 100 individuals as an initial population size. The loop is to randomly choose alleles to make up genotypes (2 copies of alleles). The rows are the individual and their corresponding columns are their two alleles making up their genotype. 


```r
Acount <-0
for (i in 1: popsize) {
  if (pop [i, 1] == "A") {
    Acount <- Acount +1
  }
  if (pop [i, 2] == "A") {
    Acount <- Acount +1
  }
} 
AFreq <- Acount/(popsize*2)
AFreq 
```

```
## [1] 0.205
```
To calculate A allele frequencies per run, a variable is defined and a loop is created asking if the allele in each position is equal to A. If it is, it will be added to the "Acount," which is then divided by the total number of allele copies present to find the frequency. 


```r
Hcount <- 0
AAcount <- 0
for (i in 1: popsize) {
  if (pop [i, 1] == "A") {
  if (pop [i, 2] == "a") {
    Hcount <- Hcount +1
  }else {
    AAcount <- AAcount +1
  }
  }
  if (pop [i, 1] == "a") {
    if (pop[i, 2] == "A") {
      Hcount <- Hcount + 1
    }
  }
}
HetFreq <- Hcount/(popsize)
AAFreq <- AAcount/(popsize)
print(c(AFreq, HetFreq, AAFreq))
```

```
## [1] 0.205 0.350 0.030
```
Genotype frequencies may also be calculated through this matrix by counting the number of heterozygotes and p-homozygotes. The loop specifies that there are a multitude of pairings it is looking for: a, A; A, a; and A, A. 

```r
plot (AFreq, HetFreq, 
      xlab = "allele frequency", ylab = " ",
      ylim = c(0,1), xlim = c(0,1), col = "blue")
par (new = TRUE)
plot (AFreq, AAFreq, xlab = " ", ylab = " ", 
      ylim = c(0,1), xlim = c(0,1), col = "green" )
par (new = TRUE)
plot (AFreq, 1-AAFreq-HetFreq, 
      xlab = " ", ylab = "genotype frequency", 
      ylim = c(0,1), xlim = c(0,1), col = "red")
curve (2*x*(1-x), 0, 1, 
       add = TRUE, ylab = NULL, lwd= 2, 
       ylim = c(0, 1), col = "darkblue")
curve (x**2, 0, 1, 
       add = TRUE, ylab = NULL, lwd= 2, 
       ylim = c(0, 1), col = "darkgreen")
curve ((1-x)**2, 0, 1, 
       add = TRUE, ylab = NULL, lwd= 2, 
       ylim = c(0, 1), col = "darkred")
text (0.5, 0.7, "Aa", col = "blue")
text (0.9, 0.7, "AA", col = "green")
text (0.1, 0.7, "aa", col = "red")
```

![](Project-1_files/figure-html/unnamed-chunk-6-1.png)<!-- -->
Single points are plotted in the above graph, which represents the same data we are working with. The curves represent Hardy-Weinberg predictions, and the points fall not far from them. 

```r
AFreq <- numeric ()
HetFreq <- numeric()
AAFreq <- numeric()
replicates <- 10 
for (j in 1: replicates) {
  for (i in 1: popsize) {
    pop [i, 1] <- sample (allele, 1)
    pop [i, 2] <- sample (allele, 1)
  }
  Acount <- 0 
  for (i in 1: popsize) {
    if (pop [i, 1] == "A") {
      Acount <- Acount +1
    }
    if (pop [i, 2] == "A") {
      Acount <- Acount +1 
  }
  }
  AFreq [j] <- Acount/(popsize*2)
  Hcount <- 0
  AAcount <- 0
  for (i in 1: popsize) {
    if (pop [i, 1] == "A") {
      if (pop [i, 2] == "a") {
        Hcount <- Hcount +1
      }else {
        AAcount <- AAcount +1
      }
    }
    if (pop [i, 1] == "a") {
      if (pop[i, 2] == "A") {
        Hcount <- Hcount + 1
      }
    }
  }
  HetFreq [j] <- Hcount/(popsize)
  AAFreq [j] <- AAcount/(popsize)
}
plot (AFreq, HetFreq, 
      xlab = "allele frequency", ylab = " ",
      ylim = c(0,1), xlim = c(0,1), col = "blue")
par (new = TRUE)
plot (AFreq, AAFreq, xlab = " ", ylab = " ", 
      ylim = c(0,1), xlim = c(0,1), col = "green" )
par (new = TRUE)
plot (AFreq, 1-AAFreq-HetFreq, 
      xlab = " ", ylab = "genotype frequency", 
      ylim = c(0,1), xlim = c(0,1), col = "red")
curve (2*x*(1-x), 0, 1, 
       add = TRUE, ylab = NULL, lwd= 2, 
       ylim = c(0, 1), col = "darkblue")
curve (x**2, 0, 1, 
       add = TRUE, ylab = NULL, lwd= 2, 
       ylim = c(0, 1), col = "darkgreen")
curve ((1-x)**2, 0, 1, 
       add = TRUE, ylab = NULL, lwd= 2, 
       ylim = c(0, 1), col = "darkred")
text (0.5, 0.7, "Aa", col = "blue")
text (0.9, 0.7, "AA", col = "green")
text (0.1, 0.7, "aa", col = "red")
```

![](Project-1_files/figure-html/unnamed-chunk-7-1.png)<!-- -->
The graph is similar as before, but this time the data was run through 10 simulations. 

```r
p <- runif (1)
replicates <- 1000
for (j in 1: replicates) {
  for (i in 1: popsize) {
    if(runif(1)<p){
      pop [i, 1] <- "A"
    } else {
      pop [i, 1] <- "a"
    }
    if(runif(1)<p){
      pop [i, 2] <- "A"
    } else {
      pop [i, 2] <- "a"
    }
  }
  Acount <- 0 
  for (i in 1: popsize) {
    if (pop [i, 1] == "A") {
      Acount <- Acount +1
    }
    if (pop [i, 2] == "A") {
      Acount <- Acount +1 
    }
  }
  AFreq [j] <- Acount/(popsize*2)
  Hcount <- 0
  AAcount <- 0
  for (i in 1: popsize) {
    if (pop [i, 1] == "A") {
      if (pop [i, 2] == "a") {
        Hcount <- Hcount +1
      }else {
        AAcount <- AAcount +1
      }
    }
    if (pop [i, 1] == "a") {
      if (pop[i, 2] == "A") {
        Hcount <- Hcount + 1
      }
    }
  }
  HetFreq [j] <- Hcount/(popsize)
  AAFreq [j] <- AAcount/(popsize)
}
plot (AFreq, HetFreq, 
      xlab = "allele frequency", ylab = " ",
      ylim = c(0,1), xlim = c(0,1), col = "blue")
par (new = TRUE)
plot (AFreq, AAFreq, xlab = " ", ylab = " ", 
      ylim = c(0,1), xlim = c(0,1), col = "green" )
par (new = TRUE)
plot (AFreq, 1-AAFreq-HetFreq, 
      xlab = " ", ylab = "genotype frequency", 
      ylim = c(0,1), xlim = c(0,1), col = "red")
curve (2*x*(1-x), 0, 1, 
       add = TRUE, ylab = NULL, lwd= 2, 
       ylim = c(0, 1), col = "darkblue")
curve (x**2, 0, 1, 
       add = TRUE, ylab = NULL, lwd= 2, 
       ylim = c(0, 1), col = "darkgreen")
curve ((1-x)**2, 0, 1, 
       add = TRUE, ylab = NULL, lwd= 2, 
       ylim = c(0, 1), col = "darkred")
text (0.5, 0.7, "Aa", col = "blue")
text (0.9, 0.7, "AA", col = "green")
text (0.1, 0.7, "aa", col = "red")
```

![](Project-1_files/figure-html/unnamed-chunk-8-1.png)<!-- -->
And here the data are simulated with 1000 replicates. A thousand possible outcomes are plotted on the above graph. 

```r
p <- (9+842/2)/23369
p
```

```
## [1] 0.01840045
```

```r
2*p*(1-p)
```

```
## [1] 0.03612374
```

```r
curve (2/x, 1e-7, 0.01, log = "y")
```

![](Project-1_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
library ("popgenr")
data ("snp")
str (snp)
```

```
## 'data.frame':	25 obs. of  6 variables:
##  $ ID        : Factor w/ 25 levels "rs1000000","rs1100000",..: 10 20 23 24 25 1 2 3 4 5 ...
##  $ p         : num  0.364 0.281 0.381 0.023 0.055 0.166 0.406 0.82 0.156 0.012 ...
##  $ hom       : num  0.145 0.08 0.157 0.002 0.005 0.027 0.167 0.734 0.026 0.002 ...
##  $ het       : num  0.439 0.403 0.448 0.042 0.099 0.278 0.478 0.171 0.261 0.02 ...
##  $ chromosome: int  5 16 8 22 6 12 21 13 14 16 ...
##  $ type      : Factor w/ 7 levels "3UTR","downstream",..: 3 3 7 7 4 7 4 6 4 3 ...
```

```r
plot (snp$type)
```

![](Project-1_files/figure-html/unnamed-chunk-9-2.png)<!-- -->
Plotted here are observational datafor different locations of SNPs in the dataset defined above. 
Allele and genotype frequencies can be plotted for this data as well. 

```r
curve (x^2, 0, 1, 
       xlab = "Allele frequencies", ylab = "Genotype frequencies",
       col = "green", lwd = 2)
text (0.6, 0.2, "Homozygotes", col = "green") 
curve (2*x* (1-x), 0, 1, add = TRUE,
       xlab = "Allele frequencies", ylab = "Genotype frequencies",
       col = "blue", lwd = 2)
text (0.25, 0.5, "Heterozygotes", col = "blue") 
points (snp$p, snp$hom, pch = 19, col = "green")
points (snp$p, snp$het, pch = 19, col = "blue")
```

![](Project-1_files/figure-html/unnamed-chunk-10-1.png)<!-- -->
Predictions and measured data appear to be mostly concordant. The deviant data behave as expected, toward homozygosity. 

Chapter 5

The purpose of the chapter was to use chi-square tests to see if observed genotype frequencies are diverging from our Hardy-Weinberg predictions. 

The first is a dataset of genotype counts from 501 people in Lagos, Nigeria of a hemoglobin-producing gene associated with sickle-cell anemia. 

```r
AA <- 366
AS <- 123
SS <- 12
n <- AA + AS + SS
print (paste ("n:", n))
```

```
## [1] "n: 501"
```

```r
p <- (SS + (AS/2))/n
print (paste ("p:", p))
```

```
## [1] "p: 0.146706586826347"
```

```r
EAA <- n*(1-p)^2
EAS <- n*2*p*(1-p)
ESS <- n*p^2
print (paste ("Expected:", EAA, EAS, ESS))
```

```
## [1] "Expected: 364.782934131737 125.434131736527 10.7829341317365"
```

```r
chi2 <- (EAA - AA)^2/EAA + 
  (EAS-AS)^2/EAS +
  (ESS-SS)^2/ESS
print (paste ("chi-square:", chi2))
```

```
## [1] "chi-square: 0.188666341317465"
```

```r
pvalue <- pchisq(chi2, df = 1, lower.tail = FALSE)
print (paste ("P-value:", pvalue))
```

```
## [1] "P-value: 0.664028935603698"
```

```r
geno <- c (AA, AS, SS)
expe <- c (EAA, EAS, ESS)
G <- 2 * sum (geno * log (geno/expe))
print (paste("G:", G))
```

```
## [1] "G: 0.184075382936793"
```

```r
pvalue <- pchisq(G, df = 1, lower.tail = FALSE)
print (paste ("P-value:", pvalue))
```

```
## [1] "P-value: 0.667894063332523"
```
Allele values are given, and allele and genotype frequencies are found. Then a chi-square test is performed with a resulting p-value to indicate whether the population is near Hardy-Weinberg expectations. 
These values indicate that our observed data is close to expected Hardy-Weinberg values. Large p-values indicate that the two values are not statistically different.

Below is the same statistical analysis, but the sample size is smaller and SS allele is zero. 

```r
AA <- 34
AS <- 100
SS <- 0
n <- AA + AS + SS
p <- (SS + (AS/2))/n
(EAA <- n*(1-p)^2)
```

```
## [1] 52.65672
```

```r
(EAS <- n*2*p*(1-p))
```

```
## [1] 62.68657
```

```r
(ESS <- n*p^2)
```

```
## [1] 18.65672
```

```r
geno <- c (AA, AS, SS)
expe <- c (EAA, EAS, ESS)
chi2 <- sum ((expe-geno)^2/expe)
print (paste ("chi square:", chi2))
```

```
## [1] "chi square: 47.4773242630385"
```

```r
pvalue <- pchisq(chi2, df = 1, lower.tail = FALSE)
print (paste ("P-value:", pvalue))
```

```
## [1] "P-value: 5.56438965556751e-12"
```
This example does not meet the expected values given with Hardy-Weinberg. This is an indicated that there is a variant in the gene that is perhaps being acted upon by an evolutionary force. It is not random, as Hardy-Weinberg suggests. 

```r
dat <- matrix (c(geno, expe), nrow = 2, byrow = T)
barplot (dat, beside = T, 
         col = c ("turquoise4", "sienna1"),
         names.arg = c ("AA", "SA", "SS"))
legend (x = "topright", legend = c ("Observed", "Expected"), 
        pch = 15, col = c ("turquoise4", "sienna1") )
```

![](Project-1_files/figure-html/unnamed-chunk-13-1.png)<!-- -->
Visually, it is clear that there is a discrepancy between observed and expected results.  

When genotypes are extended to more than two alleles, they may easily be represented in R as follows: 

```r
p1 <- 0.2
p2 <- 0.3
p3 <- 0.5
sapply (c(p1, p2, p3), function (x) x^2)
```

```
## [1] 0.04 0.09 0.25
```

```r
sum (sapply (c(p1, p2, p3), function (x) x^2))
```

```
## [1] 0.38
```
Homozygote genotype frequency may still be determined with multiple allele frequencies. 

Data below uses microsatellites, segments of DNA with short repeat blocks  that tend to be highly variable. The data is of 9 loci, with two copies from each individual at each loci. 

```r
library ("popgenr")
data (genotypes)
rownames (genotypes) <- genotypes$ID
genotypes <- genotypes [, -c (1,2)]
(num.loci <- (length(genotypes))/2)
```

```
## [1] 9
```

```r
Hom_exp <- NULL
Het_exp <- NULL
Hom_obs <- NULL
Het_obs <- NULL
for (n in 1: (num.loci)) {
  current <- n*2-1
  locus <- c(genotypes [,current], genotypes [,current+1])
  alleles <- unique (locus)
  alleles <- alleles [alleles!=-1]
  p_allele <- NULL
  for (a in 1: length (alleles)) {
    p_allele <- c(p_allele,
                  sum (alleles [a] == locus)/sum (locus!=-1))
  }
  Hom_exp <- c(Hom_exp, sum(sapply (p_allele, 
                                    function (x) x^2)))
  obs <- 0
  for (i in 1: length (genotypes [,current])){
    if (genotypes [i, current]!=-1) {
      if (genotypes [i, current] ==
          genotypes [i, current + 1]) {
        obs <- obs + 1
      }
    }
  }
  Hom_obs <- c (Hom_obs, obs/(sum (locus!=-1)/2))
}
(Het_exp <- 1 - Hom_exp)
```

```
## [1] 0.22163508 0.60444977 0.26943371 0.60637435 0.03982418 0.47140832 0.41982455
## [8] 0.81408284 0.53234569
```

```r
(Het_obs <- 1 - Hom_obs)
```

```
## [1] 0.23364486 0.70895522 0.25000000 0.58823529 0.04065041 0.45652174 0.45569620
## [8] 0.78461538 0.53378378
```
The created loop is for the purpose of performing an allele count at each loci, excluding any missing data. It is meant to find both expected heterozygous genotype frequency and observed heterozygous genotype frequency. 

Detailed plots of the heterozygous genotype frequencies are below:

```r
plot (Het_obs, Het_exp)
abline (lm (Het_exp ~ Het_obs))
reg <- summary (lm (Het_exp ~ Het_obs))
print (reg)
```

```
## 
## Call:
## lm(formula = Het_exp ~ Het_obs)
## 
## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.087500 -0.011401  0.009525  0.023184  0.049083 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) 0.007452   0.032044   0.233    0.823    
## Het_obs     0.965502   0.063609  15.179  1.3e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.04312 on 7 degrees of freedom
## Multiple R-squared:  0.9705,	Adjusted R-squared:  0.9663 
## F-statistic: 230.4 on 1 and 7 DF,  p-value: 1.296e-06
```

```r
rr <- reg$r.squared
rrlabel <- paste ("r-squared =", round (rr, digits = 3))
pv <- reg$coefficients [2,4]
pvlabel <- paste ("P-value =", pv)
text (0.6, 0.2, rrlabel)
text (0.6, 0.15, pvlabel)
```

![](Project-1_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

Blood type allele frequencies can be determined using phenotype data as well. Although more complex, the same applications as above are used. 

```r
AB <- 5
A <- 30
B <- 7
O <- 36
(N <- AB + A + B + O)
```

```
## [1] 78
```

```r
A/N
```

```
## [1] 0.3846154
```

```r
B/N
```

```
## [1] 0.08974359
```

```r
AB/N
```

```
## [1] 0.06410256
```

```r
O/N
```

```
## [1] 0.4615385
```

```r
(Pi <- sqrt(O/N))
```

```
## [1] 0.6793662
```

```r
(Pa <- sqrt ((A+O)/N)-Pi)
```

```
## [1] 0.2405
```

```r
(Pb <- sqrt ((B+O)/N)-Pi)
```

```
## [1] 0.06311748
```

```r
signif (Pa + Pb + Pi, 1)
```

```
## [1] 1
```
The function has been forced to one although the populations have deviated from Hardy-Weinberg, but this does not give us true allele frequencies. We are required to modify estimates to maximum likelihood solutions. 

```r
(Paa <- A*(Pa^2/((Pa^2)+2*(Pa*Pi))))
```

```
## [1] 4.511539
```

```r
(Pai <- A * (2*(Pa*Pi))/((Pa^2)+(2*(Pa*Pi))))
```

```
## [1] 25.48846
```

```r
(Pbb <-B*(Pb^2/((Pb^2)+(2*(Pb*Pi)))))
```

```
## [1] 0.3107377
```

```r
(Pbi <- B*(2*(Pb*Pi))/((Pb^2)+(2*(Pb*Pi))))
```

```
## [1] 6.689262
```

```r
(Pii <- O)
```

```
## [1] 36
```

```r
(Pab <- AB)
```

```
## [1] 5
```

```r
(Pa <- ((2*Paa)+Pai+Pab)/(2*N))
```

```
## [1] 0.2532791
```

```r
(Pb <- ((2*Pbb)+Pbi+Pab)/(2*N))
```

```
## [1] 0.07891499
```

```r
(Pi <- ((2*Pii)+Pai+Pbi)/(2*N))
```

```
## [1] 0.6678059
```

```r
(Paa <- A*(Pa^2/((Pa^2)+2*(Pa*Pi))))
```

```
## [1] 4.782187
```

```r
(Pai <- A * (2*(Pa*Pi))/((Pa^2)+(2*(Pa*Pi))))
```

```
## [1] 25.21781
```

```r
(Pbb <-B*(Pb^2/((Pb^2)+(2*(Pb*Pi)))))
```

```
## [1] 0.3905227
```

```r
(Pbi <- B*(2*(Pb*Pi))/((Pb^2)+(2*(Pb*Pi))))
```

```
## [1] 6.609477
```

```r
(Paa <- A*(Pa^2/((Pa^2)+2*(Pa*Pi))))
```

```
## [1] 4.782187
```

```r
(Pai <- A * (2*(Pa*Pi))/((Pa^2)+(2*(Pa*Pi))))
```

```
## [1] 25.21781
```

```r
(Pbb <-B*(Pb^2/((Pb^2)+(2*(Pb*Pi)))))
```

```
## [1] 0.3905227
```

```r
(Pbi <- B*(2*(Pb*Pi))/((Pb^2)+(2*(Pb*Pi))))
```

```
## [1] 6.609477
```
These are re-estimates of allele frequencies and genotype estimations. 

```r
AB <- 5
A <- 30
B <- 7
O <- 36
(N <- AB + A + B + O)
```

```
## [1] 78
```

```r
(Pi <- sqrt(O/N))
```

```
## [1] 0.6793662
```

```r
(Pa <- sqrt ((A+O)/N)-Pi)
```

```
## [1] 0.2405
```

```r
(Pb <- sqrt ((B+O)/N)-Pi)
```

```
## [1] 0.06311748
```

```r
SQUARE <- function(q)q^2
SQUARE (5)
```

```
## [1] 25
```

```r
Pi <- sqrt(O/N)
Pa <- sqrt ((A+O)/N)-Pi
Pb <- sqrt ((B+O)/N)-Pi
Pi0 <- 0
Pa0 <- 0
Pb0 <- 0 
counter <- 0
EM <- function (Pi, Pa, Pb) {
  while((round(Pi0, 12) == round(Pi, 12)) == FALSE &&
         (round(Pa0, 12) == round(Pa, 12)) == FALSE &&
         (round(Pb0, 12) == round(Pb, 12)) == FALSE){
    Pi0 <- Pi
    Pa0 <- Pa
    Pb0 <- Pb
    Paa <- A*(Pa0^2/((Pa0^2)+2*(Pa0*Pi0)))
    Pai <- A * (2*(Pa0*Pi0))/((Pa0^2)+(2*(Pa0*Pi0)))
    Pbb <-B*(Pb0^2/((Pb0^2)+(2*(Pb0*Pi0))))
    Pbi <- B*(2*(Pb0*Pi0))/((Pb0^2)+(2*(Pb0*Pi0)))
    Pii <- O
    Pab <- AB
    (Pa <- ((2*Paa)+Pai+Pab)/(2*N))
    (Pb <- ((2*Pbb)+Pbi+Pab)/(2*N))
    (Pi <- ((2*Pii)+Pai+Pbi)/(2*N))
    counter <- counter + 1
  }
  return(c(paste ("Pi =", Pi, ", Pa =", Pa, ", Pb = ", Pb,
                   ", Number of Loops = ", counter)))
}
c(Pi, Pa, Pb)
```

```
## [1] 0.67936622 0.24049999 0.06311748
```

```r
EM(Pi, Pa, Pb)
```

```
## [1] "Pi = 0.665226278014308 , Pa = 0.255322004475303 , Pb =  0.0794517175103893 , Number of Loops =  12"
```
file:///C:/Users/sierr/Documents/Fall%2021/Project-1.html
