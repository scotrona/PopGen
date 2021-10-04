AA <- 366
AS <- 123
SS <- 12
n <- AA + AS + SS
print (paste ("n:", n))
p <- (SS + (AS/2))/n
print (paste ("p:", p))
EAA <- n*(1-p)^2
EAS <- n*2*p*(1-p)
ESS <- n*p^2
print (paste ("Expected:", EAA, EAS, ESS))
chi2 <- (EAA - AA)^2/EAA + 
  (EAS-AS)^2/EAS +
  (ESS-SS)^2/ESS
print (paste ("chi-square:", chi2))
pvalue <- pchisq(chi2, df = 1, lower.tail = FALSE)
print (paste ("P-value:", pvalue))
geno <- c (AA, AS, SS)
expe <- c (EAA, EAS, ESS)
G <- 2 * sum (geno * log (geno/expe))
print (paste("G:", G))
pvalue <- pchisq(G, df = 1, lower.tail = FALSE)
print (paste ("P-value:", pvalue))


AA <- 34
AS <- 100
SS <- 0
n <- AA + AS + SS
p <- (SS + (AS/2))/n
EAA <- n*(1-p)^2
EAS <- n*2*p*(1-p)
ESS <- n*p^2
geno <- c (AA, AS, SS)
expe <- c (EAA, EAS, ESS)
chi2 <- sum ((expe-geno)^2/expe)
print (paste ("chi square:", chi2))
pvalue <- pchisq(G, df = 1, lower.tail = FALSE)
print (paste ("P-value:", pvalue))
dat <- matrix (c(geno, expe), nrow = 2, byrow = T)
barplot (dat, beside = T, 
         col = c ("turquoise4", "sienna1"),
         names.arg = c ("AA", "SA", "SS"))
legend (x = "topright", legend = c ("Observed", "Expected"), 
        pch = 15, col = c ("turquoise4", "sienna1") )


p1 <- 0.2
p2 <- 0.3
p3 <- 0.5
sapply (c(p1, p2, p3), function (x) x^2)
sum (sapply (c(p1, p2, p3), function (x) x^2))

library ("popgenr")
data (genotypes)
str (genotypes)
rownames (genotypes) <- genotypes$ID
genotypes <- genotypes [, -c (1,2)]
(num.loci <- (length(genotypes))/2)
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
Het_exp <- 1 - Hom_exp
Het_obs <- 1 - Hom_obs
plot (Het_obs, Het_exp)
abline (lm (Het_exp ~ Het_obs))
reg <- summary (lm (Het_exp ~ Het_obs))
print (reg)
rr <- reg$r.squared
rrlabel <- paste ("r-squared =", round (rr, digits = 3))
pv <- reg$coefficients [2,4]
pvlabel <- paste ("P-value =", pv)
text (0.6, 0.2, rrlabel)
text (0.6, 0.15, pvlabel)


# Blood types and allele frequencies
AB <- 5
A <- 30
B <- 7
O <- 36
(N <- AB + A + B + O)
A/N
B/N
AB/N
O/N
(Pi <- sqrt(O/N))
(Pa <- sqrt ((A+O)/N)-Pi)
(Pb <- sqrt ((B+O)/N)-Pi)
signif (Pa + Pb + Pi, 1)
(Paa <- A*(Pa^2/((Pa^2)+2*(Pa*Pi))))
(Pai <- A * (2*(Pa*Pi))/((Pa^2)+(2*(Pa*Pi))))
(Pbb <-B*(Pb^2/((Pb^2)+(2*(Pb*Pi)))))
(Pbi <- B*(2*(Pb*Pi))/((Pb^2)+(2*(Pb*Pi))))
(Pii <- O)
(Pab <- AB)
#allele re-estimates
(Pa <- ((2*Paa)+Pai+Pab)/(2*N))
(Pb <- ((2*Pbb)+Pbi+Pab)/(2*N))
(Pi <- ((2*Pii)+Pai+Pbi)/(2*N))
(Paa <- A*(Pa^2/((Pa^2)+2*(Pa*Pi))))
(Pai <- A * (2*(Pa*Pi))/((Pa^2)+(2*(Pa*Pi))))
(Pbb <-B*(Pb^2/((Pb^2)+(2*(Pb*Pi)))))
(Pbi <- B*(2*(Pb*Pi))/((Pb^2)+(2*(Pb*Pi))))
(Paa <- A*(Pa^2/((Pa^2)+2*(Pa*Pi))))
(Pai <- A * (2*(Pa*Pi))/((Pa^2)+(2*(Pa*Pi))))
(Pbb <-B*(Pb^2/((Pb^2)+(2*(Pb*Pi)))))
(Pbi <- B*(2*(Pb*Pi))/((Pb^2)+(2*(Pb*Pi))))
rm(list = ls())
AB <- 5
A <- 30
B <- 7
O <- 36
(N <- AB + A + B + O)
(Pi <- sqrt(O/N))
(Pa <- sqrt ((A+O)/N)-Pi)
(Pb <- sqrt ((B+O)/N)-Pi)
SQUARE <- function(q)q^2
SQUARE (5)

Pi0 <- 0
Pa0 <- 0
Pb0 <- 0 
counter <- 0
EM <- function (Pi, Pa, Pb) {
  while ((round(Pi0, 12) == round (Pi, 12)) == FALSE &&
         (round(Pa0, 12) == round (Pa, 12)) == FALSE &&
         (round(Pb0, 12) == round (Pb, 12)) == FALSE &&) {
    Pi0 <- Pi
    Pa0 <- Pa
    Pb0 <- Pb
    Paa <- A*(Pa0^2/((Pa0^2)+2*(Pa0*Pi0)))
    Pai <- A * (2*(Pa0*Pi0))/((Pa0^2)+(2*(Pa0*Pi0)))
    Pbb <-B*(Pb0^2/((Pb0^2)+(2*(Pb0*Pi0))))
    Pbi <- B*(2*(P0b*Pi0))/((Pb0^2)+(2*(Pb0*Pi0)))
    Pii <- O
    Pab <- AB
    (Pa <- ((2*Paa)+Pai+Pab)/(2*N))
    (Pb <- ((2*Pbb)+Pbi+Pab)/(2*N))
    (Pi <- ((2*Pii)+Pai+Pbi)/(2*N))
    counter <- counter + 1
  }
  return (c(paste ("Pi =", Pi, ", Pa =", Pa, ", Pb = ", Pb,
                   ", Number of Loops = ", counter)))
}
c (Pi, Pa, Pb)
EM (Pi, Pa, Pb)
