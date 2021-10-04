p <- 0.5 #setting allele frequency
p^2 #frequency of homozygotes
f <- 1/1000 #new mutant allele
sqrt (f) #frequency of p-squared homozygotes, so this gives us frequency p
1-(sqrt(f)) #probability of not inheriting allele
curve (x^2, 0, 1) #relation between x and x^2 between zero and one
curve (x^2, 0, 1, 
       xlab = "Allele Frequencies", ylab = "Genotype Frequencies", 
       col = "green", lwd = 2)
text (0.6, 0.2, "Homozygotes", col = "green")
curve (2*x*(1-x), 0, 1, add = TRUE, 
       xlab = "Allele Frequencies", ylab = "Genotype Frequencies", 
       col = "blue", lwd = 2) #add TRUE allows us to plot over what we've already done
text (0.25, 0.5, "Heterozygotes", col = "blue")

p <- 2/3
p^2==2*p*(1-p) #produces rounding error, causing false
all.equal(p^2, 2*p*(1-p)) #fixes the rounding error
points (p, p^2, lwd=2, cex=2)

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

source ("HW1.R")

allele <- c("A", "A", "a", "a", "a", "a", "a", "a", "a", "a")
allele <- c(rep ("A", 2), rep ("a", 8))
print (allele)
popsize <- 100
pop <- matrix (nrow = popsize, ncol = 2) #creates a matrix of 200 rows and 2 columns
for (i in 1: popsize) {
  pop [i, 1] <- sample (allele, 1)
  pop [i, 2] <- sample (allele, 1)
} #for loop uses variable i to keep track of position from one to popsize
pop #shows heterozygotes/homozygotes for each individual
Acount <-0 #define variable to count # As that are found
for (i in 1: popsize) {
  if (pop [i, 1] == "A") {
    Acount <- Acount +1
  }
  if (pop [i, 2] == "A") {
    Acount <- Acount +1
  }
} #loop is asking if allele in each position (1 or 2) is equal to A
AFreq <- Acount/(popsize*2)
AFreq #count of A alleles by total number of allele copies

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
print (C(AFreq, HetFreq, AAFreq))
as.factor (AFreq) 
as.factor (HetFreq)
as.factor (AAFreq)
print (C(AFreq, HetFreq, AAFreq))
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
#simulated sampling alleles into three genotypes and comparing HW predictions

#Running 10 replicates and plotting sampling alleles
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

#allele randomly determined based on p allele frequency w/random draw
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

#allele frequencies from data
p <- (9+842/2)/23369
p
2*p*(1-p)
curve (2/x, 1e-7, 0.01, log = "y")
library ("popgenr")
data ("snp")
class (snp)
head (snp)
str (snp)
plot (snp$type)
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
