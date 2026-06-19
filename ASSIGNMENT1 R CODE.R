
# --- Question 1: Chi-squared test and Cressie–Read ---

# 1) Create the 4×4 contingency table (Exposure × Grade)
table_data <- matrix(
  c(
     5,  8,  6,  7,
     9, 13, 12,  7,
    11, 20,  7,  6,
     4,  6,  6,  4
  ),
  nrow = 4,
  byrow = TRUE
)

rownames(table_data) <- c("Y1", "Y2", "Y3", "Y4")
colnames(table_data) <- c("G1", "G2", "G3", "G4")

# 2) Standard Pearson chi-squared test
res_chi <- chisq.test(table_data)
res_chi

# 3) Compare p-values using Monte Carlo (B=10, 100, 1000, 10000)
chisq.test(table_data, simulate.p.value = TRUE, B = 10)
chisq.test(table_data, simulate.p.value = TRUE, B = 100)
chisq.test(table_data, simulate.p.value = TRUE, B = 1000)
chisq.test(table_data, simulate.p.value = TRUE, B = 10000)

# 4) Manually compute the chi-squared statistic
O <- table_data
row_sums <- rowSums(O)
col_sums <- colSums(O)
total <- sum(O)
E <- outer(row_sums, col_sums, "*") / total

chi2_manual <- sum( (O - E)^2 / E )
chi2_manual

# 5) Cressie–Read statistics
cressie_read <- function(O, E, lambda){
  if(lambda == 0){
    mask <- (O > 0)
    return( 2 * sum( O[mask] * log(O[mask] / E[mask]) ) )
  } else {
    return((2 / (lambda*(lambda+1))) * sum( O * ( (O/E)^lambda - 1 ) ))
  }
}

for(lam in c(0, 0.5, 1)){
  val <- cressie_read(O, E, lam)
  cat("Lambda =", lam, " => Cressie–Read =", val, "\n")
}
