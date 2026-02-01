# Exercice 1

# ============================================================================
# Calcul de Pi avec la méthode d'Archimède
# ============================================================================
y <- 2 # valeur initiale
n <- 100 # nombre d'itérations
z <- 0
for (k in 1:n) {
  y_norm <- 2^(-k) * y
  y <- 2^k * sqrt(2 * (1 - sqrt(1 - y_norm^2)))
  z[k] <- y
  if (k == 20) cat("y_20 =", y, "\n")
  if (k == 100) cat("y_100 =", y, "\n")
}
plot(1:n,z, col="blue", pch=3,main="Evolution de y_k en fonction de n")


# ============================================================================
# Calcul de Pi avec la méthode d'Archimède en utilisant une identité remarquable
# ============================================================================

y <- 2 # valeur initiale
y_vals <- numeric(100)
for (k in 1:100) {
  y_norm <- 2^(-k) * y
  X <- sqrt(1 - y_norm^2)
  correction <- y_norm^2 / (1 + X)
  y <- 2^k * sqrt(2 * correction)
  y_vals[k] <- y
}

plot(1:100, y_vals, type = "l", col = "blue", lwd = 2,
     xlab = "k", ylab = "y_k (corrigé)", main = "y_k stable avec identité remarquable")
abline(h = pi, col = "red", lty = 2)
grid()