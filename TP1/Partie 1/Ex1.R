# Exercice 1

# ============================================================================
# CONFIGURATION: Dossier de sortie des images
# ============================================================================
output_dir <- "TP1/Partie 1/output"
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

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

png(file.path(output_dir, "Ex1_1_Evolution_y_k.png"), width = 900, height = 600)
plot(1:n, z, col="blue", pch=3, lwd = 2,
     main="Evolution de y_k en fonction de n",
     xlab="k", ylab="y_k")
grid()
dev.off()


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

png(file.path(output_dir, "Ex1_2_Evolution_identite_remarquable.png"), width = 900, height = 600)
plot(1:100, y_vals, type = "l", col = "blue", lwd = 2,
     xlab = "k", ylab = "y_k (corrigé)", main = "y_k stable avec identité remarquable")
abline(h = pi, col = "red", lty = 2, lwd = 2)
legend("bottomright", c("y_k calculé", "π exact"), col = c("blue", "red"), lwd = 2, lty = c(1, 2))
grid()
dev.off()

cat("Les graphiques ont été sauvegardés dans:", file.path(getwd(), output_dir), "\n")