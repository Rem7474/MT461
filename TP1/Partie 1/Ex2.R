# Exercice 2

# ============================================================================
# CONFIGURATION: Dossier de sortie des images
# ============================================================================
output_dir <- "TP1/Partie 1/output"
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Paramètres
a <- 10
N <- 10

# ============================================================================
# Calcul de I_n = ∫_0^1 x^n e^{ax} dx par la méthode forward
# ============================================================================

I_f <- numeric(N+1)
I_f[1] <- log((1 + a)/a) # I0
for (n in 2:(N+1)) {
  I_f[n] <- (1/(n-1)) + a * I_f[n-1]
}

# ============================================================================
# Calcul de I_n = ∫_0^1 x^n e^{ax} dx par la méthode backward
# ============================================================================
I_b <- numeric(N+1)
I_b[N+1] <- 0 # approximation pour I_N
for (n in N:1) {
  I_b[n] <- (I_b[n+1] - (1/n)) / a
}

nvals <- 0:N

# ============================================================================
# Graphique 1: Superposition Forward vs Backward
# ============================================================================
png(file.path(output_dir, "Ex2_1_Forward_vs_Backward.png"), width = 900, height = 600)
plot(nvals, I_f, type="b", col="red", pch=1, ylim=range(c(I_f, I_b)), 
     main="Forward (rouge) vs Backward (bleu)", xlab="n", ylab="I_n", lwd = 2)
lines(nvals, I_b, type="b", col="blue", pch=2, lwd = 2)
legend("topright", legend=c("Forward", "Backward"), col=c("red", "blue"), pch=c(1,2), lwd = 2)
grid()
dev.off()

# ============================================================================
# Graphique 2: Erreur absolue
# ============================================================================
err <- abs(I_f - I_b)
png(file.path(output_dir, "Ex2_2_Erreur_absolue.png"), width = 900, height = 600)
plot(nvals, err, type="b", col="purple", pch=3, lwd = 2,
     main="Erreur absolue |Forward - Backward|", xlab="n", ylab="Erreur absolue")
grid()
dev.off()

cat("Les graphiques ont été sauvegardés dans:", file.path(getwd(), output_dir), "\n")
