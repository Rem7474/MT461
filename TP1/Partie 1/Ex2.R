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

# ============================================================================
# ANALYSE THÉORIQUE: Borne supérieure de l'erreur
# ============================================================================
cat("\n=== ANALYSE THÉORIQUE ===\n\n")

cat("Encadrement: 0 ≤ I_n ≤ ∫₀¹ (x^n / a) dx = 1/(a(n+1))\n\n")

cat("Borne supérieure pour l'erreur théorique |d_n^(m)|:\n")
cat("où d_n^(m) est l'erreur accumulée à l'itération n\n\n")

# Calcul de la borne supérieure
upper_bound <- 1 / (a * (nvals + 1))

cat("Valeurs de I_n (Forward et Backward) et borne supérieure:\n")
cat(sprintf("%-4s %-15s %-15s %-15s\n", "n", "Forward", "Backward", "Borne sup"))
cat(strrep("-", 55), "\n")
for (i in 1:length(nvals)) {
  cat(sprintf("%-4d %-15.8f %-15.8f %-15.8f\n", 
              nvals[i], I_f[i], I_b[i], upper_bound[i]))
}

cat("\nObservations:\n")
cat("- Forward accumule les erreurs d'arrondi et diverge\n")
cat("- Backward reste stable dans les bornes théoriques\n")
cat("- L'erreur Forward-Backward augmente avec n\n\n")

# ============================================================================
# STABILITÉ NUMÉRIQUE DE L'ALGORITHME BACKWARD
# ============================================================================
cat("\n=== STABILITÉ NUMÉRIQUE - ALGORITHME BACKWARD ===\n\n")

cat("L'algorithme backward est numériquement stable pour a > 1:\n")
cat("Preuve: toute erreur ε à l'itération n est atténuée par le facteur 1/a\n")
cat("à l'itération n-1 (en remontant).\n\n")

cat("Chaque erreur est divisée par a à chaque remontée:\n")
cat("d_n-1 = (d_n - 0) / a  =>  erreur divisée par a à chaque pas\n\n")

cat(sprintf("Pour a = %d: facteur d'atténuation = 1/%d = %.4f\n", a, a, 1/a))
cat("=> Les erreurs sont exponentiellement atténuées\n\n")

# ============================================================================
# CALCUL DE I_20 À LA PRÉCISION MACHINE
# ============================================================================
cat("\n=== CALCUL DE I_20 À LA PRÉCISION MACHINE ===\n\n")

N_20 <- 20

# Forward (instable)
I_f_20 <- numeric(N_20 + 1)
I_f_20[1] <- log((1 + a) / a)
for (n in 2:(N_20 + 1)) {
  I_f_20[n] <- (1 / (n - 1)) + a * I_f_20[n - 1]
}

# Backward (stable)
I_b_20 <- numeric(N_20 + 1)
I_b_20[N_20 + 1] <- 0  # approximation pour I_N
for (n in N_20:1) {
  I_b_20[n] <- (I_b_20[n + 1] - (1 / n)) / a
}

cat(sprintf("Résultats pour N = %d:\n\n", N_20))
cat(sprintf("I_20 (Forward, instable):  %.16e\n", I_f_20[N_20 + 1]))
cat(sprintf("I_20 (Backward, stable):   %.16e\n", I_b_20[N_20 + 1]))
cat(sprintf("Borne supérieure théorique: %.16e\n", 1 / (a * (N_20 + 1))))
cat(sprintf("\nDifférence |Forward - Backward|: %.16e\n", abs(I_f_20[N_20 + 1] - I_b_20[N_20 + 1])))

cat("\nLe résultat du backward (stable) est fiable car dans la borne théorique.\n")
cat("Le résultat du forward (instable) a accumulé trop d'erreurs.\n\n")

# Graphique 3: Évolution comparée pour N=20
png(file.path(output_dir, "Ex2_3_Stabilite_N20.png"), width = 900, height = 600)
n_vals_20 <- 0:N_20
plot(n_vals_20, I_f_20, type = "b", col = "red", pch = 1, lwd = 2,
     main = "Stabilité numérique: Forward vs Backward (N=20)",
     xlab = "n", ylab = "I_n", ylim = range(c(I_f_20, I_b_20)))
lines(n_vals_20, I_b_20, type = "b", col = "blue", pch = 2, lwd = 2)
lines(n_vals_20, 1 / (a * (n_vals_20 + 1)), type = "l", col = "green", 
      lwd = 2, lty = 2)
legend("topright", 
       c("Forward (instable)", "Backward (stable)", "Borne sup théorique"),
       col = c("red", "blue", "green"), pch = c(1, 2, NA), lwd = 2, lty = c(1, 1, 2))
grid()
dev.off()

cat("Graphique sauvegardé: Ex2_3_Stabilite_N20.png\n")
