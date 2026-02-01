# Exercice 2 - Méthode des Trapèzes
# Stabilité numérique et équations raides

# ============================================================================
# CONFIGURATION: Dossier de sortie des images
# ============================================================================
output_dir <- "TP2/Partie 2/output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("Les images seront sauvegardées dans:", file.path(getwd(), output_dir), "\n\n")

# ============================================================================
# Problème de test: y'(x) = λy(x), y(0) = 1
# Solution exacte: y(x) = exp(λx)
# ============================================================================

lambda <- -0.5

# Solution exacte
y_exacte <- function(x, lambda) {
  exp(lambda * x)
}

# ============================================================================
# MÉTHODE 1: Euler explicite
# ============================================================================
euler_explicite <- function(lambda, x0, xf, y0, h) {
  n <- floor((xf - x0) / h)
  x <- seq(x0, x0 + n * h, by = h)
  y <- numeric(n + 1)
  y[1] <- y0
  
  for (i in 1:n) {
    y[i+1] <- y[i] * (1 + h * lambda)
  }
  
  return(list(x = x, y = y))
}

# ============================================================================
# MÉTHODE 2: Euler implicite
# ============================================================================
euler_implicite <- function(lambda, x0, xf, y0, h) {
  n <- floor((xf - x0) / h)
  x <- seq(x0, x0 + n * h, by = h)
  y <- numeric(n + 1)
  y[1] <- y0
  
  # y_{i+1} = y_i + h * lambda * y_{i+1}
  # y_{i+1} * (1 - h * lambda) = y_i
  # y_{i+1} = y_i / (1 - h * lambda)
  
  for (i in 1:n) {
    y[i+1] <- y[i] / (1 - h * lambda)
  }
  
  return(list(x = x, y = y))
}

# ============================================================================
# MÉTHODE 3: Méthode des trapèzes (implicite)
# ============================================================================
trapeze <- function(lambda, x0, xf, y0, h) {
  n <- floor((xf - x0) / h)
  x <- seq(x0, x0 + n * h, by = h)
  y <- numeric(n + 1)
  y[1] <- y0
  
  # y_{i+1} = y_i + (h/2) * (lambda * y_i + lambda * y_{i+1})
  # y_{i+1} = y_i + (h/2) * lambda * y_i + (h/2) * lambda * y_{i+1}
  # y_{i+1} * (1 - (h/2) * lambda) = y_i * (1 + (h/2) * lambda)
  # y_{i+1} = y_i * (1 + (h/2) * lambda) / (1 - (h/2) * lambda)
  
  coeff <- (1 + (h/2) * lambda) / (1 - (h/2) * lambda)
  
  for (i in 1:n) {
    y[i+1] <- y[i] * coeff
  }
  
  return(list(x = x, y = y))
}

# ============================================================================
# QUESTION 1: Euler explicite avec h = 0.5 et h = 5
# ============================================================================
cat("\n=== QUESTION 1: Euler explicite - Stabilité avec lambda = -0.5 ===\n\n")

x0 <- 0
xf <- 50
y0 <- 1

# Cas 1: h = 0.5 (stable)
h1 <- 0.5
sol_euler_h05 <- euler_explicite(lambda, x0, xf, y0, h1)
y_exact_h05 <- y_exacte(sol_euler_h05$x, lambda)

cat(sprintf("h = %.1f: Critère de stabilité |1 + h*lambda| = |1 + %.1f*(-0.5)| = |%.2f| < 1?\n", 
            h1, h1, 1 + h1 * lambda))
cat(sprintf("         %s (STABLE)\n", ifelse(abs(1 + h1 * lambda) < 1, "OUI", "NON")))

# Cas 2: h = 5 (instable)
h2 <- 5
sol_euler_h5 <- euler_explicite(lambda, x0, xf, y0, h2)
y_exact_h5 <- y_exacte(sol_euler_h5$x, lambda)

cat(sprintf("h = %.1f: Critère de stabilité |1 + h*lambda| = |1 + %.1f*(-0.5)| = |%.2f| < 1?\n", 
            h2, h2, 1 + h2 * lambda))
cat(sprintf("         %s (INSTABLE)\n", ifelse(abs(1 + h2 * lambda) < 1, "OUI", "NON")))

# Graphique 1: Euler h = 0.5 vs solution exacte
png(file.path(output_dir, "Ex2_Euler_h05_stable.png"), width = 900, height = 600)
plot(sol_euler_h05$x, sol_euler_h05$y, type = "b", col = "red", 
     main = "Euler explicite: h = 0.5 (stable)",
     xlab = "x", ylab = "y", lwd = 2, pch = 16, cex = 0.8)
lines(sol_euler_h05$x, y_exact_h05, col = "blue", lwd = 2.5, lty = 2)
legend("topright", c("Euler (h=0.5)", "Exacte (exp(-0.5x))"), 
       col = c("red", "blue"), lwd = 2, lty = c(1, 2))
grid()
dev.off()

# Graphique 2: Euler h = 5 vs solution exacte (montre l'instabilité)
png(file.path(output_dir, "Ex2_Euler_h5_instable.png"), width = 900, height = 600)
# Limitation pour la visualisation
ylim_val <- min(max(sol_euler_h5$y[1:50], na.rm = TRUE), 1e6)
plot(sol_euler_h5$x, sol_euler_h5$y, type = "b", col = "red", 
     main = "Euler explicite: h = 5 (instable)",
     xlab = "x", ylab = "y", lwd = 2, pch = 16, cex = 0.8,
     ylim = c(-ylim_val, ylim_val))
lines(sol_euler_h5$x, y_exact_h5, col = "blue", lwd = 2.5, lty = 2)
legend("topright", c("Euler (h=5) - INSTABLE", "Exacte (exp(-0.5x))"), 
       col = c("red", "blue"), lwd = 2, lty = c(1, 2))
grid()
dev.off()

cat(sprintf("\nObservations:\n"))
cat(sprintf("- Avec h=0.5: La solution approchée suit l'exacte (facteur %.4f < 1)\n", 1 + h1 * lambda))
cat(sprintf("- Avec h=5: La solution oscille et diverge (facteur %.2f, amplitude croît)\n", 1 + h2 * lambda))

# ============================================================================
# QUESTION 2: Analyse de la stabilité
# ============================================================================
cat("\n=== QUESTION 2: Condition de stabilité pour Euler explicite ===\n\n")

cat("Pour Euler explicite: y_{i+1} = y_i * (1 + h*lambda)\n")
cat("Le facteur d'amplification est: ρ(h,lambda) = 1 + h*lambda\n\n")
cat("Stabilité requiert: |ρ(h,lambda)| < 1\n")
cat("Soit: |1 + h*lambda| < 1\n\n")

# Pour lambda < 0 (équation décroissante):
cat("Pour lambda < 0 (décroissance):\n")
cat("  - |1 + h*lambda| < 1  ⟺  -1 < 1 + h*lambda < 1\n")
cat("  - Borne inférieure: 1 + h*lambda > -1  ⟹  h*lambda > -2  ⟹  h < -2/lambda = 2/|lambda|\n")
cat(sprintf("  - Ici lambda = %.1f, donc h_max = 2/|lambda| = %.1f\n", lambda, 2/abs(lambda)))

# Graphique 3: Région de stabilité
png(file.path(output_dir, "Ex2_Stabilite_Euler.png"), width = 900, height = 600)
h_vals <- seq(0.01, 10, by = 0.01)
rho <- 1 + h_vals * lambda
plot(h_vals, rho, type = "l", col = "blue", lwd = 2.5,
     main = "Région de stabilité - Euler explicite (λ = -0.5)",
     xlab = "Pas h", ylab = "Facteur d'amplification |1 + h*λ|", 
     ylim = c(-2, 2))
abline(h = 1, col = "green", lwd = 2, lty = 2)
abline(h = -1, col = "green", lwd = 2, lty = 2)
abline(v = 4, col = "red", lwd = 2, lty = 2)
abline(v = h1, col = "orange", lwd = 2.5, lty = 1)
abline(v = h2, col = "purple", lwd = 2.5, lty = 1)

# Zone stable
polygon(c(0, 4, 4, 0), c(-1, -1, 1, 1), col = rgb(0, 1, 0, 0.2), border = NA)
text(2, 0.7, "Zone stable\n|ρ| < 1", cex = 1.2, font = 2)

legend("topright", c("ρ(h)", "Limite", "h_crit = 4", "h=0.5 ✓", "h=5 ✗"), 
       col = c("blue", "green", "red", "orange", "purple"), lwd = 2, lty = c(1, 2, 2, 1, 1))
grid()
dev.off()

# ============================================================================
# QUESTION 3: Euler implicite vs Euler explicite
# ============================================================================
cat("\n=== QUESTION 3: Euler implicite (inconditionnellement stable) ===\n\n")

# Même test avec Euler implicite
sol_implicit_h05 <- euler_implicite(lambda, x0, xf, y0, h1)
sol_implicit_h5 <- euler_implicite(lambda, x0, xf, y0, h2)

cat(sprintf("Euler implicite: y_{i+1} = y_i / (1 - h*lambda)\n\n"))
cat(sprintf("h = %.1f: Facteur d'amplification = 1/(1 - %.1f*(-0.5)) = 1/%.2f = %.4f < 1 ✓ STABLE\n", 
            h1, h1, 1 - h1 * lambda, 1/(1 - h1 * lambda)))
cat(sprintf("h = %.1f: Facteur d'amplification = 1/(1 - %.1f*(-0.5)) = 1/%.2f = %.4f < 1 ✓ STABLE\n", 
            h2, h2, 1 - h2 * lambda, 1/(1 - h2 * lambda)))

# Graphique 4: Comparaison Euler explicite vs implicite (h=5)
png(file.path(output_dir, "Ex2_Euler_Explicit_vs_Implicit_h5.png"), width = 900, height = 600)
plot(sol_euler_h5$x, sol_euler_h5$y, type = "l", col = "red", 
     main = "Comparaison: Euler explicite vs implicite (h = 5, instable/stable)",
     xlab = "x", ylab = "y", lwd = 2, lty = 1,
     ylim = c(-1, 2))
lines(sol_implicit_h5$x, sol_implicit_h5$y, col = "green", lwd = 2.5, lty = 2)
lines(sol_euler_h5$x, y_exact_h5, col = "blue", lwd = 2.5, lty = 3)
legend("bottomright", c("Euler explicite (instable)", "Euler implicite (stable)", "Exacte"), 
       col = c("red", "green", "blue"), lwd = 2, lty = c(1, 2, 3))
grid()
dev.off()

# ============================================================================
# QUESTION 4: Méthode des trapèzes
# ============================================================================
cat("\n=== QUESTION 4: Méthode des trapèzes (ordre 2, inconditionnellement stable) ===\n\n")

sol_trap_h05 <- trapeze(lambda, x0, xf, y0, h1)
sol_trap_h5 <- trapeze(lambda, x0, xf, y0, h2)

cat(sprintf("Méthode des trapèzes: y_{i+1} = y_i * (1 + h*lambda/2) / (1 - h*lambda/2)\n\n"))
cat(sprintf("h = %.1f: Facteur = (1 + %.2f)/(1 - %.2f) = %.4f\n", 
            h1, h1 * lambda / 2, h1 * lambda / 2, (1 + h1 * lambda / 2)/(1 - h1 * lambda / 2)))
cat(sprintf("h = %.1f: Facteur = (1 + %.2f)/(1 - %.2f) = %.4f\n", 
            h2, h2 * lambda / 2, h2 * lambda / 2, (1 + h2 * lambda / 2)/(1 - h2 * lambda / 2)))

# Graphique 5: Comparaison des trois méthodes (h=5)
png(file.path(output_dir, "Ex2_Comparaison_methodes_h5.png"), width = 900, height = 600)
plot(sol_euler_h5$x, sol_euler_h5$y, type = "l", col = "red", 
     main = "Comparaison des trois méthodes (h = 5)",
     xlab = "x", ylab = "y", lwd = 2, lty = 1,
     ylim = c(-1, 2))
lines(sol_implicit_h5$x, sol_implicit_h5$y, col = "green", lwd = 2, lty = 2)
lines(sol_trap_h5$x, sol_trap_h5$y, col = "purple", lwd = 2, lty = 3)
lines(sol_euler_h5$x, y_exact_h5, col = "blue", lwd = 2.5, lty = 4)
legend("bottomright", c("Euler explicite", "Euler implicite", "Trapèzes", "Exacte"), 
       col = c("red", "green", "purple", "blue"), lwd = 2, lty = c(1, 2, 3, 4))
grid()
dev.off()

# Graphique 6: Comparaison (h=0.5) - toutes les méthodes convergent
png(file.path(output_dir, "Ex2_Comparaison_methodes_h05.png"), width = 900, height = 600)
sol_euler_h05 <- euler_explicite(lambda, x0, xf, y0, h1)
sol_implicit_h05 <- euler_implicite(lambda, x0, xf, y0, h1)
sol_trap_h05 <- trapeze(lambda, x0, xf, y0, h1)

plot(sol_euler_h05$x, sol_euler_h05$y, type = "l", col = "red", 
     main = "Comparaison des trois méthodes (h = 0.5)",
     xlab = "x", ylab = "y", lwd = 2, lty = 1)
lines(sol_implicit_h05$x, sol_implicit_h05$y, col = "green", lwd = 2, lty = 2)
lines(sol_trap_h05$x, sol_trap_h05$y, col = "purple", lwd = 2, lty = 3)
lines(sol_euler_h05$x, y_exact_h05, col = "blue", lwd = 2.5, lty = 4)
legend("topright", c("Euler explicite", "Euler implicite", "Trapèzes", "Exacte"), 
       col = c("red", "green", "purple", "blue"), lwd = 2, lty = c(1, 2, 3, 4))
grid()
dev.off()

# ============================================================================
# ANALYSE DE CONVERGENCE - Ordre de la méthode des trapèzes
# ============================================================================
cat("\n=== ANALYSE DE CONVERGENCE - Ordre des méthodes ===\n\n")

h_values <- 10^seq(-3, -0.5, by = 0.1)  # h de 10^-3 à ~0.3
errors_euler <- numeric(length(h_values))
errors_implicit <- numeric(length(h_values))
errors_trap <- numeric(length(h_values))

# Solution de référence avec un très petit h
h_ref <- 0.001
sol_ref <- trapeze(lambda, x0, xf, y0, h_ref)
y_exact_ref <- y_exacte(sol_ref$x, lambda)

for (i in seq_along(h_values)) {
  h <- h_values[i]
  
  # Euler explicite
  sol <- euler_explicite(lambda, x0, xf, y0, h)
  # Interpoler pour comparaison
  y_approx <- rep(NA, length(sol_ref$x))
  indices <- match(round(sol_ref$x / h) * h, sol$x, nomatch = NA)
  valid_idx <- !is.na(indices)
  y_approx[valid_idx] <- sol$y[indices[valid_idx]]
  errors_euler[i] <- max(abs(y_approx[valid_idx] - y_exact_ref[valid_idx]), na.rm = TRUE)
  
  # Euler implicite
  sol <- euler_implicite(lambda, x0, xf, y0, h)
  y_approx <- rep(NA, length(sol_ref$x))
  indices <- match(round(sol_ref$x / h) * h, sol$x, nomatch = NA)
  valid_idx <- !is.na(indices)
  y_approx[valid_idx] <- sol$y[indices[valid_idx]]
  errors_implicit[i] <- max(abs(y_approx[valid_idx] - y_exact_ref[valid_idx]), na.rm = TRUE)
  
  # Trapèzes
  sol <- trapeze(lambda, x0, xf, y0, h)
  y_approx <- rep(NA, length(sol_ref$x))
  indices <- match(round(sol_ref$x / h) * h, sol$x, nomatch = NA)
  valid_idx <- !is.na(indices)
  y_approx[valid_idx] <- sol$y[indices[valid_idx]]
  errors_trap[i] <- max(abs(y_approx[valid_idx] - y_exact_ref[valid_idx]), na.rm = TRUE)
}

# Calcul des ordres
fit_euler <- lm(log(errors_euler) ~ log(h_values))
fit_implicit <- lm(log(errors_implicit) ~ log(h_values))
fit_trap <- lm(log(errors_trap) ~ log(h_values))

order_euler <- -fit_euler$coefficients[2]
order_implicit <- -fit_implicit$coefficients[2]
order_trap <- -fit_trap$coefficients[2]

cat(sprintf("Ordre d'Euler explicite:      %.3f (théorique: 1)\n", order_euler))
cat(sprintf("Ordre d'Euler implicite:      %.3f (théorique: 1)\n", order_implicit))
cat(sprintf("Ordre des trapèzes:           %.3f (théorique: 2)\n", order_trap))

# Graphique 7: Convergence
png(file.path(output_dir, "Ex2_Ordres_convergence.png"), width = 900, height = 600)
plot(h_values, errors_euler, type = "b", pch = 1, log = "xy",
     main = "Ordre de convergence des trois méthodes",
     xlab = "Pas h", ylab = "Erreur globale max",
     lwd = 2, cex = 1.2, col = "red", ylim = c(1e-8, 1e-1))
points(h_values, errors_implicit, type = "b", col = "green", lwd = 2, cex = 1.2, pch = 2)
points(h_values, errors_trap, type = "b", col = "purple", lwd = 2, cex = 1.2, pch = 3)

# Droites de tendance
x_line <- range(h_values)
y_euler <- exp(fit_euler$coefficients[1]) * x_line^(-order_euler)
y_implicit <- exp(fit_implicit$coefficients[1]) * x_line^(-order_implicit)
y_trap <- exp(fit_trap$coefficients[1]) * x_line^(-order_trap)
lines(x_line, y_euler, col = "red", lwd = 2, lty = 2)
lines(x_line, y_implicit, col = "green", lwd = 2, lty = 2)
lines(x_line, y_trap, col = "purple", lwd = 2, lty = 2)

legend("bottomright",
       c(sprintf("Euler explicite (ordre %.2f)", order_euler),
         sprintf("Euler implicite (ordre %.2f)", order_implicit),
         sprintf("Trapèzes (ordre %.2f)", order_trap)),
       col = c("red", "green", "purple"), lwd = 2, pch = c(1, 2, 3))
grid()
dev.off()

# ============================================================================
# RÉSUMÉ ET CONCLUSIONS
# ============================================================================
cat("\n=== RÉSUMÉ ET CONCLUSIONS ===\n\n")

cat("1. STABILITÉ NUMÉRIQUE:\n")
cat("   - Euler explicite: conditionnellement stable (h < 4 pour λ = -0.5)\n")
cat("   - Euler implicite: inconditionnellement stable (stable pour tout h)\n")
cat("   - Trapèzes: inconditionnellement stable (stable pour tout h)\n\n")

cat("2. PRÉCISION:\n")
cat("   - Euler explicite et implicite: ordre 1 (O(h))\n")
cat("   - Trapèzes: ordre 2 (O(h²))\n\n")

cat("3. RECOMMANDATIONS:\n")
cat("   - Pour les équations raides (λ très négatif): préférer implicite ou trapèzes\n")
cat("   - Les trapèzes offrent le meilleur compromis (ordre 2 + stabilité)\n\n")

cat("\n Tous les graphiques ont été générés.\n")
