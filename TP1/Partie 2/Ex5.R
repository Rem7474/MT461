# Exercice 5 - Calcul de racines q-ièmes
# Résolution de F(x) = x^q - α = 0

# ============================================================================
# CONFIGURATION: Dossier de sortie des images
# ============================================================================
output_dir <- "TP1/Partie 2/output"
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# Question 1: Méthode de Newton comme point fixe
# ============================================================================
cat("\n=== QUESTION 1: Itération de Newton pour x^q - α = 0 ===\n\n")

cat("F(x) = x^q - α\n")
cat("F'(x) = q*x^(q-1)\n\n")

cat("Méthode de Newton: x_{n+1} = x_n - F(x_n)/F'(x_n)\n")
cat("                   x_{n+1} = x_n - (x_n^q - α)/(q*x_n^(q-1))\n")
cat("                   x_{n+1} = x_n - x_n^q/(q*x_n^(q-1)) + α/(q*x_n^(q-1))\n")
cat("                   x_{n+1} = x_n - x_n/q + α/(q*x_n^(q-1))\n")
cat("                   x_{n+1} = ((q-1)*x_n + α/x_n^(q-1))/q\n\n")

cat("Donc f(x) = ((q-1)*x + α/x^(q-1))/q\n\n")

# Fonction de point fixe
fixed_point_func <- function(x, alpha, q) {
  ((q - 1) * x + alpha / x^(q - 1)) / q
}

# Dérivée de la fonction de point fixe
fixed_point_derivative <- function(x, alpha, q) {
  (q - 1) / q - (q - 1) * alpha / (q * x^q)
}

cat("Convergence: f'(x*) où x* = α^(1/q)\n")
cat("f'(x) = (q-1)/q - (q-1)*α/(q*x^q)\n")
cat("f'(x*) = (q-1)/q - (q-1)*α/(q*α) = (q-1)/q - (q-1)/q = 0\n\n")
cat("Donc |f'(x*)| = 0 < 1 => convergence quadratique (Newton)\n\n")

# ============================================================================
# Question 2: Calcul de ³√5 par Newton
# ============================================================================
cat("\n=== QUESTION 2: Calcul de ³√5 par Newton ===\n\n")

# Méthode de Newton pour racine q-ième
newton_root <- function(alpha, q, x0, max_iter = 20, tol = 1e-16) {
  x <- numeric(max_iter + 1)
  x[1] <- x0
  
  for (i in 1:max_iter) {
    x[i + 1] <- fixed_point_func(x[i], alpha, q)
    
    # Critère d'arrêt
    if (abs(x[i + 1] - x[i]) < tol) {
      x <- x[1:(i + 1)]
      break
    }
  }
  
  return(x)
}

# Calcul de ³√5
alpha <- 5
q <- 3
x0 <- 2  # Condition initiale (proche de la solution)
x_true <- alpha^(1/q)

result_newton <- newton_root(alpha, q, x0)

# Calcul des erreurs
errors_newton <- abs(result_newton - x_true)

# Nombre de décimales significatives: -log10(erreur)
decimals_newton <- -log10(errors_newton)
decimals_newton[is.infinite(decimals_newton)] <- 16  # Précision machine

cat(sprintf("Valeur exacte: %.16f\n", x_true))
cat("\nÉvolution de Newton:\n")
cat(sprintf("%-4s %-20s %-15s %-15s\n", "Iter", "x_n", "Erreur", "Décimales"))
cat(strrep("-", 60), "\n")

for (i in seq_along(result_newton)) {
  cat(sprintf("%-4d %-20.16f %-15.2e %-15.1f\n", 
              i-1, result_newton[i], errors_newton[i], decimals_newton[i]))
}

cat("\nObservation: Le nombre de décimales significatives DOUBLE à chaque itération\n")
cat("=> Convergence QUADRATIQUE de la méthode de Newton\n\n")

# Graphique 1: Convergence de Newton
png(file.path(output_dir, "Ex5_2_Newton_convergence.png"), width = 900, height = 600)
par(mfrow = c(2, 1))

# Erreur en échelle log
plot(0:(length(errors_newton)-1), log10(errors_newton), type = "b", pch = 16, 
     col = "blue", lwd = 2,
     main = "Convergence de Newton pour ³√5",
     xlab = "Itération", ylab = "log10(Erreur)")
grid()

# Nombre de décimales significatives
plot(0:(length(decimals_newton)-1), decimals_newton, type = "b", pch = 16, 
     col = "darkgreen", lwd = 2,
     main = "Décimales significatives - Newton",
     xlab = "Itération", ylab = "Nombre de décimales")
grid()

dev.off()

# ============================================================================
# Question 3: Méthode de la sécante pour ³√5
# ============================================================================
cat("\n=== QUESTION 3: Méthode de la sécante pour ³√5 ===\n\n")

# Méthode de la sécante
secant_root <- function(alpha, q, x0, x1, max_iter = 20, tol = 1e-16) {
  x <- numeric(max_iter + 2)
  x[1] <- x0
  x[2] <- x1
  
  F <- function(x) x^q - alpha
  
  for (i in 2:(max_iter + 1)) {
    f_i <- F(x[i])
    f_i_minus_1 <- F(x[i - 1])
    
    # Formule de la sécante
    if (abs(f_i - f_i_minus_1) < 1e-16) break
    
    x[i + 1] <- x[i] - f_i * (x[i] - x[i - 1]) / (f_i - f_i_minus_1)
    
    # Critère d'arrêt
    if (abs(x[i + 1] - x[i]) < tol) {
      x <- x[1:(i + 1)]
      break
    }
  }
  
  return(x)
}

# Calcul avec la sécante
x0_sec <- 2
x1_sec <- 1.5
result_secant <- secant_root(alpha, q, x0_sec, x1_sec)

# Calcul des erreurs
errors_secant <- abs(result_secant - x_true)
decimals_secant <- -log10(errors_secant)
decimals_secant[is.infinite(decimals_secant)] <- 16

cat("\nÉvolution de la sécante:\n")
cat(sprintf("%-4s %-20s %-15s %-15s\n", "Iter", "x_n", "Erreur", "Décimales"))
cat(strrep("-", 60), "\n")

for (i in seq_along(result_secant)) {
  cat(sprintf("%-4d %-20.16f %-15.2e %-15.1f\n", 
              i-1, result_secant[i], errors_secant[i], decimals_secant[i]))
}

cat("\nObservation: La sécante converge plus lentement que Newton\n")
cat("Ordre de convergence ≈ φ = (1+√5)/2 ≈ 1.618 (nombre d'or)\n\n")

# Graphique 2: Comparaison Newton vs Sécante
png(file.path(output_dir, "Ex5_3_Newton_vs_Secant.png"), width = 900, height = 600)
par(mfrow = c(2, 1))

# Erreur en échelle log
max_len <- max(length(errors_newton), length(errors_secant))
plot(0:(length(errors_newton)-1), log10(errors_newton), type = "b", pch = 16, 
     col = "blue", lwd = 2, xlim = c(0, max_len-1),
     main = "Comparaison Newton vs Sécante pour ³√5",
     xlab = "Itération", ylab = "log10(Erreur)")
lines(0:(length(errors_secant)-1), log10(errors_secant), type = "b", pch = 17, 
      col = "red", lwd = 2)
legend("topright", c("Newton (ordre 2)", "Sécante (ordre ≈1.618)"), 
       col = c("blue", "red"), pch = c(16, 17), lwd = 2)
grid()

# Nombre de décimales significatives
plot(0:(length(decimals_newton)-1), decimals_newton, type = "b", pch = 16, 
     col = "blue", lwd = 2, xlim = c(0, max_len-1),
     main = "Décimales significatives - Comparaison",
     xlab = "Itération", ylab = "Nombre de décimales")
lines(0:(length(decimals_secant)-1), decimals_secant, type = "b", pch = 17, 
      col = "red", lwd = 2)
legend("topleft", c("Newton", "Sécante"), 
       col = c("blue", "red"), pch = c(16, 17), lwd = 2)
grid()

dev.off()

# ============================================================================
# Analyse de convergence théorique
# ============================================================================
cat("\n=== ANALYSE THÉORIQUE ===\n\n")

# Vérifier l'ordre de convergence empiriquement
cat("Ordre de convergence empirique:\n\n")

# Newton
if (length(errors_newton) >= 4) {
  ratios_newton <- numeric(length(errors_newton) - 3)
  for (i in 3:(length(errors_newton) - 1)) {
    if (errors_newton[i-1] > 0 && errors_newton[i] > 0) {
      ratios_newton[i-2] <- log(errors_newton[i] / errors_newton[i-1]) / 
                            log(errors_newton[i-1] / errors_newton[i-2])
    }
  }
  order_newton <- mean(ratios_newton, na.rm = TRUE)
  cat(sprintf("Newton: ordre ≈ %.3f (théorique: 2)\n", order_newton))
}

# Sécante
if (length(errors_secant) >= 4) {
  ratios_secant <- numeric(length(errors_secant) - 3)
  for (i in 3:(length(errors_secant) - 1)) {
    if (errors_secant[i-1] > 0 && errors_secant[i] > 0) {
      ratios_secant[i-2] <- log(errors_secant[i] / errors_secant[i-1]) / 
                            log(errors_secant[i-1] / errors_secant[i-2])
    }
  }
  order_secant <- mean(ratios_secant, na.rm = TRUE)
  phi <- (1 + sqrt(5)) / 2
  cat(sprintf("Sécante: ordre ≈ %.3f (théorique: φ = %.3f)\n", order_secant, phi))
}

cat("\n=== RÉSUMÉ ===\n\n")
cat("1. Newton converge QUADRATIQUEMENT (ordre 2):\n")
cat("   - Nombre de décimales correctes DOUBLE à chaque itération\n")
cat("   - Nécessite le calcul de F'(x) à chaque itération\n\n")

cat("2. Sécante converge avec ordre φ ≈ 1.618:\n")
cat("   - Plus lente que Newton mais plus rapide que linéaire\n")
cat("   - Ne nécessite PAS le calcul de dérivée\n")
cat("   - Utilise une approximation de la dérivée par différences finies\n\n")

cat("3. Pour ³√5:\n")
cat(sprintf("   - Newton atteint la précision machine en %d itérations\n", 
            length(result_newton) - 1))
cat(sprintf("   - Sécante atteint la précision machine en %d itérations\n", 
            length(result_secant) - 1))

cat("\n Tous les graphiques ont été générés.\n")
