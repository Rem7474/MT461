# Exercice 4 — Newton-Raphson pour système 2D
# Résolution d'un système d'équations non-linéaires :
# f1(x1, x2) = x1² - x2² + 2x1 = 0
# f2(x1, x2) = x1² x2 + x2 - 1 = 0

# ============================================================================
# CONFIGURATION: Dossier de sortie des images
# ============================================================================
output_dir <- "TP1/Partie 2/output"
dir.create(output_dir, showWarnings = FALSE)
cat("Les images seront sauvegardées dans:", file.path(getwd(), output_dir), "\n\n")

# ============================================================================
# QUESTION i : Localisation graphique des solutions
# ============================================================================
cat("=== QUESTION i : Localisation graphique des solutions ===\n\n")

# Définition du système d'équations F(x) = 0
systeme_equations <- function(x) {
  x1 <- x[1]
  x2 <- x[2]
  eq1 <- x1^2 - x2^2 + 2*x1
  eq2 <- x1^2 * x2 + x2 - 1
  return(c(eq1, eq2))
}

# Jacobienne J(x) = [∂fi/∂xj]
jacobienne <- function(x) {
  x1 <- x[1]
  x2 <- x[2]
  
  # Calcul des dérivées partielles
  # f1 = x1² - x2² + 2x1
  df1_dx1 <- 2*x1 + 2
  df1_dx2 <- -2*x2
  
  # f2 = x1² x2 + x2 - 1
  df2_dx1 <- 2*x1*x2
  df2_dx2 <- x1^2 + 1
  
  # Construction de la matrice jacobienne
  J <- matrix(c(df1_dx1, df1_dx2, 
                df2_dx1, df2_dx2), nrow = 2, byrow = TRUE)
  return(J)
}

# Grille pour la visualisation
x1_vals <- seq(-3, 3, length.out = 400)
x2_vals <- seq(-3, 3, length.out = 400)

# Calcul des courbes de niveau (zéro-clines)
z1 <- outer(x1_vals, x2_vals, function(a, b) a^2 - b^2 + 2*a)
z2 <- outer(x1_vals, x2_vals, function(a, b) a^2 * b + b - 1)

# Tracé des courbes de niveau zéro (intersections = solutions)
png(file.path(output_dir, "Ex4_Localisation_solutions.png"), width = 900, height = 700)
contour(x1_vals, x2_vals, z1, levels = 0, col = "blue", lwd = 2.5, 
        main = "Localisation des solutions (Intersections des zero-clines)",
        xlab = "x1", ylab = "x2", cex.main = 1.2)
contour(x1_vals, x2_vals, z2, levels = 0, col = "red", lwd = 2.5, add = TRUE)
grid()
legend("topright", legend = c("f1(x1, x2) = 0", "f2(x1, x2) = 0"), 
       col = c("blue", "red"), lty = 1, lwd = 2.5, cex = 1.1)
dev.off()

cat("Observation graphique :\n")
cat("- Une solution visible proche de (x₁, x₂) ≈ (0.3, 0.9)\n")
cat("- Cette intersection servira de point de départ X₀ pour Newton-Raphson\n\n")

# Point initial choisi d'après le graphique
x0_guess <- c(0.5, 1.0)
cat(sprintf("Point initial choisi : X₀ = (%.1f, %.1f)\n\n", x0_guess[1], x0_guess[2]))

# ============================================================================
# QUESTION ii : Algorithme de Newton-Raphson multi-dimensionnel
# ============================================================================
cat("=== QUESTION ii : Méthode de Newton-Raphson ===\n\n")

# Implémentation de la méthode de Newton pour systèmes
newton_raphson_systeme <- function(x0, tol = 1e-14, max_iter = 20) {
  x <- x0
  errors <- numeric(max_iter)
  iterations <- vector("list", max_iter + 1)
  iterations[[1]] <- x
  
  cat(sprintf("%-6s | %-12s | %-12s | %-14s\n", "Iter", "x1", "x2", "Erreur (norme)"))
  cat(strrep("-", 60), "\n")
  
  for (k in 1:max_iter) {
    # Évaluation de F(x_k)
    F_val <- systeme_equations(x)
    
    # Évaluation de J(x_k)
    J_val <- jacobienne(x)
    
    # Résolution du système linéaire : J(x_k) · Δx = -F(x_k)
    # Plus stable numériquement que d'inverser J explicitement
    delta <- solve(J_val, -F_val)
    
    # Mise à jour : x_{k+1} = x_k + Δx
    x_new <- x + delta
    
    # Calcul de l'erreur (norme euclidienne du pas)
    err <- sqrt(sum(delta^2))
    errors[k] <- err
    iterations[[k + 1]] <- x_new
    
    cat(sprintf("%4d   | %12.8f | %12.8f | %.4e\n", k, x[1], x[2], err))
    
    # Test de convergence
    if (err < tol) {
      errors <- errors[1:k]
      iterations <- iterations[1:(k + 1)]
      cat("\n✓ Convergence atteinte après", k, "itérations\n")
      break
    }
    
    x <- x_new
  }
  
  return(list(solution = x, errors = errors, iterations = iterations))
}

# Résolution du système
resultat <- newton_raphson_systeme(x0_guess)

cat("\n--- SOLUTION FINALE ---\n")
cat(sprintf("x₁ = %.15f\n", resultat$solution[1]))
cat(sprintf("x₂ = %.15f\n", resultat$solution[2]))

# Vérification : évaluation de F à la solution
F_final <- systeme_equations(resultat$solution)
cat(sprintf("\nVérification : ||F(x*)|| = %.4e\n", sqrt(sum(F_final^2))))

# ============================================================================
# Analyse de l'ordre de convergence
# ============================================================================
cat("\n=== Analyse de l'ordre de convergence ===\n\n")

err <- resultat$errors
n <- length(err)

if (n > 2) {
  # Pour vérifier l'ordre quadratique : log(e_{k+1}) vs log(e_k) doit avoir pente ≈ 2
  log_err_curr <- log10(err[1:(n - 1)])
  log_err_next <- log10(err[2:n])
  
  # Estimation de la pente par régression linéaire
  if (length(log_err_curr) > 1) {
    fit <- lm(log_err_next ~ log_err_curr)
    pente <- coef(fit)[2]
    cat(sprintf("Pente estimée (ordre de convergence) : %.3f\n", pente))
    cat("Théorique pour Newton : 2.0 (convergence quadratique)\n\n")
  }
  
  # Graphique log-log
  png(file.path(output_dir, "Ex4_Ordre_convergence.png"), width = 900, height = 700)
  plot(log_err_curr, log_err_next, type = "b", col = "purple", pch = 19, lwd = 2,
       cex = 1.3, cex.lab = 1.2, cex.main = 1.3,
       main = "Ordre de convergence de Newton-Raphson",
       xlab = "log10(Erreur_k)", ylab = "log10(Erreur_k+1)")
  grid()
  
  # Droite théorique de pente 2 (convergence quadratique)
  abline(a = 0, b = 2, col = "darkgreen", lty = 2, lwd = 2.5)
  
  # Droite de régression
  if (exists("fit")) {
    abline(fit, col = "red", lty = 1, lwd = 2)
    legend("topleft", 
           legend = c(sprintf("Données (pente = %.2f)", pente), 
                      "Théorique (pente = 2)"),
           col = c("red", "darkgreen"), 
           lty = c(1, 2), lwd = c(2, 2.5), cex = 1.1)
  } else {
    legend("topleft", legend = "Théorique (pente = 2)",
           col = "darkgreen", lty = 2, lwd = 2.5, cex = 1.1)
  }
  dev.off()
}

# ============================================================================
# Graphique de l'évolution de l'erreur
# ============================================================================
png(file.path(output_dir, "Ex4_Evolution_erreur.png"), width = 900, height = 600)
plot(1:length(err), err, type = "b", col = "steelblue", pch = 19, lwd = 2,
     log = "y", xlab = "Itération k", ylab = "Erreur ||Δx_k|| (échelle log)",
     main = "Convergence de Newton-Raphson", cex.lab = 1.2)
grid()
dev.off()

cat("\n✓ Ex4 terminé : figures enregistrées dans le dossier 'output'.\n")
cat("\nInterprétation :\n")
cat("- La convergence est très rapide (erreur divisée par ~100 à chaque itération)\n")
cat("- Le graphique log-log montre une pente proche de 2 → convergence QUADRATIQUE\n")
cat("- Cela signifie que le nombre de décimales exactes double à chaque itération\n")
