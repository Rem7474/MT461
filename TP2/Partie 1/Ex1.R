# Équation différentielle non linéaire
# Méthodes de Runge-Kutta - Exercice 1

# Définition du problème de Cauchy:
# y'(x) = 2*x*y(x) - x*y(x)^2
# y(0) = y0
# Solution exacte: y(x) = 2*y0 / (y0 + (2 - y0)*exp(-x^2))

# Fonction pour l'équation différentielle
f <- function(x, y) {
  2 * x * y - x * y^2
}

# Solution exacte
y_exacte <- function(x, y0) {
  2 * y0 / (y0 + (2 - y0) * exp(-x^2))
}

# ============================================================================
# MÉTHODE 1: Euler explicite
# ============================================================================
euler <- function(f, x0, xf, y0, h) {
  n <- floor((xf - x0) / h)
  x <- seq(x0, x0 + n * h, by = h)
  y <- numeric(n + 1)
  y[1] <- y0
  
  for (i in 1:n) {
    y[i+1] <- y[i] + h * f(x[i], y[i])
  }
  
  return(list(x = x, y = y))
}

# ============================================================================
# MÉTHODE 2: Point milieu (Runge-Kutta ordre 2)
# ============================================================================
midpoint <- function(f, x0, xf, y0, h) {
  n <- floor((xf - x0) / h)
  x <- seq(x0, x0 + n * h, by = h)
  y <- numeric(n + 1)
  y[1] <- y0
  
  for (i in 1:n) {
    k1 <- f(x[i], y[i])
    k2 <- f(x[i] + h/2, y[i] + (h/2) * k1)
    y[i+1] <- y[i] + h * k2
  }
  
  return(list(x = x, y = y))
}

# ============================================================================
# MÉTHODE 3: Runge-Kutta ordre 4
# ============================================================================
rk4 <- function(f, x0, xf, y0, h) {
  n <- floor((xf - x0) / h)
  x <- seq(x0, x0 + n * h, by = h)
  y <- numeric(n + 1)
  y[1] <- y0
  
  for (i in 1:n) {
    k1 <- f(x[i], y[i])
    k2 <- f(x[i] + h/2, y[i] + (h/2) * k1)
    k3 <- f(x[i] + h/2, y[i] + (h/2) * k2)
    k4 <- f(x[i] + h, y[i] + h * k3)
    y[i+1] <- y[i] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
  }
  
  return(list(x = x, y = y))
}

# ============================================================================
# CONFIGURATION: Dossier de sortie des images
# ============================================================================
output_dir <- "TP2/Partie 1/output"
dir.create(output_dir, showWarnings = FALSE)
cat("Les images seront sauvegardées dans:", file.path(getwd(), output_dir), "\n\n")

# ============================================================================
# QUESTION 1: Euler avec h=0.1 et h=0.01, y0=1
# ============================================================================
cat("\n=== QUESTION 1: Méthode d'Euler ===\n")
# Objectif: Étudier l'influence du pas h et de la condition initiale y0
y0 <- 1
h1 <- 0.1
h2 <- 0.01
x0 <- 0
xf <- 5

# Résolution avec h = 0.1
sol_euler_h01 <- euler(f, x0, xf, y0, h1)
y_exact_h01 <- y_exacte(sol_euler_h01$x, y0)

# Résolution avec h = 0.01
sol_euler_h001 <- euler(f, x0, xf, y0, h2)
y_exact_h001 <- y_exacte(sol_euler_h001$x, y0)

# Graphique 1: Euler h=0.1 vs solution exacte
png(file.path(output_dir, "Ex1_Euler_h01.png"), width = 800, height = 600)
plot(sol_euler_h01$x, sol_euler_h01$y, type = "b", col = "red", 
     main = "Méthode d'Euler: h = 0.1 vs Solution exacte (y0 = 1)",
     xlab = "x", ylab = "y", lwd = 2, pch = 16)
lines(sol_euler_h01$x, y_exact_h01, col = "blue", lwd = 2, lty = 2)
legend("topright", c("Euler (h=0.1)", "Exacte"), col = c("red", "blue"), 
       lwd = 2, lty = c(1, 2))
grid()
dev.off()

# Graphique 2: Comparaison h=0.1 et h=0.01
png(file.path(output_dir, "Ex1_Euler_h01_vs_h001.png"), width = 800, height = 600)
plot(sol_euler_h01$x, sol_euler_h01$y, type = "b", col = "red", 
     main = "Influence du pas h: Euler avec y0 = 1",
     xlab = "x", ylab = "y", lwd = 2, pch = 16)
lines(sol_euler_h001$x, sol_euler_h001$y, col = "green", lwd = 2, pch = 16, type = "p")
lines(sol_euler_h001$x, y_exact_h001, col = "blue", lwd = 2, lty = 2)
legend("topright", c("Euler (h=0.1)", "Euler (h=0.01)", "Exacte"), 
       col = c("red", "green", "blue"), lwd = 2, lty = c(1, 1, 2))
grid()
dev.off()

# Graphique 3: Variation de y0 de 0 à 2
png(file.path(output_dir, "Ex1_Euler_y0_variation.png"), width = 800, height = 600)
y0_values <- seq(0.1, 2, by = 0.2)
plot(NULL, xlim = c(0, 5), ylim = c(0, 2.2),
     main = "Solutions avec variation de y0 de 0 à 2 (Euler, h=0.1)",
     xlab = "x", ylab = "y", lwd = 2)
colors <- rainbow(length(y0_values))

for (i in seq_along(y0_values)) {
  sol <- euler(f, x0, xf, y0_values[i], h1)
  lines(sol$x, sol$y, col = colors[i], lwd = 2)
}

legend("topright", paste("y0 =", round(y0_values, 2)), 
       col = colors, lwd = 2, cex = 0.8)
grid()
dev.off()

# ============================================================================
# QUESTION 2: Erreur globale en fonction du pas h
# ============================================================================
cat("\n=== QUESTION 2: Analyse de convergence ===\n")
# Objectif: Étudier l'erreur globale en fonction du pas h

h_values <- 10^seq(-6, -1, by = 0.1)  # h de 10^-6 à 10^-1
errors_euler <- numeric(length(h_values))

for (i in seq_along(h_values)) {
  h <- h_values[i]
  sol <- euler(f, x0, xf, y0, h)
  y_exact <- y_exacte(sol$x, y0)
  errors_euler[i] <- max(abs(sol$y - y_exact))
}

# Graphique 4: Erreur en échelle logarithmique
png(file.path(output_dir, "Ex1_Erreur_Euler.png"), width = 800, height = 600)
plot(h_values, errors_euler, type = "b", pch = 1, log = "xy",
     main = "Erreur globale en fonction du pas h (Euler, y0 = 1)",
     xlab = "Pas h", ylab = "Erreur globale", lwd = 2, cex = 1.2)
# Estimation de l'ordre (pente)
fit <- lm(log(errors_euler) ~ log(h_values))
abline(fit, col = "red", lwd = 2, lty = 2)
order_euler <- fit$coefficients[2]
legend("bottomright", paste("Ordre estimé:", round(order_euler, 2)), 
       bty = "o", bg = "white", cex = 1)
grid()
dev.off()

cat(sprintf("Ordre de la méthode d'Euler: %.3f (théorique: 1)\n", order_euler))

# ============================================================================
# QUESTION 3: Comparaison Euler, Point milieu et RK4 (h=0.01, y0=1)
# ============================================================================
cat("\n=== QUESTION 3: Comparaison des trois méthodes ===\n")
# Objectif: Comparer les trois méthodes pour h=0.01 et y0=1
h <- 0.01
sol_euler <- euler(f, x0, xf, y0, h)
sol_midpoint <- midpoint(f, x0, xf, y0, h)
sol_rk4 <- rk4(f, x0, xf, y0, h)

y_exact <- y_exacte(sol_euler$x, y0)

# Graphique 5: Comparaison des trois méthodes
png(file.path(output_dir, "Ex1_Comparaison_methodes.png"), width = 800, height = 600)
plot(sol_euler$x, sol_euler$y, type = "l", col = "red", 
     main = "Comparaison des trois méthodes (h = 0.01, y0 = 1)",
     xlab = "x", ylab = "y", lwd = 2, lty = 1)
lines(sol_midpoint$x, sol_midpoint$y, col = "green", lwd = 2, lty = 2)
lines(sol_rk4$x, sol_rk4$y, col = "purple", lwd = 2, lty = 3)
lines(sol_euler$x, y_exact, col = "blue", lwd = 2.5, lty = 4)
legend("topright", c("Euler", "Point milieu", "RK4", "Exacte"), 
       col = c("red", "green", "purple", "blue"), lwd = 2, 
       lty = c(1, 2, 3, 4))
grid()
dev.off()

# Graphique 6: Erreurs des trois méthodes
png(file.path(output_dir, "Ex1_Erreurs_methodes.png"), width = 800, height = 600)
err_euler <- abs(sol_euler$y - y_exact)
err_midpoint <- abs(sol_midpoint$y - y_exact)
err_rk4 <- abs(sol_rk4$y - y_exact)

plot(sol_euler$x, err_euler, type = "l", col = "red", 
     main = "Erreur locale pour les trois méthodes (h = 0.01, y0 = 1)",
     xlab = "x", ylab = "Erreur |y_approchée - y_exacte|", 
     lwd = 2, lty = 1)
lines(sol_midpoint$x, err_midpoint, col = "green", lwd = 2, lty = 2)
lines(sol_rk4$x, err_rk4, col = "purple", lwd = 2, lty = 3)
legend("topleft", c("Euler", "Point milieu", "RK4"), 
       col = c("red", "green", "purple"), lwd = 2, lty = c(1, 2, 3))
grid()
dev.off()

# ============================================================================
# ORDRE DE CONVERGENCE DES TROIS MÉTHODES
# ============================================================================
cat("\n--- Analyse de l'ordre de convergence ---\n")

h_values <- 10^seq(-4, -2, by = 0.1)
errors_euler <- numeric(length(h_values))
errors_midpoint <- numeric(length(h_values))
errors_rk4 <- numeric(length(h_values))

for (i in seq_along(h_values)) {
  h <- h_values[i]
  
  # Euler
  sol <- euler(f, x0, xf, y0, h)
  y_exact <- y_exacte(sol$x, y0)
  errors_euler[i] <- max(abs(sol$y - y_exact))
  
  # Point milieu
  sol <- midpoint(f, x0, xf, y0, h)
  y_exact <- y_exacte(sol$x, y0)
  errors_midpoint[i] <- max(abs(sol$y - y_exact))
  
  # RK4
  sol <- rk4(f, x0, xf, y0, h)
  y_exact <- y_exacte(sol$x, y0)
  errors_rk4[i] <- max(abs(sol$y - y_exact))
}

# Calcul des ordres
fit_euler <- lm(log(errors_euler) ~ log(h_values))
fit_midpoint <- lm(log(errors_midpoint) ~ log(h_values))
fit_rk4 <- lm(log(errors_rk4) ~ log(h_values))

order_euler <- fit_euler$coefficients[2]
order_midpoint <- fit_midpoint$coefficients[2]
order_rk4 <- fit_rk4$coefficients[2]

cat(sprintf("Ordre de Euler:       %.3f (théorique: 1)\n", order_euler))
cat(sprintf("Ordre du point milieu: %.3f (théorique: 2)\n", order_midpoint))
cat(sprintf("Ordre de RK4:         %.3f (théorique: 4)\n", order_rk4))

# ============================================================================
# CALCUL DE L'ORDRE RK4 SUR UN INTERVALLE RÉDUIT (petits h)
# ============================================================================
cat("\n--- Calcul de l'ordre RK4 sur intervalle réduit (10^-3 à 10^-2) ---\n")

# Utiliser seulement les points avec 10^-3 <= h <= 10^-2
h_values_rk4_reduced <- h_values[h_values >= 10^(-3) & h_values <= 10^(-2)]
errors_rk4_reduced <- errors_rk4[h_values >= 10^(-3) & h_values <= 10^(-2)]

fit_rk4_reduced <- lm(log(errors_rk4_reduced) ~ log(h_values_rk4_reduced))
order_rk4_reduced <- fit_rk4_reduced$coefficients[2]

cat(sprintf("Ordre de RK4 (intervalle 10^-3 à 10^-2): %.3f (théorique: 4)\n", order_rk4_reduced))

# Graphique 8: RK4 seul sur intervalle réduit
png(file.path(output_dir, "Ex1_RK4_ordre_reduit.png"), width = 800, height = 600)
plot(h_values_rk4_reduced, errors_rk4_reduced, type = "b", pch = 3, log = "xy",
     main = "Ordre de convergence RK4 (intervalle réduit: 10^-3 à 10^-2)",
     xlab = "Pas h", ylab = "Erreur globale max", 
     lwd = 2, cex = 1.2, col = "purple")

# Droite de tendance
x_line_rk4 <- range(h_values_rk4_reduced)
y_rk4_reduced <- exp(fit_rk4_reduced$coefficients[1]) * x_line_rk4^order_rk4_reduced
lines(x_line_rk4, y_rk4_reduced, col = "red", lwd = 2.5, lty = 2)

legend("bottomright", 
       c(sprintf("RK4 (ordre %.2f)", order_rk4_reduced),
         "Droite de tendance"),
       col = c("purple", "red"), lwd = 2, lty = c(1, 2), pch = c(3, NA))
grid()
dev.off()

# Graphique 7: Convergence de tous les ordres
png(file.path(output_dir, "Ex1_Ordres_convergence.png"), width = 800, height = 600)
# Calculer les limites pour inclure tous les erreurs
ylim_min <- min(errors_euler, errors_midpoint, errors_rk4)
ylim_max <- max(errors_euler, errors_midpoint, errors_rk4)
plot(h_values, errors_euler, type = "b", pch = 1, log = "xy",
     main = "Ordre de convergence des trois méthodes",
     xlab = "Pas h", ylab = "Erreur globale max", 
     lwd = 2, cex = 1.2, col = "red",
     ylim = c(ylim_min, ylim_max))
points(h_values, errors_midpoint, type = "b", col = "green", lwd = 2, cex = 1.2, pch = 2)
points(h_values, errors_rk4, type = "b", col = "purple", lwd = 2, cex = 1.2, pch = 3)

# Droites de tendance
x_line <- range(h_values)
y_euler <- exp(fit_euler$coefficients[1]) * x_line^(-order_euler)
y_midpoint <- exp(fit_midpoint$coefficients[1]) * x_line^(-order_midpoint)
y_rk4 <- exp(fit_rk4$coefficients[1]) * x_line^(-order_rk4)
lines(x_line, y_euler, col = "red", lwd = 2, lty = 2)
lines(x_line, y_midpoint, col = "green", lwd = 2, lty = 2)
lines(x_line, y_rk4, col = "purple", lwd = 2, lty = 2)

legend("bottomright", 
       c(sprintf("Euler (ordre %.2f)", order_euler),
         sprintf("Point milieu (ordre %.2f)", order_midpoint),
         sprintf("RK4 (ordre %.2f)", order_rk4)),
       col = c("red", "green", "purple"), lwd = 2, pch = c(1, 2, 3))
grid()
dev.off()

cat("\n✓ Tous les graphiques ont été générés.\n")
