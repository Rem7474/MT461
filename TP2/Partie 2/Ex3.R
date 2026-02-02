# Exercice 3 — Oscillateur harmonique (pendule simple)
# q''(t) = - sin(q(t))
# Réduction en système du 1er ordre: y = (q, p) avec p = q'
# y' = f(y) = (p, -sin(q))

# ============================================================================
# CONFIGURATION: Dossier de sortie des images
# ============================================================================
output_dir <- "TP2/Partie 2/output"
dir.create(output_dir, showWarnings = FALSE)
cat("Les images seront sauvegardées dans:", file.path(getwd(), output_dir), "\n\n")

# Intervalle temporel et pas
h <- 0.001
T0 <- 0
TF <- 4 * pi

# ============================================================================
# OUTILS GÉNÉRIQUES
# ============================================================================
# Intégrateur Euler explicite pour systèmes y' = f(t, y)
euler_system <- function(f, t0, tf, y0, h) {
  n <- floor((tf - t0) / h)
  t <- seq(t0, t0 + n * h, by = h)
  y <- matrix(0, nrow = n + 1, ncol = length(y0))
  y[1, ] <- y0
  for (k in 1:n) {
    y[k + 1, ] <- y[k, ] + h * f(t[k], y[k, ])
  }
  list(t = t, y = y)
}

# Enveloppement d'angle: ramène q dans [-pi, pi]
wrap_angle <- function(theta) {
  ((theta + pi) %% (2 * pi)) - pi
}

# Énergie mécanique du pendule (m = l = g = 1):
# H(q, p) = p^2/2 + (1 - cos(q))
energy <- function(q, p) {
  0.5 * p^2 + (1 - cos(q))
}

# ============================================================================
# DYNAMIQUE DU PENDULE
# ============================================================================
pendulum_f <- function(t, y) {
  q <- y[1]; p <- y[2]
  c(p, -sin(q))
}

# Helper pour tracer une trajectoire en espace des phases (q, p)
plot_phase <- function(sol, main, file, add_separatrix = FALSE) {
  q <- sol$y[, 1]
  p <- sol$y[, 2]
  q_wrapped <- wrap_angle(q)
  
  png(file.path(output_dir, file), width = 900, height = 600)
  plot(q_wrapped, p, type = "l", col = "steelblue", lwd = 2,
       main = main, xlab = "q (rad) [mod 2π]", ylab = "p = dq/dt")
  grid()
  # Option: séparatrice E = 2  ->  p = ± 2 cos(q/2)
  if (add_separatrix) {
    q_grid <- seq(-pi, pi, length.out = 1000)
    lines(q_grid,  2 * cos(q_grid / 2), col = "firebrick", lwd = 2, lty = 2)
    lines(q_grid, -2 * cos(q_grid / 2), col = "firebrick", lwd = 2, lty = 2)
    legend("topright", c("Trajectoire", "Séparatrice E=2"),
           col = c("steelblue", "firebrick"), lwd = c(2, 2), lty = c(1, 2))
  }
  dev.off()
}

# ============================================================================
# QUESTION 3.2 — Deux conditions initiales: y0 = (0,1) et y0 = (0,2)
# ============================================================================
cat("=== Ex3 — Q3.2: intégration par Euler explicite (h = 0.001) ===\n")

# CI 1: y0 = (q0, p0) = (0, 1)
y0_1 <- c(0, 1)
sol1 <- euler_system(pendulum_f, T0, TF, y0_1, h)
plot_phase(sol1,
           main = "Pendule — phase (q,p), Euler explicite, y0=(0,1)",
           file = "Ex3_phase_y0_0_1.png",
           add_separatrix = TRUE)

# CI 2: y0 = (q0, p0) = (0, 2)
y0_2 <- c(0, 2)
sol2 <- euler_system(pendulum_f, T0, TF, y0_2, h)
plot_phase(sol2,
           main = "Pendule — phase (q,p), Euler explicite, y0=(0,2)",
           file = "Ex3_phase_y0_0_2.png",
           add_separatrix = TRUE)

# Interprétation (imprimée dans la console)
E1 <- energy(sol1$y[1,1], sol1$y[1,2])
E2 <- energy(sol2$y[1,1], sol2$y[1,2])
cat(sprintf("E(y0=(0,1)) = %.3f  (< 2  → oscillations bornées)\n", E1))
cat(sprintf("E(y0=(0,2)) = %.3f  (= 2  → sur la séparatrice)\n", E2))
cat("Rappel: E_separatrice = 2 correspond au sommet (q=π, p=0).\n\n")

# ============================================================================
# QUESTION 3.3 — Balayage p0 ∈ [0, 3] avec q0 = 0
# ============================================================================
cat("=== Ex3 — Q3.3: balayage de p0 dans [0,3] ===\n")

p0_values <- seq(0, 3, by = 0.25)
sols <- vector("list", length(p0_values))

png(file.path(output_dir, "Ex3_phase_sweep_p0_0_3.png"), width = 1000, height = 700)
plot(NA, xlim = c(-pi, pi), ylim = c(-3.2, 3.2),
     xlab = "q (rad) [mod 2π]", ylab = "p = dq/dt",
     main = "Pendule — familles de trajectoires pour y0=(0,p0)")
grid()

cols <- colorRampPalette(c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))(length(p0_values))

for (i in seq_along(p0_values)) {
  y0 <- c(0, p0_values[i])
  sol <- euler_system(pendulum_f, T0, TF, y0, h)
  q_wrapped <- wrap_angle(sol$y[, 1])
  lines(q_wrapped, sol$y[, 2], col = cols[i], lwd = 1.8)
}

# Séparatrice (E = 2)
q_grid <- seq(-pi, pi, length.out = 1000)
lines(q_grid,  2 * cos(q_grid / 2), col = "black", lwd = 2, lty = 2)
lines(q_grid, -2 * cos(q_grid / 2), col = "black", lwd = 2, lty = 2)

# Légende avec les valeurs de p0
legend_labels <- c(paste0("p0=", p0_values), "Séparatrice E=2")
legend_cols <- c(cols, "black")
legend_lwd <- c(rep(1.8, length(p0_values)), 2)
legend_lty <- c(rep(1, length(p0_values)), 2)
legend("topright", legend_labels, col = legend_cols, lwd = legend_lwd, 
       lty = legend_lty, bty = "n", cex = 0.8)
dev.off()

# Résumé console
cat("Interprétation physique rapide:\n")
cat("- p0 < 2  → oscillations (courbes fermées autour de (0,0)).\n")
cat("- p0 = 2  → séparatrice (passe par (q=0, p=±2) et (q=π, p=0)).\n")
cat("- p0 > 2  → rotations (le pendule fait des tours complets).\n")

cat("\n✓ Ex3 terminé: figures enregistrées dans le dossier 'output'.\n")
