# Exercice 3

output_dir <- "TP1/Partie 2/output"
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# Point fixe et méthode de Newton pour résoudre cos(x) = x
# ============================================================================
f <- function(x) cos(x)
fprime <- function(x) -sin(x)

# Point fixe approximatif (solution de cos(x) = x)
x_star <- 0.7390851332151607  # valeur connue

# Vérification que |f'(x*)| < 1
cat("f'(x*) =", fprime(x_star), "\n")
cat("|f'(x*)| =", abs(fprime(x_star)), "< 1, donc contraction locale\n")

# Test de convergence pour différents x0
x0_values <- c(-5, 0, 2, 10)
for(x0 in x0_values) {
  x <- x0
  for(i in 1:10) {
    x <- cos(x)
  }
  cat("x0 =", x0, "=> x_10 =", x, "\n")
}

# ============================================================================
# Question 2
# ============================================================================

# Analyse de convergence - Point fixe
point_fixe_convergence <- function(x0, max_iter = 30) {
  x <- numeric(max_iter + 1)
  x[1] <- x0
  x_star <- 0.7390851332151607
  
  for(i in 2:(max_iter + 1)) {
    x[i] <- cos(x[i-1])
  }
  
  # Calcul des erreurs
  errors <- abs(x - x_star)
  
  # Analyse de l'ordre de convergence
  ratios <- numeric(max_iter - 5)
  for(i in 6:max_iter) {
    if(errors[i] > 0 && errors[i+1] > 0) {
      ratios[i-5] <- errors[i+1] / errors[i]
    }
  }
  
  return(list(x = x, errors = errors, ratios = ratios))
}

# Test avec x0 = 1
result_pf <- point_fixe_convergence(1)

# Graphiques
png(file.path(output_dir, "Q2_PointFixe.png"))
par(mfrow=c(2,1))
plot(1:length(result_pf$errors), log10(result_pf$errors), type="b",
     main="Erreur (échelle log) - Point fixe", xlab="Itération", ylab="log10(erreur)")

plot(6:30, result_pf$ratios, type="b", 
     main="Ratio d'erreurs consécutives", xlab="Itération", ylab="e_(n+1)/e_n")
dev.off()

cat("Convergence linéaire avec facteur ≈", mean(result_pf$ratios, na.rm=TRUE), "\n")



# ============================================================================
# Question 3
# ============================================================================

# Méthode de Newton
newton_convergence <- function(x0, max_iter = 10) {
  x <- numeric(max_iter + 1)
  x[1] <- x0
  x_star <- 0.7390851332151607
  
  for(i in 2:(max_iter + 1)) {
    fx <- cos(x[i-1]) - x[i-1]
    fpx <- -sin(x[i-1]) - 1
    x[i] <- x[i-1] - fx/fpx
  }
  
  # Calcul des erreurs
  errors <- abs(x - x_star)
  
  # Test de convergence quadratique
  ratios_quad <- numeric(max_iter - 2)
  for(i in 3:max_iter) {
    if(errors[i-1]^2 > 0 && errors[i] > 0) {
      ratios_quad[i-2] <- errors[i] / errors[i-1]^2
    }
  }
  
  return(list(x = x, errors = errors, ratios_quad = ratios_quad))
}

# Test Newton
result_newton <- newton_convergence(1)
png(file.path(output_dir, "Q3_Newton.png"))
par(mfrow=c(2,1))
plot(1:length(result_newton$errors), log10(result_newton$errors), type="b",
  main="Erreur (échelle log) - Newton", xlab="Itération", ylab="log10(erreur)")

plot(3:10, result_newton$ratios_quad, type="b",
  main="Test convergence quadratique: e_(n+1)/e_n^2", 
  xlab="Itération", ylab="e_(n+1)/e_n^2")
dev.off()

cat("Convergence quadratique, constante ≈", mean(result_newton$ratios_quad, na.rm=TRUE), "\n")


# ============================================================================
# Question 4
# ============================================================================

# Point fixe pour x = 10*cos(x)
point_fixe_10cos <- function(x0, max_iter = 20) {
  x <- numeric(max_iter + 1)
  x[1] <- x0
  
  for(i in 2:(max_iter + 1)) {
    x[i] <- 10 * cos(x[i-1])
  }
  
  return(x)
}
png(file.path(output_dir, "Q4_10cos.png"))
x0_values <- c(0, 1, 2, 5)
par(mfrow=c(2,2))

for(x0 in x0_values) {
  result <- point_fixe_10cos(x0)
  plot(0:20, result, type="b", 
       main=paste("x0 =", x0), xlab="Itération", ylab="x_n")
}
dev.off()

cat("Comportement observé: divergence ou oscillation selon x0\n")
cat("Cause: |f'(x)| = 10|sin(x)| peut être > 1, pas de contraction\n")

# ============================================================================
# Question 5 : Accélération d'Aitken
# ============================================================================
aitken_acceleration <- function(x_seq) {
  n <- length(x_seq)
  x_acc <- numeric(n-2)
  
  for(i in 1:(n-2)) {
    delta1 <- x_seq[i+1] - x_seq[i]
    delta2 <- x_seq[i+2] - x_seq[i+1]
    
    if(abs(delta2 - delta1) > 1e-15) {
      x_acc[i] <- x_seq[i] - delta1^2 / (delta2 - delta1)
    } else {
      x_acc[i] <- x_seq[i]
    }
  }
  return(x_acc)
}

# Application à cos(x) = x
x0 <- 1
x_pf <- numeric(20)
x_pf[1] <- x0
for(i in 2:20) {
  x_pf[i] <- cos(x_pf[i-1])
}

png(file.path(output_dir, "Q5_Aitken.png"))
par(mfrow=c(1,2))
plot(1:20, log10(abs(x_pf - x_star)), type="b", col="red",
  main="Point fixe vs Aitken", xlab="Itération", ylab="log10(erreur)")
lines(1:18, log10(abs(x_aitken - x_star)), type="b", col="blue")
legend("topright", c("Point fixe", "Aitken"), col=c("red", "blue"), lty=1)
dev.off()

# Pour x = 10*cos(x) - généralement diverge, Aitken ne peut pas aider
cat("Aitken accélère la convergence pour cos(x)=x mais ne peut pas\n")
cat("stabiliser une suite divergente comme x=10*cos(x)\n")

# ============================================================================
# Question 6 : Méthode d'Aitken-Steffensen
# ============================================================================
steffensen <- function(f, x0, max_iter = 10, tol = 1e-12) {
  x <- numeric(max_iter + 1)
  x[1] <- x0
  
  for(i in 2:(max_iter + 1)) {
    y0 <- x[i-1]
    y1 <- f(y0)
    y2 <- f(y1)
    
    # Formule de Steffensen
    denom <- y2 - 2*y1 + y0
    if(abs(denom) > tol) {
      x[i] <- y0 - (y1 - y0)^2 / denom
    } else {
      x[i] <- y1  # retombe sur point fixe classique
    }
    
    if(abs(x[i] - x[i-1]) < tol) break
  }
  
  return(x[1:i])
}

# Application à cos(x) = x
f_cos <- function(x) cos(x)
result_steff <- steffensen(f_cos, 1)

png(file.path(output_dir, "Q6_Methods_Comparison.png"))
par(mfrow=c(1,1))
plot(1:length(result_newton$errors), log10(result_newton$errors), 
     type="b", col="green", main="Comparaison des méthodes",
     xlab="Itération", ylab="log10(erreur)")
lines(1:20, log10(abs(x_pf - x_star)), type="b", col="red")
lines(1:length(errors_steff), log10(errors_steff), type="b", col="blue")
legend("topright", c("Newton", "Point fixe", "Steffensen"), 
  col=c("green", "red", "blue"), lty=1)
dev.off()

cat("Steffensen: convergence quadratique sans calcul de dérivée\n")
cat("Plus rapide que point fixe, comparable à Newton\n")
