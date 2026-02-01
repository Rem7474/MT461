# Exercice 2

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



# Superposition des deux méthodes
plot(nvals, I_f, type="b", col="red", pch=1, ylim=range(c(I_f, I_b)), 
     main="Forward (rouge) vs Backward (bleu)", xlab="n", ylab="I_n")
lines(nvals, I_b, type="b", col="blue", pch=2)
legend("topright", legend=c("Forward", "Backward"), col=c("red", "blue"), pch=c(1,2))

# Tracé de l'erreur absolue
err <- abs(I_f - I_b)
plot(nvals, err, type="b", col="purple", pch=3, 
     main="Erreur absolue |Forward - Backward|", xlab="n", ylab="Erreur absolue")
