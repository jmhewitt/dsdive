# set coefficients for the polynomial f(x) = coefs[1] * x^n + ... + coefs[n+1]
coefs = c(1,-.1,1,0)

# set plotting range for polynomial
xmin = -1
xmax = 1

# evaluate f(x) via polyval
xseq = seq(from = xmin, to = xmax, length.out = 100)
y = polyval(coef = coefs, x = xseq)

# evaluate and plot f(x) via standard methods
curve(coefs[1] * x^3 + coefs[2] * x^2 + coefs[3] * x + coefs[4], 
      from = xmin, to = xmax, xlab = expression(x), ylab = expression(f(x)))

# overlay polyval results for comparison
lines(xseq, y, col = 2, lty = 2)