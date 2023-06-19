###########################################################
##### Compute AR reciprocal roots given the AR coefficients
###########################################################
# Assume the folloing AR coefficients for an AR(8)
phi=c(0.27, 0.07, -0.13, -0.15, -0.11, -0.15, -0.23, -0.14)
roots=1/polyroot(c(1, -phi)) # compute reciprocal characteristic roots
r=Mod(roots) # compute moduli of reciprocal roots
lambda=2*pi/Arg(roots) # compute periods of reciprocal roots
# print results modulus and frequency by decreasing order
print(cbind(r, abs(lambda))[order(r, decreasing=TRUE), ][c(2,4,6,8),]) 
