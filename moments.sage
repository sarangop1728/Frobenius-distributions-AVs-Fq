prec = 20 # precision
CC = ComplexField(prec)
R.<x> = PolynomialRing(QQ)

"""
Given a q-Weil polynomial (over QQ) 'poly' and a length 'N', moments(poly, n) calculates the first N normalized traces of Frobenius sequence.

"""

def moments(poly, N):
    g = ZZ(poly.degree()/2)
    q = ZZ((poly.coefficients()[0])**(1/g))
    polyC = poly.base_extend(CC) 
    q_weil_numbers = polyC.roots(ring = CC, multiplicities = False)
    F = diagonal_matrix(q_weil_numbers)
    return [CC((F^r).trace()/(2*g*sqrt(q)^r)) for r in range(1,N+1) ]

# EXAMPLE:
# 1. Take the first 5-weil poly of degree 4:
weil(4,5) = R.weil_polynomaials(4,5)
f = weil(4,5)[0]
moments_4_5_0 = moments(f, 1000)

# 2. copy the output
# 3. paste into the 'data' variable in the histograms.py file
# 4. run: python3 histograms.py in the terminal 
