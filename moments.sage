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
    return [ (F^r).trace()/(2*g*sqrt(q)^r) for r in range(1,N+1) ]
