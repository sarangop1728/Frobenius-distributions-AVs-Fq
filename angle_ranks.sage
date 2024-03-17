# The original version of this code was shared to us by David Zureick-Brown
from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials

def num_angles(poly, prec=500):
    '''Numerical algorithm that calculates the Frobenius angles'''
    myroots = poly.radical().roots(ComplexField(prec))
    angles = [z[0].argument()/RealField(prec)(pi) for z in myroots]
    return [angle for angle in angles if angle>0]

def significant(rel,prec=500):
    '''Checks significance'''
    m = min(map(abs,rel))
    if (m+1).exact_log(2)>=sqrt(prec):
        return False
    else:
        if (max(map(abs, rel))+1).exact_log(2) >= sqrt(prec):
            raise RuntimeError('Mixed significance')
        return True

def sage_lindep(angles):
    '''Calls lindep function from PARI'''
    rel = gp.lindep(angles)
    return [Integer(rel[i]) for i in range(1,len(angles)+1)]

def compute_rank(numbers, prec=500):
    '''Numerical computation of angle rank'''
    r = len(numbers)
    if r == 1:
        return 1
    else:
        rels = sage_lindep(numbers)
        if significant(rels, prec):
            i=0
            while i<len(rels):
                if rels[i] != 0:
                    numbers.pop(i)
                    return compute_rank(numbers, prec)
                else:
                    i+=1
        else:
            return len(numbers)

def num_angle_rank(mypoly,prec=500):
    '''We added 1 to the span of the normalized angles then subtract 1 from the result.'''
    angles = num_angles(mypoly, prec)
    angles = angles + [1]
    return compute_rank(angles,prec)-1
