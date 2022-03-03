import csv
from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials

# Global objects
prec = 200 # precision
CC = ComplexField(prec)
R.<x> = PolynomialRing(QQ)

'''
Potentially usefull functions
'''

# This function takes a list and kills repeats. 
def kill_repeats(mylist):
    newlist=[]
    for l in mylist:
        if newlist.count(l)==0:
            newlist.append(l)
    return newlist

# This is the S-unit algorithm for calculating angle rank of a q-weil polynomial 'poly'
def angle_rank(poly):
    #p is the prime dividing q
    p = radical(poly.coefficients()[0])
    K.<a> = poly.splitting_field()
    l = [p] + [i[0] for i in poly.roots(K)]
    S = K.primes_above(p)
    UGS = UnitGroup(K, S = tuple(S), proof=False)
    ## Even with proof=False, it is guaranteed to obtain independent S-units; just maybe not the fully saturated group.
    d = K.number_of_roots_of_unity()
    gs = [K(i) for i in UGS.gens()]
    l2 = [UGS(i^d).exponents() for i in l] #For x = a^1b^2c^3 exponents are (1,2,3)
    for i in range(len(l)):
        assert(l[i]^d == prod(gs[j]^l2[i][j] for j in range(len(l2[i]))))
    M = Matrix(l2)
    return M.rank()-1

# Numerical algorithm that calculates the Frobenius angles
def num_angles(poly, prec=500):
    myroots = poly.roots(ComplexField(prec))
    angles = [z[0].argument() for z in myroots]
    return [angle for angle in angles if angle>0]

# I'm not quite sure what this does
def significant(rel,prec=500):
    m = min(map(abs,rel))
    if (m+1).exact_log(2)>=sqrt(prec):
        return False
    else:
        if (max(map(abs, rel))+1).exact_log(2) >= sqrt(prec):
            raise RuntimeError("Mixed significance")
        return True

# Calls lindep function from PARI
def sage_lindep(angles):
    rel = gp.lindep(angles)
    return [ Integer(rel[i]) for i in range(1,len(angles)+1)]

# Numerical computation of angle rank
def compute_rank(numbers, prec=500):
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

# Numerical calculation of the angle rank of a q-weil polynomial
def num_angle_rank(mypoly, prec=500):
    angles = num_angles(mypoly, prec)
    angles = angles + [1]
    return compute_rank(angles,prec)-1

def is_easy_sn_an(f, num_trials = 50, assume_irreducible = False):
   """
   Use the Davenport-Smith test to attempt to certify that `f` has Galois group A_n or S_n.
   
   Return 1 if the Galois group is certified as S_n, 2 if A_n, or 0 if no conclusion is reached.
   """
   if not assume_irreducible and not f.is_irreducible():
      return 0
   d = f.degree()
   confirm_sn = False
   for p in primes_first_n(num_trials):
       fp = f.change_ring(GF(p))
       l = fp.factor()
       if (len(l)-d)%2 == 1:
           confirm_sn = True
       g = l[-1][0]
       d1 = g.degree()
       if (d1 <= 7 and (d,d1) in ((1,1),(2,1),(3,2),(3,3),(4,3),(5,3),(5,4),(6,5),(7,5))) or\
           (d1 > d/2 and d1 < d-2 and d1.is_prime()):
           return (2 if (not confirm_sn and f.disc().is_square()) else 1)
   return 0

'''
This might be useful later:
'''
#for g in [2]:
#  for q in [2,3]:
#    R.<T> = ZZ[]
#    L = R.weil_polynomials(2*g, q) # iterates over q-weil polynomials of degree 2*g
#    it = iter(L)
#    for f in it:
#      bool = f.is_irreducible() 
#      if bool:
#       h = f.trace_polynomial()[0]
#       G = magma.GaloisGroup(f)
#       H = magma.GaloisGroup(h)
#       sl=f.newton_slopes(radical(q))
#       codeSize = magma.Order(G)/magma.Order(H)       
#       bool2 = codeSize<2^g       
#       labelG = magma.TransitiveGroupIdentification(G)
#       labelH = magma.TransitiveGroupIdentification(H)        
#       if True:
#        r = num_angle_rank(f)
#        print([g,q, magma.GroupName(G), r])


'''
The poly_data class
'''

class qPoly_Data:
    """ A class which attaches data to a polynomial"""

    def __init__(self, poly):
        """Initialise"""
        self.poly = poly
        # self.is_irreducible = poly.is_irreducible()
        # self.galois_group = poly.galois_group()
        self.dimension = ZZ(poly.degree()/2)
        # self.code_size = []
        # self.q =ZZ(CC((poly.coefficients()[0])^{1/self.dimension}))
        self.angle_rank = num_angle_rank(poly)
        # self.trace_poly = poly.trace_polynomial()[0]
        
    def __repr__(self):
        """Repr"""
        return "q-weil data for the polynomial {p}".format(p = self.poly)


# List with entries given by Poly_Data of q-Weil polynomaisl of degree d.
def poly_data_list(d,q):
    l = R.weil_polynomials(d,q)
    return [qPoly_Data(p) for p in l]


"""
Given a q-Weil polynomial (over QQ) 'poly' and a length 'N', moments(poly, n) calculates the first N normalized traces of Frobenius sequence.

"""

def moments(poly, N):
    g = ZZ(poly.degree()/2)
    q = ZZ((poly.coefficients()[0])**(1/g))
    polyC = poly.base_extend(CC) 
    q_weil_numbers = polyC.roots(ring = CC, multiplicities = False)
    F = diagonal_matrix(q_weil_numbers)
    return [CC((F^r).trace()/(2*g*sqrt(q)^r)) for r in range(1,N+1)]

 # sort class list by `attribute` 
 #   class_list.sort(key=lambda x: x.attribute, reverse=False)

# Example:
# Lets calculate the moment distributions for all degree 4 3-weil polynomials.
# 1. Define the list
weil_4_3 = poly_data_list(4,3)
# 2. Sort the list by increasing order of angle rank
weil_4_3.sort(key = lambda x : x.angle_rank)
[x.angle_rank for x in weil_4_3]
# 3. Define list with all moment sequences (up to some bound N)
N = 1000
moments_4_3 = [moments(x.poly,N) for x in weil_4_3]
# 4. write csv file
with open ('weil_4_3.csv', 'w') as F:
    writer = csv.writer(F)
    for m in moments_4_3:
        writer.writerow(m)
# 5. run python3 histograms.py

