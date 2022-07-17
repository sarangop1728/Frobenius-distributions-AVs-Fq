import csv
import os
from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials

# Global objects
prec = 20 # precision
CC = ComplexField(prec)
R.<x> = PolynomialRing(QQ)

'''
Potentially usefull functions
'''

def kill_repeats(mylist):
    """This function takes a list and kills repeats. """
    newlist=[]
    for l in mylist:
        if newlist.count(l)==0:
            newlist.append(l)
    return newlist


def angle_rank(poly):
    """This is the S-unit algorithm for calculating angle rank of a q-weil polynomial 'poly'"""
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


def num_angles(poly, prec=500):
    """Numerical algorithm that calculates the Frobenius angles"""
    myroots = poly.radical().roots(ComplexField(prec))
    angles = [z[0].argument()/RealField(prec)(pi) for z in myroots]
    return [angle for angle in angles if angle>0]


def significant(rel,prec=500):
    """I'm not quite sure what this does"""
    m = min(map(abs,rel))
    if (m+1).exact_log(2)>=sqrt(prec):
        return False
    else:
        if (max(map(abs, rel))+1).exact_log(2) >= sqrt(prec):
            raise RuntimeError("Mixed significance")
        return True


def sage_lindep(angles):
    """Calls lindep function from PARI"""
    rel = gp.lindep(angles)
    return [ Integer(rel[i]) for i in range(1,len(angles)+1)]


def compute_rank(numbers, prec=500):
    """Numerical computation of angle rank"""
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
    """
        We actually have enough functionality at this point to compute the entire group!
        we added 1 to the span of the normalized angles then subtract 1 from the result
    """
    angles = num_angles(mypoly, prec)
    angles = angles + [1]
    #print angles
    return compute_rank(angles,prec)-1


def moments(poly, N):
    """Calculates sequence of the first 'N' normalized traces of Frobenius of the q-weil polynomial 'poly'"""
    g = ZZ(poly.degree()/2)
    q = ZZ((poly.coefficients()[0])**(1/g))
    polyC = poly.base_extend(CC) 
    q_weil_numbers = polyC.roots(ring = CC, multiplicities = False)
    F = diagonal_matrix(q_weil_numbers)
    return [CC((F^r).trace()/(2*g*sqrt(q)^r)) for r in range(1,N+1)]


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
        self.galois_group = poly.splitting_field('z').galois_group()
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

def trace_sequence(poly, N):
    g = ZZ(poly.degree()/2)
    q = ZZ((poly.coefficients()[0])**(1/g))
    polyC = poly.base_extend(CC) 
    q_weil_numbers = polyC.roots(ring = CC, multiplicities = False)
    F = diagonal_matrix(q_weil_numbers)
    return [CC((F^r).trace()/(2*g*sqrt(q)^r)) for r in range(1,N+1)]

 # sort class list by `attribute` 
 #   class_list.sort(key=lambda x: x.attribute, reverse=False)

'''
 a0_to_csv(d,q,n) : creates a .csv file with the sequence a_0 (up to 10^n) for every polynomial in the list of (d,q) Weil polynomials, sorted by angle rank.

'''

def a0_to_csv(d,q,n):
    
    # create list of qPolyData
    poly_list = poly_data_list(d,q)

    # sort by angle rank
    poly_list.sort(key = lambda x : x.angle_rank)

    # trace sequence
    a0 = [trace_sequence(P.poly,10^n) for P in poly_list]

    file_name = 'a0_' + str(d) + '_' + str(q) + '_10^' + str(n) + '.csv'
    new_dir = "./stats/d=" + str(d) + "/q=" + str(q) 
    if not os.getcwd()[-9:] == "q-moments":
        print("Oh you duffer, you should change to the correct directory or I have no idea where you are!!")
    else:
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)

    path = new_dir + '/' + file_name

    # write csv file
    with open(path, 'w') as F:
        writer = csv.writer(F)
        for row in a0:
            writer.writerow(row)

'''
polys_to_txt(d,q,n) : creates a .txt file with the list of (d,q) Weil polynomials, sorted by angle rank.

'''

def polys_to_txt(d,q):
    
    # create list of qPolyData
    poly_list = poly_data_list(d,q)

    # sort by angle rank
    poly_list.sort(key = lambda x : x.angle_rank)

    # trace sequence
    string_poly = [str(P.poly) for P in poly_list]

    # write csv file
    file_name = 'polys_' + str(d) + '_' + str(q) + '.txt'

    new_dir = "./stats/d=" + str(d) + "/q=" + str(q) 
    if not os.getcwd()[-9:] == "q-moments":
        print("Oh you duffer, you should change to the correct directory or I have no idea where you are!!")
    else:
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)

    path = new_dir + '/' + file_name
    
    with open(path, 'w') as F:
        for row in string_poly:
            F.write(row)
            F.write('\n')
