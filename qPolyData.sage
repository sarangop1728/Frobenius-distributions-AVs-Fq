import csv # for exporting .csv files
import os # file management


# global objects
prec = 100 
CC = ComplexField(prec)
RR = RealField(prec)
R.<x> = PolynomialRing(QQ)

# load angle ranks functions
load("angle_ranks.sage")

#________________________________________________________________________

'''
qPolyClass: Python class with all the data associated to a Weil polynomial we can get.

INPUT:  a q-Weil polynomial qpoly.

OUTPUT: a Python class with the following list of attributes:

    - poly:           the q-Weil polynomial in R.
    - dimension:      integer g, the dimension of a corresponding abelian variety.
    - q:              corresponding prime power.
    - name:           string conversion of qpoly.
    - is_irreducible: boolean, true if qpoly is irreducible.
    - galois_group:   Galois group of the Galois closure of qpoly.
    - angle_rank:     integer, angle rank of qpoly.

'''

class qPolyClass:
    def __init__(self, qpoly):
        """Initialise"""
        self.poly = qpoly
        self.dim = self.poly.degree() // 2
        self.q = ZZ(RR((self.poly).coefficients()[0])**(1/self.dim))
        self.name = str(self.poly)
        self.is_irreducible = self.poly.is_irreducible()
        # we call 'z' the generator of the Galois closure
        self.galois_group = self.poly.splitting_field('z').galois_group()
        self.angle_rank = num_angle_rank(self.poly)
        self.trace_poly = self.poly.trace_polynomial()[0]
        # self.code_size = []
        
        
    def __repr__(self):
        """Repr"""
        return "q-Weil data for the polynomial {p}".format(p = self.poly)

#________________________________________________________________________

'''
qPolyClass_list: list of all qPolyClass data of dimension d

INPUT: d = dimension, q = prime power.

OUTPUT: list of all dimension d qPolyClass, sorted by angle rank.

'''
    
def qPolyClass_list(d,q):
    l = R.weil_polynomials(d,q)
    poly_list =  [qPolyClass(p) for p in l]
    # sort by angle rank
    poly_list.sort(key = lambda x : x.angle_rank)
    return poly_list

#________________________________________________________________________

'''
a0_sequence: sequence of geometric normalized traces of Frobenius

INPUT: qpolyClass, q-Weil poly class. N, the length of the sequence.

OUTPUT: list (x_1, x_2, x_3, ..., x_N) of normalized traces of Frobenius. N=10^4 by default.

'''

def a0_sequence(qpolyClass, N = 10^4):
    g = qpolyClass.dim
    q = qpolyClass.q
    polyC = qpolyClass.poly.base_extend(CC) 
    q_weil_numbers = polyC.roots(ring = CC, multiplicities = False)
    F = diagonal_matrix(q_weil_numbers)
    return [RR((F^r).trace()/(2*g*sqrt(q)^r)) for r in range(1,N+1)]

#________________________________________________________________________

'''
k-moments: given a sequence (x_1, x_2, ..., x_N) calculates the average of
           the sequence (x_1^k, x_2^k, ..., x_N^k).

INPUT: sequence, a list with numerical entries. k, an integer.

OUTPUT: the average of the k-power sequence, in  RR.
'''

def k_moments(sequence, k):
    k_pwr_sequence = [RR(s)^k for s in sequence]
    return sum(k_pwr_sequence)/len(sequence)

#________________________________________________________________________

'''
moments: calculates the first N moments of a sequence. 

INPUT: sequence, a sequence with numerical values. N = number of moments to display (default is 10)

OUTPUT: list with RR entries of length N. We display only 3 digits after the decimal point.

'''

def moments(sequence, N = 10):
    return ["{:.3f}".format(k_moments(sequence,k)) for k in range(N)]

#________________________________________________________________________

'''
all_moments: list of moments for (d,q)-Weil polynomials.

INPUT: d = dimension, q = prime power.

OUTPUT: list of (the first 10) moments, for every (d,q)-Weil polynomial.

'''

def all_moments(d,q):
    poly_list = qPolyClass_list(d,q)
    return [moments(a0_sequence(P)) for P in poly_list]
    
#________________________________________________________________________

'''
 a0_to_csv(d,q,n) : creates a .csv file with the sequence a_0 (up to 10^n) for every polynomial in the list of (d,q) Weil polynomials, sorted by angle rank.

'''

def a0_to_csv(d,q,n):
    
    # create list of qPolyData
    poly_list = qPolyClass_list(d,q)

    # trace sequence
    a0 = [a0_sequence(P,10^n) for P in poly_list]

    # save file in correct directory
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

#________________________________________________________________________

'''
polys_to_txt(d,q,n) : creates a .txt file with the list of (d,q)-Weil polynomials, sorted by angle rank.

'''

def polys_to_txt(d,q):
    
    # create list of qPolyData
    poly_list = qPolyClass_list(d,q)

    # trace sequence
    string_poly = [P.name for P in poly_list]

    # write txt file
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

#________________________________________________________________________

'''
moments_to_txt(d,q,n) : creates a .txt file with the list of (d,q)-Weil moments, sorted by angle rank.

'''

def moments_to_txt(d,q):

    # trace sequence
    string_moments = [str(moments) for moments in all_moments(d,q)]

    # write txt file
    file_name = 'moments_' + str(d) + '_' + str(q) + '.txt'

    new_dir = "./stats/d=" + str(d) + "/q=" + str(q) 
    if not os.getcwd()[-9:] == "q-moments":
        print("Oh you duffer, you should change to the correct directory or I have no idea where you are!!")
    else:
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
    path = new_dir + '/' + file_name
    
    with open(path, 'w') as F:
        for row in string_moments:
            F.write(row)
            F.write('\n')
