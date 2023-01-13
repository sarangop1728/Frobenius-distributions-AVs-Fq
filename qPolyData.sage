import csv
import os
import time


# Global objects.

prec = 500
CC = ComplexField(prec)
RR = RealField(prec)
R.<x> = PolynomialRing(QQ)

import csv
import os

# Load angle ranks functions.
load('angle_ranks.sage')


# LMFDB labels.
# ______________________________________________________________________________

# Writting numbers in base 26
a_lowercase = ord('a') # Gives the unicode number of the character.
alphabet_size = 26

def _decompose(number):
    '''Returns a generator object with the digits from `number` in base alphabet'''

    while number:
        number, remainder = divmod(number, alphabet_size)
        yield remainder

def base_10_to_alphabet(number):
    '''Converts base 10 integer to its base alphabet representation. If <= 0; the word starts with an 'a'.'''

    if number <= 0:
        word = 'a'
        number = -number
    else:
        word = ''

    return word + ''.join(
            chr(a_lowercase + part)
            for part in _decompose(number)
    )[::-1] # string[::-1] reverses the string.

def base_alphabet_to_10(letters):
    '''Convert an alphabet number (string) to its decimal representation'''
    sgn = 1
    if letters[0] == 'a' and len(letters)>1:
        letters = letters[1:]
        sgn = -1

    number = sgn*sum((ord(letter) - a_lowercase)*alphabet_size**i for i, letter in enumerate(reversed(letters)))

    return number

# q-Weil polynomial to LMFDB label.

''''
    Recall that the label format is g.q.iso, where
    g : the dimension of the AV contained in the isogeny class
    q : the cardinality of the base field
    iso : the integer coefficients a_1,...,a_g of the L-polynomial polynomial in base 26 in the symbols a,b,c,...,z with a=0 and separated by underscores
'''

def lmfdb_label(poly, q):

    g = poly.degree()//2 # Dimension.
    coefficients = poly.coefficients(sparse=False)[2*g-1:g-1:-1] # [a_1,...,a_g].
    iso = '_'.join(base_10_to_alphabet(Integer(coeff)) for coeff in coefficients)

    return str(g) + '.' + str(q) + '.' + iso

# Returns Frobenius polynomial from label.
def poly_from_label(label):
    g = ZZ(label.split(".")[0])
    q = ZZ(label.split(".")[1])
    iso = label[4:]
    coeffs = [base_alphabet_to_10(letters) for letters in list(iso.split('_'))]
    all_coeffs = [1] + coeffs
    reverse = all_coeffs[-2::-1]
    all_coeffs += [q*q^i*reverse[i] for i in range(g)]
    return R(all_coeffs[::-1])



# qPolyClass
# ______________________________________________________________________________

'''
qPolyClass: Python class with all the data associated to a Weil polynomial we can get.

INPUT:  a q-Weil polynomial qpoly.

OUTPUT: a Python class with the following list of attributes:

- poly: the q-Weil polynomial in \mathbf{Q}[x].
- dimension: integer g, the dimension of an abelian variety in the corresponding isogeny class.
- q: corresponding prime power.
- name: string conversion of qpoly.
- is_irreducible: boolean, true ifqpolyis irreducible.
- galois_group: Galois group of the Galois closure of qpoly.
- roots: complex roots of the Weil polynomial.
- angle_rank: integer, angle rank of qpoly.
- label: LMFDB label.
'''

class qPolyClass:
    def __init__(self, qpoly):
        '''Initialise'''

        self.poly = qpoly
        self.p, self.a = ZZ((self.poly).coefficients()[0]).perfect_power()
        # q^g = p^a -> a = rg.
        self.g = self.poly.degree() // 2
        self.r = self.a // self.g
        self.q = (self.p)^(self.r)
        self.name = str(self.poly)
        self.is_irreducible = self.poly.is_irreducible()
        self.galois_group = self.poly.splitting_field('z').galois_group() # 'z' the generator of the Galois closure.
        self.roots = self.poly.base_extend(CC).roots(ring = CC, multiplicities = False)
        self.angle_rank = num_angle_rank(self.poly)
        self.trace_poly = self.poly.trace_polynomial()[0]
        self.label = lmfdb_label(self.poly, self.q)

    def __repr__(self):
        '''Repr'''
        return 'q-Weil data for {p}'.format(p = self.label)


def poly_class_list(d,q, with_extremal_angle_ranks = True):
    '''list of all qPolyClass data of 
    ension d and prime power q.'''

    poly_list =  [qPolyClass(p) for p in R.weil_polynomials(d,q)]

    # Filter out the extremal angle ranks.
    if not with_extremal_angle_ranks:
        poly_list = list(filter(lambda x : x.angle_rank != 0 and x.angle_rank != x.g, poly_list))

    poly_list.sort(key = lambda x : x.angle_rank) # Sort by angle rank.

    return poly_list


# a1 sequence & moments
# ______________________________________________________________________________

def a1_sequence(qpolyClass, N = 10^6):
    '''Generator sequence of normalized traces of frobenius. Default length = 10^6.'''

    g = qpolyClass.g
    q = qpolyClass.q
    label = qpolyClass.label
    q_weil_numbers = qpolyClass.roots
    F = diagonal_matrix(q_weil_numbers)

    return (RR((F^r).trace()/(g*sqrt(q)^r)) for r in range(1,N+1))

def moments(sequence, N=10):
    '''Returns the first N k-moments of the sequence, N=10 by default.'''

    seq_list = list(sequence)
    moments = []
    if len(seq_list) > 0:
        for k in range(N):
            k_pwr_sequence = [s^k for s in seq_list]
            approx = round((sum(k_pwr_sequence)/len(seq_list)),3)
            if abs(approx) < 0.001:
                approx = round(0.0,1)
            moments.append(approx)
    return moments


# Exporting data.
# ______________________________________________________________________________

def a1_to_csv(d,q,n,extremal=True):

    # Create list of qPolyData.
    poly_list = poly_class_list(d,q,extremal)

    # Dictionary of a_1 sequences with lmfdb label as keys.
    a1_dict = {P.label : list(a1_sequence(P,10^n)) for P in poly_list}

    # Header.
    fieldnames = ['Label', 'a_1']

    # Filename
    file_name = 'a1_' + str(d) + '_' + str(q) + '_10^' + str(n) + '.csv'

    # Directory
    new_dir = './stats/d=' + str(d) + '/q=' + str(q) 
    if not os.path.exists(new_dir):
            os.makedirs(new_dir)
    path = new_dir + '/' + file_name

    # Write csv file.
    with open(path, 'w') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()

        for label in a1_dict.keys():
            writer.writerow({'Label' : label, 'a_1': a1_dict[label]})
            timestamp = time.time()
            print('label ' + label + ' written at ' + str(timestamp))

# single label sequence
def a1_label_csv(label, n):

    # Create qPolyData.
    poly = qPolyClass(poly_from_label(label))
    d = 2*ZZ(label.split(".")[0])
    q = ZZ(label.split(".")[1])

    # Dictionary of a_1 sequences with lmfdb label as keys.
    a1_dict = {label : list(a1_sequence(poly,10^n))}

    # Header.
    fieldnames = ['Label', 'a_1']

    # Filename
    file_name = 'a1_' + label + '_10^' + str(n) + '.csv'

    # Directory
    new_dir = './stats/d=' + str(d) + '/q=' + str(q) 
    if not os.path.exists(new_dir):
            os.makedirs(new_dir)
    path = new_dir + '/' + file_name

    # Write csv file.
    with open(path, 'w') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()

        for lab in a1_dict.keys():
            writer.writerow({'Label' : lab, 'a_1': a1_dict[lab]})



def moments_to_csv(d,q,extremal=True):

    # Create list of qPolyData.
    poly_list = poly_class_list(d,q,extremal)

    # Dictionary of trace sequences.
    moments_dict = {P.label : moments(a1_sequence(P,10^5)) for P in poly_list}

    # Header.
    fieldnames = ['Label', 'a_1 moments']

    # Filename.
    file_name = str(d) + '_' + str(q) + '_moments.csv'

    # Directory.
    new_dir = './stats/d=' + str(d) + '/q=' + str(q) 
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    path = new_dir + '/' + file_name

    # Write csv file.
    with open(path, 'w') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()

        for label in moments_dict.keys():
            writer.writerow({'Label' : label, 'a_1 moments': moments_dict[label]})

def labels_to_txt(d,q,extremal=True):

    # Create list of qPolyData.
    poly_list = poly_class_list(d,q,extremal)

    # List of trace sequences.
    labels = [P.label for P in poly_list]

    # Filename.
    file_name = str(d) + '_' + str(q) + '_labels.txt'

    # Directory.
    new_dir = './stats/d=' + str(d) + '/q=' + str(q) 
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    path = new_dir + '/' + file_name

    # Writting the file.
    with open(path, 'w') as F:
        for row in labels:
            F.write(row)
            F.write('\n')


# Generate files.
# ______________________________________________________________________________


labels = ['1.2.ab']

for label in labels:
    for n in [2,3,4,5,6]:
        a1_label_csv(label,n)

