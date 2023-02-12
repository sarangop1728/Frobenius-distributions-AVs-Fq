
#############----------
##Create a set of objects corresponding to each simple case
##Both in dimension 1 and dimension 2.
#############----------

#####Class########
class m_case:
    def __init__(self, label, modulus, cong_class):
        self.label=label
        self.modulus = modulus
        self.cong_class=cong_class
    
    def dim(self, dim):
        self.dim=dim #dim = 1 or 2
        
    def m_list(self, m_list):
        self.m_list=m_list;
        
    def q_is_square(self, q_is_square):
        self.q_is_square=q_is_square;
        
    
########Dimension 1 cases for q square   
ec_not_1_mod4_sq = m_case("ec_not_1_mod4_sq", 4, [1]);
ec_not_1_mod4_sq.dim=1
ec_not_1_mod4_sq.m_list=[4]
ec_not_1_mod4_sq.q_is_square=True
##
ec_not_1_mod3_sq = m_case("ec_not_1_mod3_sq", 3, [0, 2])
ec_not_1_mod3_sq.dim=1
ec_not_1_mod3_sq.m_list=[3,6]
ec_not_1_mod3_sq.q_is_square=True;
##
ec_all_sq = m_case("ec_all_sq", 1, [0])
ec_all_sq.dim=1
ec_all_sq.m_list=[1]
ec_all_sq.q_is_square=True

##############
###Dimension 2 cases for q-square

ss_1_mod4_sq = m_case("ss_1_mod4_sq", 4, [1])
ss_1_mod4_sq.dim=2
ss_1_mod4_sq.m_list=[4]
ss_1_mod4_sq.q_is_square=True
##
ss_1_mod3_sq = m_case("ss_1_mod3_sq", 3, [1])
ss_1_mod3_sq.dim=2
ss_1_mod3_sq.m_list=[3, 6]
ss_1_mod3_sq.q_is_square=True
##
ss_not_1_mod5_sq = m_case("ss_not_1_mod5_sq", 5, [0, 2, 3, 4])
ss_not_1_mod5_sq.dim=2
ss_not_1_mod5_sq.m_list=[5, 10]
ss_not_1_mod5_sq.q_is_square=True
##
ss_not_1_mod8_sq = m_case("ss_not_1_mod8_sq", 8, [0, 2, 3, 4, 5, 6, 7])
ss_not_1_mod8_sq.dim=2
ss_not_1_mod8_sq.m_list=[8]
ss_not_1_mod8_sq.q_is_square=True
##
ss_not_1_mod12_sq = m_case("ss_not_1_mod12_sq", 12, [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
ss_not_1_mod12_sq.dim=2
ss_not_1_mod12_sq.m_list=[12]
ss_not_1_mod12_sq.q_is_square=True
################
dim1_cases=[ec_not_1_mod4_sq, ec_not_1_mod3_sq, ec_all_sq]
dim2_cases=[ss_1_mod4_sq, ss_1_mod3_sq, ss_not_1_mod5_sq, ss_not_1_mod8_sq, ss_not_1_mod12_sq];

###########
#FUNCTIONS
###########


def crt_test(A,M): #input:([a b], [m n])#tests if a crt system 1) has a solution; 2) has a prime solution
    if (A[0]-A[1])%gcd(M[0], M[1])==0:
        x_0 = crt(A,M);
        g=gcd(x_0, lcm(M));
        if x_0.is_prime()==True: #x_0 itself is prime
            return True;
        elif gcd(x_0, g)==1: #in this case, there are primes in the AP
            return True;
        elif x_0==0 and lcm(M).is_prime()==True: #in this case, crt returns 0 mod M, but if M is prime, M is also a valid solution.
            return True;
        else: return False;
    else: return False;
        
        
def crt_array_from_class(case_1, case_2): #create input for above function
    A1=case_1.cong_class;
    M1=[case_1.modulus]*len(A1);
    A2=case_2.cong_class;
    M2=[case_2.modulus]*len(A2);
    pairwise_classes=[];
    for i in range(len(A1)):
        for j in range(len(A2)):
            pairwise_classes.append([[A1[i], A2[j]], [M1[i], M2[j]]]);
    return pairwise_classes;

def lcm_builder(case_1, case_2):
    L1=case_1.m_list;
    L2=case_2.m_list;
    L=[];
    for i in range(len(L1)):
        for j in range(len(L2)):
            L.append(lcm(L1[i], L2[j]));
    return L;

#########
##CALCULATING LCMS FOR ExS case
#########

all_lcms=[];
for case1 in dim1_cases:
    for case2 in dim2_cases:
        pair_classes=crt_array_from_class(case1,case2);
        for k in range(len(pair_classes)):
            if crt_test(pair_classes[k][0],pair_classes[k][1])==True:
                #print(case1.label, case2.label);
                #print(pair_classes[k]);
                L_new=lcm_builder(case1, case2);
                all_lcms=all_lcms+L_new;

set_of_lcms=set(all_lcms)

#########
##OUTPUT for q-square
#########
#{3, 4, 5, 6, 8, 10, 12, 15, 20, 24, 30}


