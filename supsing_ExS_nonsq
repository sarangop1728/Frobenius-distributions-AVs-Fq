
#############----------
##Create a set of objects corresponding to each simple case
##Both in dimension 1 and dimension 2.
#############----------

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
        
    
########Dimension 1 cases for q not square   
ec_2 = m_case("ec_2", 2, [0]); #the only prime that's 0 mod 2 is 2
ec_2.dim=1
ec_2.m_list=[8]
ec_2.q_is_square=True
##
ec_3 = m_case("ec_3", 3, [0]) #the only prime that's 0 mod 3 is 3
ec_3.dim=1
ec_3.m_list=[12]
ec_3.q_is_square=True;
##
ec_all= m_case("ec_all", 1, [0])
ec_all.dim=1
ec_all.m_list=[4]
ec_all.q_is_square=False

##############
###Dimension 2 cases for q not square

ss_not_2 = m_case("ss_not_2", 4, [1, 3])
ss_not_2.dim=2
ss_not_2.m_list=[8]
ss_not_2.q_is_square=False
##
ss_all=m_case("ss_all", 1,[0])
ss_all.dim=2
ss_all.m_list=[2,6]
ss_all.q_is_square=False
##
ss_not_3 = m_case("ss_not_3", 3, [1, 2])
ss_not_3.dim=2
ss_not_3.m_list=[12]
ss_not_3.q_is_square=False
##
ss_5 = m_case("ss_5", 5, [0])
ss_5.dim=2
ss_5.m_list=[10]
ss_5.q_is_square=False
##
ss_2 = m_case("ss_2", 2, [0])
ss_2.dim=2
ss_2.m_list=[24]
ss_2.q_is_square=False

################
dim1_cases=[ec_2, ec_3, ec_all]
dim2_cases=[ss_not_2, ss_all, ss_not_3, ss_5, ss_2];

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

#############
###OUTPUT
####
#{4, 8, 12, 20, 24}
