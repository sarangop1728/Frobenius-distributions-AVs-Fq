load("qPolyData.sage")
# -----------------------------------------------------------
# Example 6.4.2
# m = 12
# link: https://www.lmfdb.org/Variety/Abelian/Fq/3/3/ac_e_ad
# -----------------------------------------------------------
# LMFDB label:
label = "3.3.ac_e_ad"
# Frobenius polynomial
P = poly_from_label(label)
P_class = qPolyClass(P)
q = P_class.q
K.<a> = NumberField(P)
# Check that sqrt(q) is not in K
sqrt_q_in_K = len((x^2-q).roots(K)) > 0
# Number field K = Q(a)
K.<a> = NumberField(P)
QuadraticSubfields = K.subfields(2)
# There is only one quadratic subfield B in K
len(QuadraticSubfields)
B, B_into_K, _ = QuadraticSubfields[0]
# Tell SAGE to consider K as a relative number field over B
K_over_B.<aa,bb> = K.relativize(B_into_K)
u1_squared = aa^2/q
# As expected, the relative norm of u1^2 is a primitive 6-th root of 1
nrm = u1_squared.relative_norm()
print([nrm^n for n in [1,2,3,4,5,6]])
