load("qPolyData.sage")
# -----------------------------------------------------------
# Example 6.4.2 (E)
# link: https://www.lmfdb.org/Variety/Abelian/Fq/3/4/ad_j_au
# -----------------------------------------------------------
# LMFDB label:
label = "3.4.ad_j_au"
# Frobenius polynomial
P = poly_from_label(label)
P_class = qPolyClass(P)
q = P_class.q
sqrt_q = ZZ(sqrt(q))
# Number field K = Q(a)
K.<a> = NumberField(P)
QuadraticSubfields = K.subfields(2)
# There is only one quadratic subfield B in K
len(QuadraticSubfields)
B, B_into_K, _ = QuadraticSubfields[0]
# Tell SAGE to consider K as a relative number field over B
K_over_B.<aa,bb> = K.relativize(B_into_K)
u1 = aa/sqrt_q
# As expected, u1 is a primitive 6-th root of 1
nrm = u1.relative_norm()
print([nrm, nrm^2, nrm^3, nrm^4, nrm^5, nrm^6])