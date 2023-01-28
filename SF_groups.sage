load('qPolyData.sage')

'''
In this script we calculate the Serre-Frobenius groups of isogeny classes of
abelian varieties in dimension less than or equal to 3, using the results of 
our paper.
'''
#_____________________________________________________________________________________
# Supersingular elliptic curves
#_____________________________________________________________________________________

def SF_ss_simple_1(karg):
    '''
    Input: LMFDB label or Frobenius poly of a (supersingular) elliptic curve E over Fq.
    Output: m; the order of the component subgroup of SF(E).
    '''
    if type(karg) == str:
        label = karg
        poly = poly_from_label(label)
    else:
        assert karg.is_weil_polynomial(), f'{poly} is not a Weil polynomial'
        poly = karg

    Poly = qPolyClass(poly)
    g = Poly.g
    p = Poly.p
    a1 = poly[1]
    assert g == 1 and a1 % p == 0, f'{label} is not a supersingular elliptic curve.'
    q = Poly.q
    r = Poly.r # Recall that q = p^r.
    q_is_square = (r % 2 == 0)
    if q_is_square:
        sqrt_q = ZZ(sqrt(q))
        if a1 == 2*sqrt_q or a1 == -2*sqrt_q:
            return 1
        elif a1 == -sqrt_q:
            assert p % 3 != 1, f'Something is wrong.. p = {p} is 1 mod 3'
            return 3
        elif a1 == sqrt_q:
            assert p % 3 != 1, f'Something is wrong.. p = {p} is 1 mod 3'
            return 6
        else:
            assert a1 == 0 and p % 4 != 1 , f'Something is wrong.. either p = {p} is 1 mod 3 or {a1} is not zero.'
            return 4
    else:
        if a1 == 0:
            return 4
        elif p == 2:
            assert a1 == ZZ(sqrt(p*q)) or a1 == -ZZ(sqrt(p*q)), 'Something is wrong.. {a1} should be \pm {ZZ(sqrt(p*q))}.'
            return 8
        else:
            assert p == 3
            assert a1 == ZZ(sqrt(p*q)) or a1 == -ZZ(sqrt(p*q)), 'Something is wrong.. {a1} should be \pm {ZZ(sqrt(p*q))}.'
            return 12

# Dictionary with torision order set M = M(p,d) data.
SS_1_M_data_d_even = {
    1:['(-)',{1}],
    2:['(p not 1 mod 3)',{1,3,6}],
    3:['(p not 1 mod 4)',{1,4}]
    }

SS_1_M_data_d_odd = {
    1:['(-)',{4}],
    2:['(p=2)',{4,8}],
    3:['(p=3)',{4,12}]
    }
    
#_____________________________________________________________________________________

# Supersingular non-simple surfaces
#_____________________________________________________________________________________

def SF_ss_nonsimple_2(karg):
    '''
    Input: LMFDB label or Frobenius poly of a nonsimple SS surface S over Fq.
    Output: m; the order of the component subgroup of SF(S).
    '''
    if type(karg) == str:
        label = karg
        poly = poly_from_label(label)
    else:
        assert karg.is_weil_polynomial(), f'{poly} is not a Weil polynomial.'
        poly = karg

    F = poly.factor()

    if len(F) == 2: # If poly = (poly1)x(poly2) is not a square.
        poly1 = F[0][0]
        poly2 = F[1][0]
        return lcm(SF_ss_simple_1(poly1), SF_ss_simple_1(poly1))
    elif len(F) == 1 and F[0][1] == 4: # If poly = h^4.
        poly1 = F[0][0]^2
        return SF_ss_simple_1(poly1)
    else: # If poly=(poly1)^2 is a square.
        poly1 = F[0][0]
        return SF_ss_simple_1(poly1)

#_____________________________________________________________________________________

# Supersingular simple surfaces
#_____________________________________________________________________________________

def SF_ss_simple_2(karg):
    '''
    Input: LMFDB label or Frobenius poly of a simple SS surface S over Fq.
    Output: m; the order of the component subgroup of SF(S).
    '''
    if type(karg) == str:
        label = karg
        poly = poly_from_label(label)
    else:
        assert karg.is_weil_polynomial(), f'{poly} is not a Weil polynomial.'
        poly = karg

    Poly = qPolyClass(poly)
    g = Poly.g
    q = Poly.q
    p = Poly.p
    a1 = poly[1]
    a2 = poly[2]
    r = Poly.r # Recall that q = p^r.
    assert g == 2 and Poly.angle_rank == 0, f'{label} is not a supersingular surface.'
    q_is_square = (r % 2 == 0)
    if q_is_square:
        sqrt_q = ZZ(sqrt(q))
        if [a1,a2] == [0,0]:
            assert p % 8 != 1, f'Something is wrong.. p = {p} is 1 mod 8.'
            return 8
        elif [a1,a2] == [0,-q]:
            assert p % 12 != 1, f'Something is wrong.. p = {p} is 1 mod 12.'
            return 12
        elif [a1,a2] == [sqrt_q,q]:
            assert p % 5 != 1, f'Something is wrong.. p = {p} is 1 mod 5.'
            return 5
        elif [a1,a2] == [-sqrt_q,q]:
            assert p % 5 != 1, f'Something is wrong.. p = {p} is 1 mod 5.'
            return 10
        elif [a1,a2] == [0,2*q]:
            assert p % 4 == 1, f'Something is wrong.. p = {p} is not 1 mod 4.'
            return 4
        elif [a1,a2] == [2*sqrt_q,3*q]:
            assert p % 3 == 1, f'Something is wrong.. p = {p} is not 1 mod 3.'
            return 3
        else:
            assert [a1,a2] == [-2*sqrt_q,3*q], f'Something is wrong..'
            assert p % 3 == 1, f'Something is wrong.. p = {p} is not 1 mod 3.'
            return 6
    else:
        if [a1, a2] == [0,0]:
            assert p != 2, f'Something is wrong.. p = 2.'
            return 8
        elif [a1,a2] == [0,q]:
            return 6
        elif [a1,a2] == [0,-q]:
            assert p != 3, f'Something is wrong.. p = 3.'
            return 12
        elif [a1,a2] == [0,-2*q]:
            return 2
        elif p == 5:
            assert [a1,a2] == [ZZ(sqrt(p*q)),3*q] or [a1,a2] == [-ZZ(sqrt(p*q)),3*q], f'Something is wrong..'
            return 10
        else:
            assert p == 2, f'Something is wrong.. p should be 2.'
            assert [a1,a2] == [ZZ(sqrt(p*q)),3*q] or [a1,a2] == [-ZZ(sqrt(p*q)),3*q], f'Something is wrong..'
            return 10
        

def create_nonsimple_SS_2_M_data_d_even():
    '''
    Input: No input.
    Output: Dictionary with all torsion order sets.
    '''
    data_dict = dict()
    counter = 1
    for i in SS_1_M_data_d_even:
        for j in SS_1_M_data_d_even:
            category = ''
            if i == j: # When S ~ E1xE2 of same (p,d)-type.
                category = SS_1_M_data_d_even[i][0]
                M = SS_1_M_data_d_even[i][1]
                data_dict[counter] = [category,M]
                counter += 1
            elif i < j and i > 1: # When S ~ E1xE2 of different (p,d)-type.
                category = SS_1_M_data_d_even[i][0] + 'x' + SS_1_M_data_d_even[j][0]
                I = SS_1_M_data_d_even[i][1]
                J = SS_1_M_data_d_even[j][1]
                M = {lcm(x) for x in cartesian_product([I,J])}                
                data_dict[counter] = [category,M]
                counter += 1
    return data_dict


def create_nonsimple_SS_2_M_data_d_odd():
    '''
    Input: No input.
    Output: Dictionary with all torsion order sets.
    '''
    data_dict = dict()
    counter = 1
    for i in SS_1_M_data_d_odd:
        category = ''
        category = SS_1_M_data_d_odd[i][0]
        M = SS_1_M_data_d_odd[i][1]
        data_dict[counter] = [category,M]
        counter += 1
    return data_dict

# non-simple surface data
SS_2_nonsimple_M_data_d_even = create_nonsimple_SS_2_M_data_d_even()
SS_2_nonsimple_M_data_d_odd = create_nonsimple_SS_2_M_data_d_odd()

# ----

def create_simple_SS_2_M_data_d_even():
    '''
    Input: No input.
    Output: Dictionary with all torsion order sets.
    '''
    data_dict = dict()
    categories = ['(p not 1 mod 8)', '(p not 1 mod 12)', '(p not 1 mod 5)', '(p 1 mod 4)', '(p 1 mod 3)']
    torsion_sets = [{8}, {12}, {5,10}, {4}, {3,6}]
    for i in range(len(categories)):
        data_dict[i+1] = [categories[i], torsion_sets[i]]

    return data_dict

def create_simple_SS_2_M_data_d_odd():
    '''
    Input: No input.
    Output: Dictionary with all torsion order sets.
    '''
    data_dict = dict()
    categories = ['(p not 2)', '(-)', '(p not 3)', '(p = 5)', '(p = 2)']
    torsion_sets = [{8}, {2,6}, {12}, {10}, {24}]
    for i in range(len(categories)):
        data_dict[i+1] = [categories[i], torsion_sets[i]]

    return data_dict



# simple surface data
SS_2_simple_M_data_d_even = create_simple_SS_2_M_data_d_even()
SS_2_simple_M_data_d_odd = create_simple_SS_2_M_data_d_odd()
#_____________________________________________________________________________________

# Supersingular simple threefolds
#_____________________________________________________________________________________

def create_simple_SS_3_M_data_d_even():
    '''
    Input: No input.
    Output: Dictionary with all torsion order sets.
    '''
    data_dict = dict()
    categories = ['(7 nmid p^3-1)', '(p not 1 mod 3)']
    torsion_sets = [{7,14}, {9,18}]
    for i in range(len(categories)):
        data_dict[i+1] = [categories[i], torsion_sets[i]]
    
    return data_dict

def create_simple_SS_3_M_data_d_odd():
    '''
    Input: No input.
    Output: Dictionary with all torsion order sets.
    '''
    data_dict = dict()
    categories = ['(p=7)', '(p=3)']
    torsion_sets = [{28}, {36}]
    for i in range(len(categories)):
        data_dict[i+1] = [categories[i], torsion_sets[i]]
    
    return data_dict

# simple threefold data
SS_3_simple_M_data_d_even = create_simple_SS_3_M_data_d_even()
SS_3_simple_M_data_d_odd = create_simple_SS_3_M_data_d_odd()

#_____________________________________________________________________________________

# Supersingular non-simple threefolds
# (In this case X ~ SxE for supersingular S and E.)
#_____________________________________________________________________________________


def create_nonsimple_SS_3_M_data_d_even():
    '''
    Input: No input.
    Output: Dictionary with all torsion order sets.
    '''
    data_dict = dict()
    counter = 1
    # S simple.
    for i in SS_2_simple_M_data_d_even:
        for j in SS_1_M_data_d_even:
            category = SS_2_simple_M_data_d_even[i][0] + 'x' + SS_1_M_data_d_even[j][0]
            I = SS_2_simple_M_data_d_even[i][1]
            J = SS_1_M_data_d_even[j][1]
            M = {lcm(x) for x in cartesian_product([I,J])}
            data_dict[counter] = [category, M]
            counter += 1
    return data_dict





# non-simple threefold data
SS_3_nonsimple_M_data_d_even = create_nonsimple_SS_3_M_data_d_even()
