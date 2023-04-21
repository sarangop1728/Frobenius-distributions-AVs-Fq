load('qPolyData.sage')

'''
Use this script to generate data, using the functions from qPolyData.sage.
'''

# Generate files.
# ______________________________________________________________________________

all_paper_labels = ['1.2.ab', '2.2.ab_a', '2.5.a_ab', '2.25.ac_bz', '2.2.ab_b', '2.2.a_ad', '2.2.ab_ab', '2.3.ac_c', '2.2.ad_f', '2.2.ac_f', '2.2.a_d', '2.7.af_s', '2.5.ag_s', '2.7.aj_bi', '3.2.a_a_ac', '3.2.a_a_ad', '3.2.a_a_af', '3.2.ab_b_b', '3.2.ab_f_ad', '3.2.ad_f_ah', '3.2.ad_j_an', '3.2.ae_j_ap', '3.2.ae_k_ar', '3.3.ad_d_ac', '3.3.af_p_abg', '3.3.af_r_abi', '3.5.ak_bv_afc', '3.7.ao_di_alk', '3.8.ag_bk_aea', '3.8.ai_bk_aeq']

labels =  ['3.8.ag_bk_aea']

for label in all_paper_labels:
    for n in [6]:
        a1_label_csv(label,n)
