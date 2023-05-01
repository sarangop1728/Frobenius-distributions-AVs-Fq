load('qPolyData.sage')

'''
Use this script to generate data, using the functions from qPolyData.sage.
'''


# SAGE histograms.
# ______________________________________________________________________________
from sage.plot.histogram import Histogram

def str_moments(seq):
    return ', '.join([str(m) for m in seq])

def q_histogram(data, n, g):
    list_data = list(data)
    buckets = 10*2^(n-1)
    title = '$a_1$-distribution for ' + label
    subtitle = '\n\n' + '( 10^' + str(n) + ' data points and ' + str(buckets) + ' buckets )'
    str_mom = str_moments(moments(list_data))
    H = histogram(list_data, bins=buckets, density=True, range=[-2*g, 2*g], ticks=[[],[]],
                      color = 'hotpink', axes = False, edgecolor='hotpink', title=title,
                      frame = False, axes_labels=['Moments: '+ str_mom + subtitle, ''],
                      fontsize=18, axes_labels_size = 0.6, ymax=1/g)
    P = plot(1/(4*g), (x, -2*g, 2*g), color='gray', thickness = 2)
    return H+P
    

def histogram_from_label(label, n, extension = '.pdf'):
    '''
    label: LMFDB label of isogeny class.
    n: integer >= 2. 10^n is the number of data points.
    '''
    
    # Create qPolyData.
    poly = qPolyClass(poly_from_label(label))
    g = poly.g
    d = 2*g
    q = poly.q

    data = a1_sequence(poly, 10^n)
    figure_name = 'a1_' + label + '_10^' + str(n) + extension 
    path =  './stats/d=' + str(d) + '/q=' + str(q) + '/' + figure_name
    q_histogram(data, n, g).save(path) 
     
    

# Generate histograms.
# ______________________________________________________________________________

all_paper_labels = ['1.2.ab', '2.2.ab_a', '2.5.a_ab', '2.25.ac_bz', '2.2.ab_b', '2.2.a_ad', '2.2.ab_ab', '2.3.ac_c', '2.2.ad_f', '2.2.ac_f', '2.2.a_d', '2.7.af_s', '2.5.ag_s', '2.7.aj_bi', '3.2.a_a_ac', '3.2.a_a_ad', '3.2.a_a_af', '3.2.ab_b_b', '3.2.ab_f_ad', '3.2.ad_f_ah', '3.2.ad_j_an', '3.2.ae_j_ap', '3.2.ae_k_ar', '3.3.ad_d_ac', '3.3.af_p_abg', '3.3.af_r_abi', '3.5.ak_bv_afc', '3.7.ao_di_alk', '3.8.ag_bk_aea', '3.8.ai_bk_aeq']

paper_labels = [ '3.2.a_a_ac', '3.2.a_a_ad', '3.2.a_a_af', '3.2.ab_b_b', '3.2.ab_f_ad', '3.2.ad_f_ah', '3.2.ad_j_an', '3.2.ae_j_ap', '3.2.ae_k_ar', '3.3.ad_d_ac', '3.3.af_p_abg', '3.3.af_r_abi', '3.5.ak_bv_afc', '3.7.ao_di_alk', '3.8.ag_bk_aea', '3.8.ai_bk_aeq']

labels =  ['3.2.a_a_ad', '3.2.a_a_af', '3.2.ab_b_b', '3.2.ab_f_ad', '3.2.ad_f_ah', '3.2.ad_j_an', '3.2.ae_j_ap', '3.2.ae_k_ar', '3.3.ad_d_ac', '3.3.af_p_abg', '3.3.af_r_abi', '3.5.ak_bv_afc', '3.7.ao_di_alk', '3.8.ag_bk_aea', '3.8.ai_bk_aeq']

for label in all_paper_labels:
    for n in [2,3,4,5,6]:
        histogram_from_label(label, n, extension = '.pdf')
        print('histogram ' + label + '10^' + str(n) + ' has been generated...')

# to run this, clone the q-moments directory form github
# 1. open the terminal and go to the q-moments directory that you just cloned
# 2. run the command: sage generate_data.sage
