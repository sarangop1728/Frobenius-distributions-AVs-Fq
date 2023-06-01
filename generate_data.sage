load('qPolyData.sage')

'''
Use this script to generate data, using the functions from qPolyData.sage.
'''

# SAGE histograms.
# ______________________________________________________________________________
from sage.plot.histogram import Histogram

def str_moments(seq):
    return ', '.join([str(m) for m in seq])

def q_histogram(data, n, g, label):
    list_data = list(data)
    buckets = 4^n
    title = '$a_1$-distribution for ' + label
    subtitle = '\n\n' + '(' + str(16^n) + ' data points and ' + str(buckets) + ' buckets )'
    str_mom = str_moments(moments(list_data))
    H = histogram(list_data, bins=buckets, density=True, range=[-2*g, 2*g], ticks=[[],[]],
                      color = 'hotpink', axes = False, edgecolor='hotpink', title=title,
                      frame = False, axes_labels=['Moments: '+ str_mom + subtitle, ''],
                      fontsize=18, axes_labels_size = 0.6, ymax=1/g)
    P = plot(1/(4*g), (x, -2*g, 2*g), color='gray', thickness = 2)
    return H+P
    

def histogram_from_label(label, n, paper_hist = True, extension = '.pdf'):
    '''
    label: LMFDB label of isogeny class.
    n: integer >= 2. 10^n is the number of data points.
    '''
    
    # Create qPolyData.
    poly = qPolyClass(poly_from_label(label))
    g = poly.g
    d = 2*g
    q = poly.q
    
    if paper_hist:
        directory = './stats/paper_histograms/'
    else:
        directory = './stats/g=' + str(g) + '/q=' + str(q)

    # Directory
    if not os.getcwd()[-9:] == "q-moments":
        print("Oh you duffer, you should change to the q-moments directory!")
    else:
        if not os.path.exists(directory):
            os.makedirs(directory)
    
    data = a1_sequence(poly, 16^n)
    figure_name = 'a1_' + label + '_16^' + str(n) + extension 
    path =  directory + '/' + figure_name
    q_histogram(data, n, g, label).save(path) 
     
    

# Generate histograms.
# ______________________________________________________________________________

all_paper_labels = ['1.2.ab', '2.2.ab_a', '2.5.a_ab', '2.25.ac_bz', '2.2.ab_b', '2.2.a_ad', '2.2.ab_ab', '2.3.ac_c', '2.2.ad_f', '2.2.ac_f', '2.2.a_d', '2.7.af_s', '2.5.ag_s', '2.7.aj_bi', '3.2.a_a_ac', '3.2.a_a_ad', '3.2.a_a_af', '3.2.ab_b_b', '3.2.ab_f_ad', '3.2.ad_f_ah', '3.2.ad_j_an', '3.2.ae_j_ap', '3.2.ae_k_ar', '3.3.ad_d_ac', '3.3.af_p_abg', '3.3.af_r_abi', '3.5.ak_bv_afc', '3.7.ao_di_alk', '3.8.ag_bk_aea', '3.8.ai_bk_aeq']

paper_labels = [ '3.2.a_a_ac', '3.2.a_a_ad', '3.2.a_a_af', '3.2.ab_b_b', '3.2.ab_f_ad', '3.2.ad_f_ah', '3.2.ad_j_an', '3.2.ae_j_ap', '3.2.ae_k_ar', '3.3.ad_d_ac', '3.3.af_p_abg', '3.3.af_r_abi', '3.5.ak_bv_afc', '3.7.ao_di_alk', '3.8.ag_bk_aea', '3.8.ai_bk_aeq']

labels =  ['2.5.a_ab', '2.25.ac_bz']
    

def paper_histograms(labels, exponents=[2,3,4,5]):
    for label in labels:
        for n in exponents:
            histogram_from_label(label, n, paper_hist = True, extension = '.pdf')
            print('histogram ' + label + '_16^' + str(n) + ' has been generated...')
    print('That is all, folks!')


def all_histograms(g,q, exponents=[2,3,4,5], extension='.pdf'):
    d = 2*g
    poly_list = poly_class_list(d,q)
    ext = extension
    for label in [poly.label for poly in poly_list]:
        for n in exponents:
            histogram_from_label(label, n, paper_hist = False, extension= ext)
            print('histogram ' + label + '_16^' + str(n) + ' has been generated...')
    print('That is all, folks!')
        
# paper_histograms(labels)
# all_histograms(2,2,exponents=[5],extension='.png')
# labels_to_txt(2,2)

for n in [2,3,4,5,6]:
    for label in labels:
        histogram_from_label(label,n,paper_hist=False,extension='.png')
        print('histogram ' + label + '_16^' + str(n) + ' was calculated!')
        
    

# to run this, clone the q-moments directory form github
# 1. open the terminal and go to the q-moments directory that you just cloned
# 2. run the command: sage generate_data.sage
