import matplotlib.pyplot as plt
import numpy as np
import csv
import os
import imageio
import sys
from statistics import mean

csv.field_size_limit(sys.maxsize) # csv files contain very large fields. 

# Read data.
# ______________________________________________________________________________

def read_sequence_data(d,q,n):
    '''Reads file: a1_d_q_10^n.csv '''

    data = {}

    file_name = './stats/d=' +str(d) + '/q=' + str(q) + '/a1_' + str(d) + '_' + str(q) + '_10^' + str(n) + '.csv'

    with open(file_name, 'r') as F:

        reader = csv.DictReader(F)

        for line in reader:

            temp_string = line['a_1'].lstrip('[').rstrip(']') # Removes [ and ].
            temp_list = temp_string.split(',') # Create list of x_i's as strings.
            a_1 = list(map(float, temp_list)) # Coerce the x_i from strings to floats.
            data[line['Label']] = a_1

    return data


def read_sequence_from_label(label,n):
    '''Reads file: a1_g.q.iso_10^n.csv '''
    d = 2*int(label.split(".")[0])
    q = int(label.split(".")[1])

    data = {}

    file_name = './stats/d=' +str(d) + '/q=' + str(q) + '/a1_' + label + '_10^' + str(n) + '.csv'

    with open(file_name, 'r') as F:

        reader = csv.DictReader(F)

        for line in reader:

            temp_string = line['a_1'].lstrip('[').rstrip(']') # Removes [ and ].
            temp_list = temp_string.split(',') # Create list of x_i's as strings.
            a_1 = list(map(float, temp_list)) # Coerce the x_i from strings to floats.
            data[line['Label']] = a_1

    return data

def read_labels(d,q):

    file_name = './stats/d=' +str(d) + '/q=' + str(q) + '/' +str(d) +'_' + str(q) + '_labels.txt'

    with open(file_name) as F:

        data = F.read().splitlines()

    return data

# Moments.
# ______________________________________________________________________________

def k_moment(sequence, k):

    l = len(sequence)

    if l > 0:
        if k == 0:
            return 1
        else:
            return abs(mean(s**k for s in sequence))
    else:
        return -10^10

def moments(d,q,N=5,M=10):

    a1 = read_sequence_data(d,q,N)

    values = {}

    for label in a1.keys():
        values[label] = ['{:.2f}'.format(k_moment(a1[label], k)) for k in range(M)]

    return {key : value for key,value in values.items()}

def moments_from_label(label,n,M=10):

    a1 = read_sequence_from_label(label,n)

    values = {}

    for label in a1.keys():
        values[label] = ['{:.2f}'.format(k_moment(a1[label], k)) for k in range(M)]

    return {key : value for key,value in values.items()}

# Histograms.
# ______________________________________________________________________________


def histogram(label,n):

    # Read d and q.
    d = 2*int(label.split(".")[0])
    q = int(label.split(".")[1])

    # Number of buckets.
    buckets = (2**(n-1))*10

    # Read data
    sequence = read_sequence_data(d,q,n)
    moment_data = moments(d,q,n)

    # Plot
    title = 'a1 distribution for ' + label
    subtitle = '( 10^' + str(n) + ' data points and ' + str(buckets) + ' buckets )'
    figure_name = 'a1_' + label + '_10^' + str(n) + '.png'
    bins = np.arange(-2,2,2/buckets)
    fig , ax = plt.subplots()
    fig.suptitle(title, fontsize=11)
    ax.set_title(subtitle, fontsize=9)
    ax.hist(sequence[label], bins, density = True, color = 'hotpink')
    ax.set_yticks([])
    ax.set_xticks([-2,0,2])
    ax.set_ylim([0,1.6])
    ax.set_xlim([-2.1,2.1])
    momts = 'Moments: ' + str(moment_data[label])
    momts_new = momts.replace('[','').replace(']','').replace("'",'').replace(',',' ')
    ax.set_xlabel(momts_new, fontsize=10)
    path =  './stats/d=' + str(d) + '/q=' + str(q) + '/' + figure_name
    plt.savefig(path, dpi=300)
    plt.close('all')


def histogram_from_label(label,n):

    # Read d and q.
    g = int(label.split(".")[0])
    d = 2*g
    q = int(label.split(".")[1])

    # Number of buckets.
    buckets = (2**(n-1))*10

    # Read data
    sequence = read_sequence_from_label(label,n)
    moment_data = moments_from_label(label,n)

    # Plot
    title = 'a1 distribution for ' + label
    subtitle = '( 10^' + str(n) + ' data points and ' + str(buckets) + ' buckets )'
    figure_name = 'a1_' + label + '_10^' + str(n) + '.pdf'
    bins = np.arange(-d,d,2*d/buckets)
    fig , ax = plt.subplots()
    fig.suptitle(title, fontsize=12)
    ax.set_title(subtitle, fontsize=9)
    ax.hist(sequence[label], bins, density = True, color = 'hotpink')
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylim([0,2/d]) # Vertical scaling.
    ax.set_xlim([-d,d])
    ax.axhline(y = 1/(2*d), color = 'lightgray', linestyle = '-') # Uniform distribution baseline.
    momts = 'Moments: ' + str(moment_data[label])
    momts_new = momts.replace('[','').replace(']','').replace("'",'').replace(',',' ')
    ax.set_xlabel(momts_new, fontsize=10)
    path =  './stats/d=' + str(d) + '/q=' + str(q) + '/' + figure_name
    plt.savefig(path, dpi=300)
    plt.close('all')

def d_q_histograms(d,q,n):
    '''Generates all (d,q)-histograms for a1-sequences of length 10^n.'''

    labels = read_labels(d,q)

    for label in labels:
        histogram(label,n)





# GIFs (n=2,3,4,5,6)
# ______________________________________________________________________________

def gif_histogram(label):

    # Build list with images to cicle through.
    images = []

    # Get d and q.
    d = 2*int(label[0])
    q = int(label[2])

    file_name =  './stats/d=' + str(d) + '/q=' + str(q) + '/a1_' + label
    for n in [2,3,4,5,6]:
        for i in range(10): # Repeat same image 10 times.
            images.append(file_name + '_10^' + str(n) + '.png')

        if (n == 6): # Pause on the last frame.
            for i in range(30): # Repeat same image 30 times.
                images.append(file_name + '_10^' + str(n) + '.png')

    # Build gif.
    gif_name = file_name + '.gif'
    with imageio.get_writer(gif_name, mode='I') as writer:
        for filename in images:
            image = imageio.imread(filename)
            writer.append_data(image)


def gif_all_single_histograms(d,q):

    labels = read_labels(d,q)

    # loop over all files and delete them
    for label in labels:
        gif_histogram(label)


# DELETE GARBAGE HISTOGRAMS (n=2,3,4)
# ______________________________________________________________________________
def delete_histograms(d,q,n):
    # read poly data
    labels = read_labels(d,q)
    # loop over all files and delete them
    for label in labels:
        directory =  './stats/d=' + str(d) + '/q=' + str(q) + '/'
        file_name = 'a1_' + label + '_10^' + str(n) + '.png'
        if os.path.exists(directory + file_name):
            os.remove(directory + file_name)
        else:
            print('The file ' + file_name + ' does not exist!')

'''            
for n in [2,3,4,5,6]:
    for q in [2,3]:
        delete_histograms(4,q,n)
'''

# GENERATE DATA
# ______________________________________________________________________________

all_paper_labels =  ['1.2.ab', '2.2.ab_a', '2.5.a_ab', '2.25.ac_bz', '2.2.ab_b', '2.2.a_ad', '2.2.ab_ab', '2.3.ac_c', '2.2.ad_f', '2.2.ac_f', '2.2.a_d', '2.7.af_s', '2.5.ag_s', '2.7.aj_bi', '3.2.a_a_ac', '3.2.a_a_ad', '3.2.a_a_af', '3.2.ab_b_b', '3.2.ab_f_ad', '3.2.ad_f_ah', '3.2.ad_j_an', '3.2.ae_j_ap', '3.2.ae_k_ar', '3.3.ad_d_ac', '3.3.af_p_abg', '3.3.af_r_abi', '3.5.ak_bv_afc', '3.7.ao_di_alk', '3.8.ag_bk_aea', '3.8.ai_bk_aeq']
            
labels = ['3.8.ag_bk_aea']

for label in labels:
    for n in [6]:
        histogram_from_label(label,n)
#    gif_histogram(label)
