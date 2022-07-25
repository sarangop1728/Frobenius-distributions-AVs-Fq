import matplotlib.pyplot as plt
import numpy as np
import csv
import os
import imageio


#_____________________________________________________________
# read data
#_____________________________________________________________


'''
read_sequence_data

'''

def read_sequence_data(d,q,n):
    data = []
    file_name = "./stats/d=" + str(d) + "/q=" + str(q) + '/a0_' + str(d) + '_' + str(q) + '_10^' + str(n) + '.csv'
    with open(file_name, 'r') as F:
        reader = csv.reader(F, quoting=csv.QUOTE_NONNUMERIC)
        for line in reader:
            data.append(line)
    return data
#_____________________________________________________________

def read_poly_data(d,q):
    file_name = './stats/d=' +str(d) + '/q=' + str(q) + '/polys_' +str(d) +'_' + str(q) + '.txt'
    with open(file_name) as F:
        data = F.read().splitlines()
    return data

#_____________________________________________________________

def read_moments_data(d,q):
    file_name = './stats/d=' +str(d) + '/q=' + str(q) + '/moments_' +str(d) +'_' + str(q) + '.txt'
    with open(file_name) as F:
        data = F.read().splitlines()
    return data
        
#_____________________________________________________________
# HISTOGRAMS
#_____________________________________________________________

'''
all_histograms: gif for all (d,q)-a0 histograms

*** looks bad when there are too many polynomials.

'''
def all_histograms(d,q,n):
    buckets = 2**(n-1)*10
    title = 'a0 dist. for ' + str(q) + '-weil polynomials of degree ' + str(d) + ', \n 10^' + str(n) + ' data points and ' + str(buckets) + ' buckets.'
    figure_name = 'a0_' + str(d) + '_' + str(q) + '_10^' + str(n) + '.png'
    data = read_sequence_data(d,q,n)
    r = len(data)
    bins = np.arange(-1,1,2/buckets)
    cols = 7
    rows = int(np.ceil(r/cols))
    fig, axs = plt.subplots(rows, cols, sharex = True, sharey = False, tight_layout = True)
    fig.suptitle(title)
    for i in range(r):
        axs[int(np.floor(i/cols)), i%cols].hist(data[i], bins, density = True)
        axs[int(np.floor(i/cols)), i%cols].set_yticks([])
        plt.xticks([-1,0,1])
        axs[int(np.floor(i/cols)), i%cols].set_ylim([0,1.6])
    path =  "./stats/d=" + str(d) + "/q=" + str(q) + '/' + figure_name
    plt.savefig(path, dpi=300)


#_____________________________________________________________
'''
histogram: histogram for a given (d,q)-Weil polynomial

INPUT: i, index in the list of all (d,q)-Weil polynomials. n, corresponding to 10^n data points.

OUTPUT:

REQUIRES: files 
'''
# 'x^4 + 4*x^3 + 8*x^2 + 8*x + 4'
# './stats/d=4/q=2/a0_4_2_10^2.csv'
def histogram(polyname,n):

    # get d and q
    d = polyname[2] # the third letter is always the degree
    q = str(int(float(polyname[-1])**(2/float(d)))) # the last letter is always the prime power
    
    # number of buckets
    buckets = 2**(n-1)*10

    # get index
    polydata = read_poly_data(d,q)
    index = polydata.index(polyname)

    # a0 sequence and moments
    sequence = read_sequence_data(d,q,n)[index]
    moments = read_moments_data(d,q)[index]
    moments = moments.replace("'",'')
    
    # remove the '*' from name
    polyname = polyname.replace('*','')

    title = 'a0 dist. for ' + polyname
    subtitle = ' 10^' + str(n) + ' data points and ' + str(buckets) + ' buckets.'
    figure_name = 'a0_' + polyname.replace(' ','') + '_10^' + str(n) + '.png'
    
    
    bins = np.arange(-1,1,2/buckets)
    fig , ax = plt.subplots()
    fig.suptitle(title, fontsize=15)
    ax.set_title(subtitle, fontsize=10)
    ax.hist(sequence, bins, density = True)
    ax.set_yticks([])
    ax.set_xticks([-1,0,1])
    ax.set_ylim([0,1.6])
    ax.set_xlim([-1.1,1.1])
    ax.set_xlabel('moments: ' + moments, fontsize=10)
    path =  "./stats/d=" + str(d) + "/q=" + str(q) + '/' + figure_name
    plt.savefig(path, dpi=300)
    plt.close('all')

#_____________________________________________________________
def d_q_histograms(d,q,n):
    polydata = read_poly_data(d,q)
    for polyname in polydata:
        histogram(polyname,n)



#_____________________________________________________________
# GIFs
#_____________________________________________________________

'''
gif_all_histograms: gif for all (d,q)-a0 histograms

*** looks bad when there are too many polynomials.

'''
def gif_all_histograms(d,q):
    
    # build list with images to cicle through
    images = []
    file_name =  "./stats/d=" + str(d) + "/q=" + str(q) + '/a0_' + str(d) + '_' + str(q) 
    for n in [2,3,4,5]:
        for i in range(10): #repeat same image 10 times
            images.append(file_name + '_10^' + str(n) + '.png')
        if (n == 5): # pause on the last frame
            for i in range(20): #repeat same image 20 times
                images.append(file_name + '_10^' + str(n) + '.png')

    # build gif
    gif_name = file_name + '.gif'
    with imageio.get_writer(gif_name, mode='I') as writer:
        for filename in images:
            image = imageio.imread(filename)
            writer.append_data(image)

#_____________________________________________________________

'''
gif_histogram: .gif for one (d,q)-a0 histogram

'''
def gif_histogram(polyname):
    
    # build list with images to cicle through
    images = []

    # get d and q
    d = polyname[2] # the third letter is always the degree
    q = str(int(float(polyname[-1])**(2/float(d)))) # the last letter is always the prime power

    
    # remove the '*' from name and delete spaces
    polyname = polyname.replace('*','')
    polyname = polyname.replace(' ','')
    
    file_name =  "./stats/d=" + d + "/q=" + q + '/a0_' + polyname
    for n in [2,3,4,5]: 
        for i in range(10): #repeat same image 10 times
            images.append(file_name + '_10^' + str(n) + '.png')

        if (n == 5): # pause on the last frame
            for i in range(20): #repeat same image 20 times
                images.append(file_name + '_10^' + str(n) + '.png')

    # build gif
    gif_name = file_name + '.gif'
    with imageio.get_writer(gif_name, mode='I') as writer:
        for filename in images:
            image = imageio.imread(filename)
            writer.append_data(image)

#_____________________________________________________________

'''
gif_all_single_histograms:

'''
def gif_all_single_histograms(d,q):
    # read poly data
    polydata = read_poly_data(d,q)
    # loop over all files and delete them
    for polyname in polydata:
        gif_histogram(polyname)
    

#_____________________________________________________________
# DELETE GARBAGE HISTOGRAMS (n=2,3,4)
#_____________________________________________________________

def delete_histograms(d,q,n):
    # read poly data
    polydata = read_poly_data(d,q)
    # loop over all files and delete them
    for polyname in polydata:
        # remove the '*' from name and delete spaces
        polyname = polyname.replace('*','')
        polyname = polyname.replace(' ','')
        directory =  "./stats/d=" + str(d) + "/q=" + str(q) + '/'
        file_name = 'a0_' + polyname + '_10^' + str(n) + '.png'
        if os.path.exists(directory + file_name):
            os.remove(directory + file_name)
        else:
            print("The file " + file_name + " does not exist!")
        

