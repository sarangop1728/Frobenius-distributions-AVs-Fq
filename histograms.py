import matplotlib.pyplot as plt
import numpy as np
import csv
import os
import imageio

"""
Copy and paste the output of the moments() function here.
Coming soon: automatization of this caveman approach.
"""

''' HISTOGRAM '''

# read data

# a0_data = './stats/d=4/q=2/a0_4_2_10^3.csv'
def read_sequence_data(file_name):
    data = []
    with open(file_name, 'r') as F:
        reader = csv.reader(F, quoting=csv.QUOTE_NONNUMERIC)
        for line in reader:
            data.append(line)
    return data

# polys = './stats/d=4/q=2/polys_4_2.txt'
def read_poly_data(file_name):
    with open(file_name) as F:
        data = F.read().splitlines()
    return data
        

# buiild all a_0 histograms (10^n data points)  of all q-weil polys of degree d
def all_histograms(d,q,n):
    file_name = "./stats/d=" + str(d) + "/q=" + str(q) + '/a0_' + str(d) + '_' + str(q) + '_10^' + str(n) + '.csv'
    title = 'a_0 dist. for ' + str(q) + '-weil polynomials of degree ' + str(d) + ', 10^' + str(n) + ' data points.'
    figure_name = 'a0_' + str(d) + '_' + str(q) + '_10^' + str(n) + '.png'
    data = read_sequence_data(file_name)
    r = len(data)
    bins = np.arange(-1,1,0.02)
    cols = 7
    rows = int(np.ceil(r/cols))
    fig, axs = plt.subplots(rows, cols, sharex = True, sharey = False, tight_layout = True)
    fig.suptitle(title)
    for i in range(r):
        axs[int(np.floor(i/cols)), i%cols].hist(data[i], bins, density = True)
        axs[int(np.floor(i/cols)), i%cols].set_yticks([])
        plt.xticks([-1,0,1])
    path =  "./stats/d=" + str(d) + "/q=" + str(q) + '/' + figure_name
    plt.savefig(path, dpi=300)


# build GIF
def gif_all_histograms(d,q):
    
    # build list with images to cicle through
    images = []
    file_name =  "./stats/d=" + str(d) + "/q=" + str(q) + '/a0_' + str(d) + '_' + str(q) 
    for n in [2,3,4,5]:
        images.append(file_name + '_10^' + str(n) + '.png')
        images.append(file_name + '_10^' + str(n) + '.png')
        images.append(file_name + '_10^' + str(n) + '.png')
        images.append(file_name + '_10^' + str(n) + '.png')
        images.append(file_name + '_10^' + str(n) + '.png')
        images.append(file_name + '_10^' + str(n) + '.png')
        if (n == 5): # pause on the last frame
            images.append(file_name + '_10^' + str(n) + '.png')
            images.append(file_name + '_10^' + str(n) + '.png')
            images.append(file_name + '_10^' + str(n) + '.png')
            images.append(file_name + '_10^' + str(n) + '.png')
            images.append(file_name + '_10^' + str(n) + '.png')
            images.append(file_name + '_10^' + str(n) + '.png')
            images.append(file_name + '_10^' + str(n) + '.png')
            images.append(file_name + '_10^' + str(n) + '.png')
            images.append(file_name + '_10^' + str(n) + '.png')
            images.append(file_name + '_10^' + str(n) + '.png')

    # build gif
    gif_name = file_name + '.gif'
    with imageio.get_writer(gif_name, mode='I') as writer:
        for filename in images:
            image = imageio.imread(filename)
            writer.append_data(image)
