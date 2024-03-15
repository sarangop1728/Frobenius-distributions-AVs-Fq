import matplotlib.pyplot as plt
import numpy as np
import csv
import os
import imageio as iio
import sys


# Read data.
# ______________________________________________________________________________


def read_labels(g,q):
    file_name = './stats/g=' +str(g) + '/q=' + str(q) + '/' +str(g) +'_' + str(q) + '_labels.txt'
    with open(file_name) as F:
        data = F.read().splitlines()
    return data


# GIFs (n=2,3,4,5,6)
# ______________________________________________________________________________

def gif_histogram(label, exponents = [2,3,4,5], extension = '.png'):
    # Build list with images to cicle through
    images = []
    # Get d and q
    g = int(label[0])
    q = int(label[2])
    file_path =  './stats/g=' + str(g) + '/q=' + str(q) + '/a1_' + label
    for n in exponents:
        for i in range(10): # Repeat same image 10 times.
            file_name = file_path  + '_16^' + str(n) + extension
            images.append(file_name)
        if (n == exponents[-1]): # Pause on the last frame.
            file_name = file_path + '_16^' + str(n) + extension
            images.append(file_name)
            for i in range(30): # Repeat same image 30 times.
                images.append(file_name)
    # Build gif
    gif_name = file_path + '_16^' + str(exponents[-1]) + '.gif'
    # iio.mimsave(gif_name, images, fps=55)
    with iio.get_writer(gif_name, mode='I') as writer:
        for filename in images:
            image = iio.imread(filename)
            writer.append_data(image)

def gif_all_single_histograms(g,q):
    labels = read_labels(g,q)
    for label in labels:
        gif_histogram(label)
        print('gif of ' + label + ' has been generated...')

# DELETE GARBAGE HISTOGRAMS (n=2,3,4)
# ______________________________________________________________________________
def delete_histograms(g,q,exponents=[2,3,4,5],extension='.png'):
    # read poly data
    labels = read_labels(g,q)
    # loop over all files and delete them
    for label in labels:
        directory =  './stats/g=' + str(g) + '/q=' + str(q) + '/'
        for n in exponents:
            file_name = 'a1_' + label + '_16^' + str(n) + extension
            if os.path.exists(directory + file_name):
                os.remove(directory + file_name)
                print(file_name + ' removed!')
            else:
                print('The file ' + file_name + ' does not exist!')


# GENERATE DATA
# ______________________________________________________________________________
# g = 2
# q = 2
# gif_all_single_histograms(g,q)
# delete_histograms(g,q)

