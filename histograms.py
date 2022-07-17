import matplotlib.pyplot as plt
import numpy as np
import csv

"""
Copy and paste the output of the moments() function here.
Coming soon: automatization of this caveman approach.
"""



''' HISTOGRAM '''

# read moments data
data = []
with open('weil_4_2.csv', 'r') as F:
    reader = csv.reader(F, quoting=csv.QUOTE_NONNUMERIC)
    for line in reader:
        data.append(line)

# plot histograms
r = len(data)
bins =  np.arange(-1, 1, 0.01)
cols = 7
rows = int(np.ceil(r/cols))
fig, axs = plt.subplots(rows, cols, sharex = True, sharey = False, tight_layout = True)
fig.suptitle('Distribution of NTF for 2-weil polynomials of degree 4')
for i in range(r):
    axs[int(np.floor(i/cols)), i%cols].hist(data[i], bins, density = True)
    axs[int(np.floor(i/cols)), i%cols].set_yticks([])
plt.xticks([-1,0,1])
plt.show()
