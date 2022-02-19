import matplotlib.pyplot as plt

"""
Copy and paste the output of the moments() function here.
Coming soon: automatization of this caveman approach.
"""

data = 


# subdivision of interval [-1,1] into 100 bins
bins = [-1.0 + i/100 for i in range(0,201)] 
plt.hist(data, bins)
plt.show()
