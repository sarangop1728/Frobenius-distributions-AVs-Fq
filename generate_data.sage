load('qPolyData.sage')

'''
Use this script to generate data, using the functions from qPolyData.sage.
'''

# Generate files.
# ______________________________________________________________________________


labels = ['1.2.ab']

for label in labels:
    for n in [2,3,4,5,6]:
        a1_label_csv(label,n)