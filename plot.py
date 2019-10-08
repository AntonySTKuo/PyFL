import os, sys, fnmatch
dirname = os.path.dirname(os.path.abspath(__file__))
sys.path.append(dirname)
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result
fileList = find('*.npy', dirname)
fileList = [x.replace('.npy','') for x in fileList]


import numpy as np
import matplotlib.pylab as plt
# print("\n###########################################################")
# print("5. Plot results")

for file in fileList:
	gen_fit_plot = np.load(file+'.npy')
	n = 10
	for i in range(n):
	    plt.plot(gen_fit_plot[:,i])
	plt.show()


