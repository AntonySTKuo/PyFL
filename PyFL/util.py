import numpy as np 
from numpy.random import multinomial

def neighbors(seq):
    '''
    Generate a list of 1-mutation neighbors
    i.e. neighbors('AAA') --> ['UAA', 'CAA', 'GAA', 'AUA', 'ACA', 'AGA', 'AAU', 'AAC', 'AAG']
    '''
    base = ['A','C','G','U']
    neighbors = [seq[:i]+ch+seq[i+1:] for i in range(len(seq)) if seq[i] != 'N' for ch in base if ch != seq[i]]
    return neighbors


def seq2index(seq):
    v = 0
    for idx, char in enumerate(seq):
        k = 4**(len(seq)-idx-1)
        if char == "A":
            v += k * 0
        elif char == "C":
            v += k * 1
        elif char == "G":
            v += k * 2
        elif char == "U":
            v += k * 3
        else:
            print("error")
    return v


def to_matrix_index(seq):
    return [seq2index(s) for s in neighbors(seq)]



def weighted_variance(samples, weight):
    
    # Calculate weighted average
    average = numpy.average(values, weights=weights)
    # Calculate (biased) weighted variance
    variance = numpy.average((values-average)**2, weights=weights)
    return (average, variance)



    return 


def plot_clonal_interference(matrix):
    """
    matrix is a 2D numpy array
    """
    num_variant, time_length = matrix.shape 
    
    x = np.arange(0,time_length)
    #x = np.arange(0,10)
    #c = np.array([np.random.dirichlet(np.ones(10)*1000.,size=1).reshape(-1) for i in range(10)])
    for i in range(num_variant-1):
        matrix[i+1,:] = matrix[i+1,:] + matrix[i,:]
    plt.plot(matrix.T, lw = 0.1)
    for i in range(num_variant-1):
        plt.fill_between(x, matrix.T[:,i], matrix.T[:,i+1], alpha='0.8')
    plt.savefig("test.jpg")


def constant_sum_array_generator(N = 100000, total_length = 262144):
    """
    generate an array with length = 'total_length' and the sum of the elements = 'N'
    """
    a = multinomial(N, [1/total_length] * total_length)
    return a


