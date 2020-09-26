# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 15:16:02 2018

@author: MateusMarcuzzo

Some function prototypes related to Blind Source Separation with ICA 
Like: creation of the mixingMatrix
Permutation Matrices and so on

TO-DO
    create an inverse mod P, it will be used to check if the sep_matrix is
    A^-1 with permutated rows and scale adjustment
    
    11/11/2019
    The above was already done somewhere...

"""

import numpy as np
from itertools import product
from numpy.linalg import det
from sympy import Matrix
from numpy.random import random
from numpy.random import shuffle
from collections import Counter,defaultdict
from operator import itemgetter

import time




def round_det(the_matrix):
    """
    gives the rounded determinant of a matrix.
    Necessary, because python creates wrong determinants given an integer
    elements matrix.
    """
    return np.round(det(the_matrix))

def round_det_modP(the_matrix,prime=2):
    """
    It does round_det procedure and then 
    return the result mod P
    """
    return round_det(the_matrix) % prime

## new implementation! Based on Yeredor and D.G Silva code.
## It is faster and guarantee a full-rank matrix.
    
def time_of_generateMixingMatrix(n_sources,prime = 2):
    """
    This function is used for time evaluation of the generation of a mixing
    matrix
    
    It returns the total elapsed time
    """
    start_time = time.perf_counter()
    
    generateMixingMatrix(n_sources,prime)
    
    return  time.perf_counter() - start_time 

def time_of_generating_inverting(n_sources,prime=2): 
    """
    This function is used for time evaluation of finding the inverse matrix
    mod P of sympy in python
    
    It returns the total elapsed time
    """
    start_time = time.perf_counter()
    
    the_mix_matrix = generateMixingMatrix(n_sources,prime)
    Matrix_form = Matrix(the_mix_matrix)
    the_inv = Matrix_form.inv_mod(prime)
    
    the_inv = np.array(the_inv).astype(np.int64)
    
    identity_we_expect = the_mix_matrix.dot(the_inv) % prime
    assert( (identity_we_expect == np.eye(n_sources)).all() )
    
    print("n_sources: ",n_sources)
    print("prime: ",prime)
    print("generate mix matrix and inverse in: ")
    total_time = time.perf_counter()-start_time
    print(total_time," seconds")
    
    return total_time    
     
        

def generateMixingMatrix(n_sources, prime = 2):
    """
    It returns a random mixing matrix, given the n_sources and which P(prime)
    is given. Such that the mixing matrix is mod P
    """
    the_matrix = np.random.choice(prime,size=(n_sources,n_sources))
     
    the_matrix += np.diag(np.diag(the_matrix) == 0)
     
    lower_part = np.tril(the_matrix)
    upper_part = np.eye(n_sources) + np.triu(the_matrix,1)
     
    the_matrix = lower_part.dot(upper_part)
    the_matrix = the_matrix % prime
    the_matrix = np.array(the_matrix,dtype=np.int64)

#    print(the_matrix)
#    print(det(the_matrix))
#    print(round_det(the_matrix))
#    assert( round_det_modP(the_matrix,prime) != 0 )
    
    # The above assertion may fail for big matrices. since it may overflow 
    # or underflow float parts and give incorrect rounding values for giant 
    # floating numbers.
    # Since we have no direct det_in_GF version of it (where every operation
    # is mod P), we just could do mod P at the end, which gives this problem.
    
    # Since the algorithm is mathematically correct, we may comment the above 
    # line 
    
    
    return the_matrix

def numberOfUniqueHitsInSeparation(mix_matrix_x_sep):
    """ 
    Gives the number of separated sources of a separation problem.
    It works in the A*W matrix.
    
    it skips the step of checking if its a diagonal * permutation
    Matrix.
    
    OBS:(17/03/2019)
    The accurate method of hits is the rank (in modulo Prime) of the matrix
    The methods in numpy gives the rank in Real numbers, not Modulo P.
    (edit (11/03/2019)):
        it may not. since many results in row echelon form o mix_x_sep
        gives the identity, which is full rank.
        example:
            GF(5), mix x sep:
                0 1 0 0 0
                3 0 0 0 0
                1 0 0 3 0
                0 0 4 0 0
                0 0 0 0 1
                
                gives the identity in reduced row echelon form.
                OBS: some fail cases may give identity in row echelon form.
                
    
    This method is inspired in Yeredor's code and have shown some possible
    inconsistent results before for when not all the hits happen
    
    """
    
    return (((np.asarray(mix_matrix_x_sep) > 0).sum(axis=0)) == 1).sum()

## numpy.linalg import matrix_rank


#
# ## Not being used...
#    
#def generatePermutationMatrix(n_sources):
#    """
#    Generates a random Permutation Matrix given n_sources
#    A Permutation Matrix is the identity with its rows permutated
#    """
#    
#    my_identity = np.identity(n_sources)
#    
#    return np.random.permutation(my_identity)


def generateOccurrenceTensorShape(prime,n_sources):
    """
    Creates the Tuple with n_sources dimensions, with each one prime possible 
    values
    It is used in the building of a OcurrenceTensor, which later will become 
    a Probability Tensor
    """
    
    a_list = []
    for each in range(n_sources):
        a_list.append(prime)
        
    return tuple(a_list)


# Python does not use function overloading. So...we did this:
def generateOccurrenceTensor(prime,n_sources):
    """
    Returns a zero valued Occurrence Tensor given the prime and n_sources
    """
    shape = generateOccurrenceTensorShape(prime,n_sources)
    return generateOccurrenceTensorFrom(shape)

def generateOccurrenceTensorFrom(tensor_shape):
    """
    Returns a zero valued Occurrence Tensor given the tensor_shape.
    The format of the tensor is a tuple.
    """
    return np.full(tensor_shape,0)

# Changes the tensor in-place
def incrementOccurrenceOnTensor(tensor, position_tuple):
    """
    When reading an observation, increment its value on corresponding
    Occurrence Tensor
    """
    a_tuple = tuple(position_tuple)
    tensor[a_tuple] += 1
    
def generateRandomProbabilityMassFunction(pmf_size):
    """
    Returns an np.array which positions designate symbols probablities.
    These positions have random values.
    
    These random values are the probabilities of symbols 0, until (pmf_size-1)
    included. Like a Probability Mass Function
    """
    the_pmf = []
    remaining = 1.0
    
    if(pmf_size == 1):
        return np.array([1.0])
    
    for position in range(pmf_size):
        
        if(position == 0):
            prob = random()
            the_pmf.append(prob)
            remaining = 1.0 - prob
            
            
        elif(position < (pmf_size-1)):
            prob = remaining * random()
            the_pmf.append(prob)
            remaining = remaining - prob
            
            
        else:
            the_pmf.append(remaining)

            
    the_pmf = np.array(the_pmf)
    shuffle(the_pmf)

    return the_pmf

def generateUniformDistribution(pmf_size):
    scalar = 1.0/pmf_size
    the_pmf = [scalar]*pmf_size
    
    the_pmf = np.array(the_pmf)
    
    return the_pmf
        
        
    
def generateRandomTupleFrom(pmf_list):
    """ 
    Creates a random Tuple from given pmf_list
    
    It infers by the size of the pmf length, the order of the field
    So it knows the possible outcomes.
    
    All pmf must be the same size. So pmf_list is a 2d array.
    """
    
    
    possible_outcomes = np.arange(len(pmf_list[0]))
    
    the_tuple = []
    
    for pmf in pmf_list:
        the_tuple.append(np.random.choice(possible_outcomes,p=pmf))
        
    the_tuple = tuple(the_tuple)
    
    return the_tuple

def OccurToProbTensor(occur_tensor):
    """
    Transforms an Occurrence Tensor in a Probability Tensor
    """
    return occur_tensor/np.sum(occur_tensor)


def createIndexVectorSize(k_length, prime):
    """
    The (n-1)-th index is the (n-1) number in P-ary form:
    example:
        
        3 is when prime=2 , k_length = 2:
        (2,1) 
        which is : 2* 2 + 1*1 = 3.
        
        5 when prime=3 k_length = 3:
        (0,1,2) = 1* 3 + 2*1 = 5.
        
    ####
    
    This functions create a K_length index-vector which is a 1d representation
    of a Tensor. It organizes in an orderly manner the elements of the Tensor.
    
    It depends on k_length and prime, such as the above example.
            
            
    """
    the_index_vector = list(np.full(prime**k_length,0))
    
    for index,value in enumerate(product(range(prime),repeat=k_length)):
        list_value = list(value)
    
        the_index_vector[index] = tuple(list_value)
    
    return the_index_vector


def createReversedIndexVectorSize(k_length, prime):
    """
    Creates the index_vector like original article from Yeredor (2011)
    The first (0 index) is always an all-zero tuple.
    The following is (1,0,...,0) tuple vector
    The last is (P-1,P-1,...,P-1)
    
    So instead of growing the P-ary tuple from the right to the left 
    It grows from the left-most algarism to the rightmost
    
    ###
    
    In reversed order of the traditional way, like the above function
    """
    the_index_vector = list(np.full(prime**k_length,0))
    
    for index,value in enumerate(product(range(prime),repeat=k_length)):
        list_value = list(value)
    

        #This is just to get the order that the original article uses
        list_value.reverse()
        the_index_vector[index] = tuple(list_value)
    
    return the_index_vector

#########################################################
#########################################################
#########################################################
#################    ENTROPY PART      ##################
#########################################################
#########################################################
#########################################################
#########################################################

def onlineEntropyLauradouxEstimator_Initialization(the_samples,r,log_base=2):
    """
    Its the initialization step for the online entropy estimator of
    Cedric Lauradoux et. al. It's on the page 9 of the original article.
    
    It returns the empirical entropies. 
    
    The mean of these entropies gives the true estimative of
    entropy of the sources. 
    
    Each value in the empirical entropies is like an
    "instantaneous velocity" of entropy.
    """
    assert (r < len(the_samples))
    
    log_2 = np.log(2)
    
    # r==1 is a special case. It's descripted in the paper
    if(r == 1):
        adjust = 2
    else:
        adjust = 1/log_2
        
    entropy_empirical = np.full(len(the_samples)-r,0,dtype=np.float64)


    for i in np.arange(r,len(the_samples)):
        # don't have a name at all
        w = 0
        
        # don't have a name too
        j = 1
        
        while(j <= r and the_samples[i] != the_samples[i-j]):
            w = w + 1/j
            j = j + 1
          
         # H_P(X) = H_2(X)/log(P)_base(2)
        entropy_empirical[i-r] =  w *  adjust /np.log2(log_base)
        
    return entropy_empirical


def onlineEntropyLauradouxEstimator_LatestSample(the_samples,r,log_base):
    """
    This function is like the "start" (Initialization) one.
    But we don't want to recalculate all
    the empirical entropies everytime. 
    
    So we just estimate a new entropy to add into the entropies array 
    
    with the latest added sample in "the_samples"
    
    obs:it seems it DOES NOT change entropy_empirical in-place.
     so it returns a single new instantaneous entropy.
        
    """
    assert (r < len(the_samples))
    
    log_2 = np.log(2)
    
    # r==1 is a special case. It's descripted in the paper
    if(r == 1):
        adjust = 2
    else:
        adjust = 1/log_2
     
        
    i = -1
    w = 0
    
    # don't have a name too
    j = 1
    
    while(j <= r and the_samples[i] != the_samples[i-j]):
        w = w + 1/j
        j = j + 1
        
    
    the_new_entropy_empirical =  w *  adjust/np.log2(log_base)
    
    
    return the_new_entropy_empirical



def entropyRandomVariable(randomvariable_probabilities,log_base = 2):
    """     
    Given a list with the probabilities, calculates its entropy.
    Doesn't check if it is a true Probability Mass Function (PMF)
    It will transform in a PMF based on occurrences and use the Maximum
    Likelihood Estimation for the probabilities.
    
    """
    
    entropy=0
    denominator = np.asarray(randomvariable_probabilities,dtype='float64')\
    .sum()
    for value in randomvariable_probabilities:
        if(value > 0.0):
            entropy += (value/denominator) *  np.log2(value/denominator)
            
        
    return -entropy/np.log2(log_base)



def entropyFromArray(the_array,log_base=2):
    """
    Given an array of occurrences, and the possible_events of the occurrences,
    calculate it's shannon entropy of predefined log_base
    
    04/11/2019
    we may change such that we don't need "possible_events"
    done.
    
    Also works for non-numerical values
    """
    
    event_counter = Counter(the_array)
    
    pmf = []
    for each in event_counter.values():
        pmf.append(each)
        
    pmf = np.array(pmf)
    pmf = pmf/pmf.sum()

    
    return entropyRandomVariable(pmf,log_base)


def jointEntropyRandomVariable(probability_tensor, log_base = 2):
    """
    Given a probability tensor, calculates its entropy.
    Since a joint-entropy is an entropy of a random variable too,
    we can interpret the joint random variable as a singular one. We just need
    To format the initial input in something that resembles traditional
    entropy and then calculate its entropy as usual.
    """
    
    # here we just make the tensor 1-D so we treat as a normal RV 
    randomvariable_probabilities = np.ravel(probability_tensor)
    randomvariable_probabilities /= randomvariable_probabilities.sum()
    return entropyRandomVariable(randomvariable_probabilities, log_base)


def entropiesSum(probability_tensor, log_base = 2):
    """
    Given the Probability Tensor, gives the Sum of the individual entropies.
    Supose you have a joint RV XWYZ. (4 sources) in GF(2).
    So the tensor looks like:
    array(
        [[[[1, 1],
         [1, 1]],

        [[1, 1],
         [1, 1]]],


       [[[1, 1],
         [1, 1]],

        [[1, 1],
         [1, 1]]]]
    
    
    We can see that each line is uniform (base 2 entropy = 1) 
    so the sum must be 8 
    
    OBSERVATION:
    (17/03/2019) - não deveria ser 4?
    visto que são 4 v.as aleatórias com dist uniforme cada uma
    
    The sum of entropies should be:
        H(X) + H(W) + H(Y) + H(Z)
        
        
        is it really correct? (17/03/2019)
        No, I don't think so (04/11/2019)
        to extract the marginal entropies, you should something like
        the new function in list_of_marg...
        
        I think it's correct now
    """
    marg_probs_list = list_of_marg_probs_in(probability_tensor)
    
    entropies_sum = 0
    
    for pmf in marg_probs_list:
        entropies_sum += entropyRandomVariable(pmf,log_base)   
    
    return entropies_sum

def totalCorrelation(probability_tensor, log_base = 2):
    """
    Gives the total Correlation of a Random Variable Vector 
    (A Vector of Random Variables)
    Given by:
        C(y) = Sum H(y_i) - H(y)
    
    Also known as generalized Mutual Information, or simply Mutual Information
    
    Total Correlation name came from Satoshi Watanabe.
    
    totalCorrelation relies on entropiesSum, take care if it's not correct.
    """
    
    return entropiesSum(probability_tensor, log_base) -\
 jointEntropyRandomVariable(probability_tensor, log_base)

def tuplify(array_array):
    """
    04_11_2019
    We need this because list is a not hashable type.
    So we can't do totalCorr and such, we need to transform lists into tuples
    inside of the array_of_arrays
    """
    tuplified = []
    for line in array_array:
        tuplified.append(tuple(line))
        
    return tuplified

 
def totalCorr_FromArray(the_array, log_base=2):
    """
    04_11_2019
    I'll try to do here a TotalCorr directly from Array.
    04_11_2019
    Continuing...
    08_03_2020
    I don't know if it works.
    """
    the_array = tuplify(the_array)
    joint_entropy = entropyFromArray(the_array)
    
    marg_entropies_sum = 0
    for line in np.array(the_array).T:
        marg_entropies_sum+= entropyFromArray(line)
        
    return marg_entropies_sum - joint_entropy
    
    
 
def k_L_Divergence(pmf_1,pmf_2 = 'uniform', log_base= np.e):
    """
    Gives the Kullback-Leibler Divergence from pmf_1,pmf_2
    given:
        sum P(x)ln(P(x)/Q(x))
        x in X
        
        P here represents pmf_1, and Q pmf_2
        
        log_base can be used to switch base 'e' to any other valid base
    """
    
    if(pmf_2 == 'uniform'):
        pmf_2 = generateUniformDistribution(len(pmf_1))
        
    total = 0.0
    for i in range(len(pmf_1)):
        if(pmf_1[i] > 0.0):
            # This if uses the same scheme from entropy. Don't know
            # if it's valid for this case, btw.
            total += pmf_1[i] *  np.log(pmf_1[i] / pmf_2[i])
        
    return total/np.log(log_base)


def gini(array):
    """
    Calculate the Gini coefficient of an array-like type.
    
    OBS: entropy evaluation is still faster than gini's 
     calculation for small values (in python)
    Gini seems to scales better (also in python)
    
    for resemblance as entropy, use one_less_gini.
    
    # Internet:
    # https://github.com/oliviaguest/gini  
    # Had to do some adaption
    """
    
    # All values are treated equally, arrays must be 1d:
    temp_array = np.array(array).flatten()
    temp_array = np.array(temp_array,dtype=np.float64)
    
    if np.amin(array) < 0:
        # Values cannot be negative:
        temp_array -= np.amin(temp_array)
        
    # Values cannot be 0:
    
    ## Adds the most minimum positive number after zero of numpy
    #temp_array += np.finfo(np.float64).eps
    
    # Values must be sorted:
    temp_array = np.sort(temp_array)
    
    # Number of array elements:
    n = temp_array.shape[0]
    
    # Index per array element:
    index = np.arange(1,n+1)
    
    
    
    # Gini coefficient:
    return ((np.sum((2 * index - n  - 1) * temp_array)) / (n * np.sum(temp_array)))



def one_less_gini(array):
    """ 
    Which is really close or similar to entropy when log_base = len(pmf_list)
    """
    return 1 - gini(array)

def tuple_rev(a_tuple):
    """
    reverses a tuple
    """
    return tuple(reversed(a_tuple))

def tuple_given(number,n_sources,prime):
    """
    returns a tuple given which number it is
    and n_sources and prime
    
    example:
        number = 10
        n_sources = 3
        prime = 3
        
        The number base = prime,
        the maximum slots for fullfilling is n_sources, and the number in this
        form is:
            (1,0,1)
            
            3**2 + 3**0 = 10, as expected
    """
    remainders = []
    quotient = number
    
    while(True):
        quotient, r = (divmod(quotient,prime))
        remainders.append(r)
        
        if quotient == 0:
            break
        
    while(len(remainders) < n_sources):
        remainders.append(0)
        
        
    remainders.reverse()
    return tuple(remainders)
    


# 31/10/2019
# https://en.wikipedia.org/wiki/Tsallis_entropy
def Tsallis_Entropy(pmf,q=2):
    """
    Gives the Tsallis entropy of a probability mass function
    """
    
    pmf = np.array(pmf)
    if(q==1):
        return entropyRandomVariable(pmf)
    
    
    the_sum = 1- (pmf**q).sum()
    return (1/(q-1))*the_sum

def extract_marg_prob_of_value_in_dim_(value,dimension,the_tensor):
    """
    A support method for list_of_marg_probs_in
    """
    assert(dimension<= the_tensor.ndim)
    index_like=dimension-1
    assert(value < the_tensor.shape[index_like])
    
    marg_prob_of_value  = the_tensor[(slice(None), )*index_like + (value, )]
    marg_prob_of_value = marg_prob_of_value.sum()
    
    return marg_prob_of_value

def list_of_marg_probs_in(the_tensor):
    """
    Gives a list of the marginal probabilities in a tensor
    """
    all_marg_probs = []

    for dimension in range(the_tensor.ndim):
        marg_probs = []
        for values in range(the_tensor.shape[dimension]):
            marg_probs.append(extract_marg_prob_of_value_in_dim_(values,dimension+1,the_tensor))
        all_marg_probs.append(marg_probs)
        
    return all_marg_probs
    

def show_marg_probs_in(prob_tensor):
    """
    Show the marginal probababilities in a tensor
    """
    for dimension in range(pmf_conj.ndim):
        for values in range(pmf_conj.shape[dimension]):
            print('dimension {}, value {}'.format(dimension+1,values))
            print(extract_marg_prob_of_value_in_dim_(values,dimension+1,pmf_conj))
    
def generate_random_pmf_conj_in(alphabet_size,n_dim):  
    """ 
    Generates a random probability mass function, which is joint,given an
    alphabet size and number of dimensions
    """
    pmf_conj = generateRandomProbabilityMassFunction(alphabet_size**n_dim)
    
    shape = [alphabet_size]*n_dim
    
    pmf_conj = pmf_conj.reshape(tuple(shape))
    
    return pmf_conj

def generate_random_pmf_from_independent_sources(alphabet_size,n_dim):
    """
    Given a bunch of pmfs of fixed size,
    you merge than such that they generate a joint pmf with 
    a totalCorrelation 0
    """
    pmfs_list = []
    if n_dim < 2:
        raise ValueError('n_dim must be 2 or more'
                         )
    for _ in range(n_dim):
        pmfs_list.append(generateRandomProbabilityMassFunction(alphabet_size))
    
    pmf_conj = np.kron(pmfs_list[0],pmfs_list[1])
    
    for i in range(2,n_dim):
        pmf_conj = np.kron(pmf_conj,pmfs_list[i])
        
    the_shape = generateOccurrenceTensorShape(alphabet_size,n_dim)
    return pmf_conj.reshape(the_shape)

    
def sort_pmf_and_index(the_pmf):
    """
    Sorts a 1-D pmf and returns the indices permutation
    """
    # inspired in:
    # https://stackoverflow.com/a/7851186/1644727
    the_pmf = list(the_pmf)
    indices, pmf_sorted = zip(*sorted(enumerate(the_pmf), key=itemgetter(1)))
    
    return pmf_sorted,indices

def permutate_pmf_given_indices(the_pmf,indices):
    """
    Given a pmf: e.g.
    [0.2,0.1,0.7]
    
    permutate given indices:
        
    [0,2,1]
    
    which would result in:
        
    [0.2,0.7,0.1]
        
    
    """
    # inspired in
    # https://stackoverflow.com/a/6459056/1644727
    
    the_pmf = list(the_pmf)
    permuted_pmf = [the_pmf[i] for i in indices]
    
    return permuted_pmf
    
def order_permutation_algorithm(pmf_conj):
    """
    The binary order permutation algorithm from Painsky's
    """
    line_pmf_conj = pmf_conj.ravel()
    line_pmf_conj,indices = sort_pmf_and_index(line_pmf_conj)
   
    
    ndim = pmf_conj.ndim
    new_pmf_conj = np.zeros(shape=shape)
    for index,p in enumerate(product(range(pmf_conj.shape[0]),repeat=ndim)):
        new_pmf_conj[p] = line_pmf_conj[index]
        
    return new_pmf_conj
    
if __name__ == '__main__':
    ## Testing something 31/10/2019 about number of permutations
    import itertools as it
    from time import sleep
    from collections import Counter,defaultdict
    from scipy.special import factorial
    # Uma PMF conjunta de duas v.as x1,x2 tais que elas resultam nisso, independentes.
    # Por inspeção (caderno de ideias), verificamos que existem 3 combinações das
    # conjuntas, que geram três ordenações possíveis das marginais em ordem.
    
    # se tivéssemos agora, 3 v.as binárias...
    
    #pmf1 = 0.7 e 0.3
    #pmf2 = 0.6 e 0.4
    # quantas combinações possíveis das marginais em ordem teríamos?
    #pmf_conj = np.array([0.42,0.28,0.18,0.12])
    #Algo mais louco...descomente o de cima e comente o debaixo
    pmf_conj = generate_random_pmf_conj_in(3,2)
    
    shape = pmf_conj.shape
    
    combinations_dict = defaultdict(int)
    
    ## Rapaz, isso aqui nunca vai dar certo pra P**K maior que
    ## 8
    
    ## ou seja, os casos admissíveis são:
    ## 2^2, 2^3, 3^2
    ## E SÓ
    

    alphabet_size=5
    dimension=3
    pmf_conj = generate_random_pmf_from_independent_sources(alphabet_size,dimension)
    print(pmf_conj)
    line_pmf_conj = pmf_conj.ravel()
    shuffle(line_pmf_conj)
    pmf_conj = line_pmf_conj.reshape(pmf_conj.shape)
    print(pmf_conj)
    
    #pmf_conj = generate_random_pmf_conj_in(3,3)
    line_pmf_conj = pmf_conj.reshape(-1)
    shape = pmf_conj.shape
    fac_val = factorial(len(line_pmf_conj),exact=True)
    print("factorial of number_of_joint_probs is ",fac_val)
    
    print(f'alphabet size {alphabet_size}')
    print(f'dimensions : {dimension}')
    min_total_corr = totalCorrelation(pmf_conj)
    max_actual_total_corr = min_total_corr
    print('actual max total_corr',max_actual_total_corr)
    idx_max = 0 
    idx_min = 0
    the_max_conj = pmf_conj
    the_min_conj = []
#    for idx,_ in enumerate(it.permutations(line_pmf_conj)):
#        
#        #Just to check when cardinality is too big
#        print(idx,fac_val)
#        _ = np.array(_)
#        
#        pmf_conj_temp = _.reshape(shape)
#        total_corr = totalCorrelation(pmf_conj_temp)
#        
#        if  (total_corr > max_actual_total_corr):
#            max_actual_total_corr = total_corr
#            idx_max = idx
#            the_max_conj = pmf_conj_temp
#            
#        if (total_corr < min_total_corr):
#            min_total_corr = total_corr
#            idx_min = idx
#            the_min_conj = pmf_conj_temp
#            
#            
#        concat = np.concatenate(list_of_marg_probs_in(pmf_conj_temp))
#    
#        the_comb= np.sort(concat)
#        
#        # the sort is ascending order
#        the_comb= np.flip(the_comb)
#        
#        #dicts do not accept numpy arrays
#        the_comb = str(the_comb)
#       # print(pmf_conj,"\n-----------------",the_comb)
#        
#        combinations_dict[the_comb]+=1
#    
#        #sleep(1)
#    print('max actual total_corr and the index',max_actual_total_corr,idx_max)
#    print(the_max_conj) 
#    print('min total corr and the index',min_total_corr,idx_min)
#    print(the_min_conj)   
#        
#    #print(combinations_dict.items())
#    #print(len(combinations_dict))
#    print(Counter(combinations_dict.values()))
    
    
     ## Test code in terminal to see plots showing were o_gini(pmf) > entropy(pmf,log_base = len(pmf))
     ## 11/03/2019
     ## 22/03/2019 - It does not hold on Octave/Matlab, entropy is still faster.
     ## Some other tests in entropy_estimators.py
    #prime = 5
    #pmf=generateRandomProbabilityMassFunction(prime)
    #while (one_less_gini(pmf) <= entropyRandomVariable(pmf,len(pmf))):
    #    pmf=generateRandomProbabilityMassFunction(prime)
    #    plt.plot(pmf)
    #    
    #plt.show()
    #plt.plot(pmf)
    #plt.show()
    #print(pmf)
    #print(max(pmf)/min(pmf))
    #print(max(pmf)/np.mean(pmf))
    
    
    # not working for 2^3
    print('now with order permutation')
    print(f'pmf conj shape {pmf_conj.shape}')
    print(totalCorrelation(order_permutation_algorithm(pmf_conj)))
     
     
     
     
     ##Using pandas to store the values (11/03/2019)
    import pandas as pd
    #my_df = pd.DataFrame(columns=['prime','o_gini','entropy','k_l_div'])
    #df_pmf = pd.DataFrame(columns=['pmf'])
    
    #### Some experiments in 30/10/2019############################################
    
    ## THese tries weren't good... 31/10/2019
    #def mixed_sort(an_array):
    #    temp_array = list(np.sort(an_array))
    #    
    #    final_array = []
    #    
    #    while len(final_array) != len(an_array):
    #
    #            final_array.append(temp_array.pop())
    #            
    #            if(len(temp_array) != 0):
    #                final_array.append(temp_array.pop(0))
    #                
    #    return np.array(final_array)
    #
    #
    #def menos_mais(an_array):
    #    temp_array = list(an_array)
    #    
    #    for idx,_ in enumerate(temp_array):
    #        temp_array[idx] *= (-1)**idx 
    #     
    #    return np.abs(np.array(temp_array).sum())
    #
    #def menos_mais_mult(an_array):
    #    temp_array = list(an_array)
    #    
    #    for idx,_ in enumerate(temp_array):
    #        temp_array[idx] **= (-1)**idx 
    #     
    #    return np.abs(np.prod(np.array(temp_array)))
    #    
    
    #def nova_f_custo(an_array):
    #    temp_array = mixed_sort(an_array)
    #    
    #    return menos_mais(temp_array)
    #    
    #    
    #def nova_f_custo2(an_array):
    #    temp_array = mixed_sort(an_array)
    #    
    #    return  menos_mais_mult(temp_array)
    #    
    #def nova_f_custo3(a_prob_array):
    #    temp_array = np.array(a_prob_array)
    #    return (np.e**(temp_array.prod()))
    
    
    #def implies_relation(x1,x2):
    #    return not(x1) or x2
    #
    #########################################################################
    #prime = 7
    #
    #n_trials = 10**5
    #rows_list = []
    #for _ in range(n_trials):
    #
    #    pmf = generateRandomProbabilityMassFunction(prime)
    #    pmf2 = generateRandomProbabilityMassFunction(prime)
    #    o_gini1 = one_less_gini(pmf)
    #    o_gini2 = one_less_gini(pmf2)
    #    
    #    entropy1 = entropyRandomVariable(pmf,2)
    #    entropy2 = entropyRandomVariable(pmf2,2)
    #    
    #
    #    temp_dict = {
    #    'prime':prime,
    #    'o_gini_1':o_gini1,
    #    'o_gini_2':o_gini2,
    #    'entropy1':entropy1,
    #    'entropy2':entropy2,
    #    'o_gini1_less_o_gini2':o_gini1 < o_gini2,
    #    'entropy1_less_entropy2':entropy1 < entropy2,
    #    'entropy_less_equivalence_o_gini_less': (o_gini1 < o_gini2) == (entropy1 < entropy2),
    #
    #    'o_gini_implies_entropy': implies_relation(o_gini1< o_gini2,entropy1 < entropy2) 
    #    }
    #    
    #    
    #    rows_list.append(temp_dict.copy())
    #
    #my_df = pd.DataFrame(rows_list)
    ##df_pmf = df_pmf.append([{'pmf':pmf}],ignore_index=True)
    #    
    #    
    #print(my_df.head())
    #print(my_df.mean())
    ## Bussabi e Casellla
    #print('----')
    #
    #my_df_erro_padrao_media = my_df.std()/np.sqrt(my_df.count())
    #print("erros padrões da média de df")
    #print(my_df_erro_padrao_media)
    # 
    #
    #
    #
    #before = time.time()
    #print('entropy times')
    #for _ in range(n_trials):
    #
    #    pmf = generateRandomProbabilityMassFunction(prime)
    #    entropyRandomVariable(pmf)
    #    
    #print("elapsed ",time.time()-before," seconds")
    #print("mean: {} seconds".format((time.time()-before)/n_trials))
    #
    #
    #before = time.time()
    #print('o_gini times')
    #for _ in range(n_trials):
    #
    #    pmf = generateRandomProbabilityMassFunction(prime)
    #    one_less_gini(pmf)
    #    
    #print("elapsed ",time.time()-before," seconds")
    #print("mean: {} seconds".format((time.time()-before)/n_trials))
    ######### end 30/10/2019###################################################
    
            
    #prime = 5
    #n_trials = 50000
    #before = time.time()
    #for _ in range(n_trials):
    #    pmf=generateRandomProbabilityMassFunction(prime)
    #    my_df=my_df.append([{
    #
    #'prime':prime,
    #'pmf':pmf,
    #'o_gini':one_less_gini(pmf),
    #'entropy':entropyRandomVariable(pmf,2),
    #'k_l_div':k_L_Divergence(pmf)}],ignore_index=True)
    #print("elapsed ",time.time()-before," seconds")
    # 
     
     
    ## for performance consideration (11/03/2019)
    