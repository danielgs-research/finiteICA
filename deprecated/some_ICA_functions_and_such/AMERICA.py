# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 00:55:28 2018

@author: MateusMarcuzzo

Here we implement and provide methods and interface
for using the Algorithm AMERICA (Yeredor) for BSS in ICA.

Original Article:
    ICA over Galois Fields of Prime Order (Yeredor)

Done:
    Implementation of Online Entropy Estimator K-ICA algorithm, without
       update step
    
Implemented functions:
    generateObservations_UpdateTensor
    generateObservations_AppendToSamples
    executeAMERICAAlgorithm
    trial_AMERICA_EqualPMFs
    execute_AMERICA_online
    trial_AMERICA_EqualPMFs_OnlineEntropy
"""


import BSS_functions as bss

from BSS_functions import entropyRandomVariable
from BSS_functions \
import onlineEntropyLauradouxEstimator_Initialization as online_Init

from BSS_functions \
import onlineEntropyLauradouxEstimator_LatestSample as online_LatestSample

from BSS_functions \
import one_less_gini as o_gini


import numpy as np
from numpy.linalg import det
from numpy.fft import fftn, ifft

from Galois_Fields import getGFarray,getINTarray, GF, \
row_echelon_form

import time


class k_ica_algorithm():
    def __init__(self,min_criteria):
        self.min_criteria
        
    
class trial():
    """
    This class is meant to abstract the K-ICA trials 
     with different algorithms and samples generation
     
    It may generate metrics and different results
    """
    
    def __init__(self,algorithm,sample_generator):
        """
        Initialize a trial setup from two functions:
            -> an algorithm
            -> a sample_generator
        """
        self.algorithm = algorithm
        self.sample_generator = sample_generator
        
        #(n_samples,n_sources,src_0_pmf,mix_matrix_type = 'random'):
    def execute_trial(self,n_samples,pmf_list,mix_matrix,\
                      function_criteria="entropy",verbosity=True):
        """
        You must find before executing the trial its parameters:
            n_samples, pmf_list and mix_matrix.
            
            n_samples is a positive integer
            pmf_list is dimension-2 matrix with each line with pmf
             respective to each component.
             
            to generate a pmf_list, you could use
            -> bss.generateRandomProbabilityMassFunction
            -> bss.generateUniformDistribution
            -> put individually the pmf_list
            
            and finally, the mix_matrix.
            to generate a mix_matrix, you could use
            -> bss.generateMixingMatrix
            -> put a custom mixing_matrix
            
            OBS: the mix_matrix dimensions are correlated to n_sources.
             one should may wonder why we don't use n_sources. 
                first: DRY principle.
                 we use n_sources to generate a mix_matrix before inputing it
                 to the function.
                two: It may be extracted from the mix_matrix dimensions or 
                 the tensor dimension.
             
             
        
        """
        
        start_time = time.time()
        finite_field = len(pmf_list[0])
        prime = finite_field
        
        n_sources = mix_matrix.shape[0]
        
        
        mix_matrix = getGFarray(mix_matrix,finite_field)
        occTensor = bss.generateOccurrenceTensor(finite_field,n_sources)
        self.sample_generator(mix_matrix,occTensor,pmf_list,n_samples)
       
        
        sep_matrix = self.algorithm(occTensor,function_criteria)
        
        W_x_A = getINTarray(sep_matrix.dot(mix_matrix))
        
        n_unique_sep_sources = bss.numberOfUniqueHitsInSeparation(W_x_A)
        
        if(verbosity):
            print("Also: W*A")
            
            print(W_x_A)
            
            print("the number of unique hits:")
            
            print(n_unique_sep_sources)
            
            if(n_unique_sep_sources == n_sources):
                print("Fully separated!")
            else:
                print("Not Totally separated...")
            
            
            print("The mix matrix")
            print(getINTarray(mix_matrix))
            print("the_mixing_matrix in row_echelon_form")
            mixing_in_row_echelon = getINTarray(row_echelon_form(mix_matrix))
            print(mixing_in_row_echelon)
            
            
            print("The sep_matrix")
            print(sep_matrix)
            print("Must be full rank: det != 0")
            print("sep_matrix in row_echelon_form")
            sep_in_row_echelon = getINTarray(row_echelon_form(getGFarray(sep_matrix,prime)))
            print(sep_in_row_echelon)
            
            W_x_A_row_echelon = W_x_A.copy()
            W_x_A_row_echelon = getGFarray(W_x_A_row_echelon,prime)
            W_x_A_row_echelon = getINTarray(row_echelon_form(W_x_A_row_echelon))
            print("W_x_A in echelon form")
            print(W_x_A_row_echelon)
            
            
            
            assert(bss.round_det_modP(sep_matrix,finite_field) != 0)
            
            print(bss.round_det_modP(sep_matrix,finite_field))
        
        print("run n_samples = ",n_samples," in: ","--- %s seconds ---" % (time.time() - start_time))
    

    
        return n_unique_sep_sources == n_sources
    
    
###




############ STARTED ABOVE in 14/03/2019 ##############
############ TRYING TO USE OOP ##################

def generateObservations_UpdateTensor(mix_matrix,occ_tensor,pmf_list,n_samples):
    """
    Changes the occ_tensor in-place, modyfing the occurrences values it stores
    Given a pmf_list, which indicates the probability of a given tuple 
     to happen
     
    The mix_matrix must be in GF form. See the GaloisFields Package
    """
    for sample in range(n_samples):
        the_tuple = bss.generateRandomTupleFrom(pmf_list)
        sample_observed = (mix_matrix).dot(np.asarray(the_tuple))
        sample_observed = getINTarray(sample_observed)

        bss.incrementOccurrenceOnTensor(occ_tensor,sample_observed)
        
def generateObservations_AppendToSamples(mix_matrix,the_samples,pmf_list):
    """
    Given the Mixing Matrix, generates an observation and appends to the
    current observed samples_list.
    
    The mix_matrix must be in GF form
    """
    the_tuple = bss.generateRandomTupleFrom(pmf_list)
    sample_observed = (mix_matrix).dot(np.asarray(the_tuple))
    sample_observed = getINTarray(sample_observed)

    the_samples.append(tuple(sample_observed) )
    
def executeAMERICAAlgorithm(occur_tensor,function = "entropy"):
    """
    The description on original article is dubious and may direct
    the programmer to an error.

    We will describe it later (12/11/2018)
    """
    
    prob_tensor = bss.OccurToProbTensor(occur_tensor)
    
    prob_char_tensor = fftn(prob_tensor)


    n_sources = len(occur_tensor.shape)
    prime = occur_tensor.shape[0]



    index_vector = bss.createReversedIndexVectorSize(n_sources,prime)
    ## just made to get some improvements on processing
    index_vector_nparray = np.array(index_vector)
    
    ## See Table I from original article
    y_char_tensor = np.full((prime**n_sources,prime),1.0,dtype=np.complex128)
    y_empirical_prob = np.full((prime**n_sources,prime),0,dtype=np.float64)

    
    possible_n_values = range(prime**n_sources)
    
    ### THis part needs optimization.
    
    for n in possible_n_values:
        
#        y_char_tensor[n][0] = 1.0
        y_char_tensor[n,1] = prob_char_tensor[index_vector[n]].copy()

        if prime > 2:
            y_char_tensor[n,prime-1] = np.conj(prob_char_tensor[index_vector[n]])
            
        #When prime is 3, it does not execute the code below.
        for m in np.arange(2,(prime)/2,dtype=int):
            
            temporary_m_index = (m * index_vector_nparray % prime).astype(int)
            
            m_index_vector = []
    
            ## pra acessar como matlab faz, sem o índice completo:
            ## an_array.flat[0] <- primeira posição
            ## dá pra modificar elementos fazendo:
            ## an_array.flat[posicao] = novo_valor
            
            ## flat também pode acessar vários elementos:
            
            ## an_array.flat[[pos1,pos2,...]] = novo_valor
            for a_tuple in temporary_m_index:
                m_index_vector.append(tuple(a_tuple))
         
            y_char_tensor[n,m] = prob_char_tensor[m_index_vector[n]]
    
    
    
## The code above is greatly inspired in this section of Yeredor's Matlab
#  original code
##if P>2
#    qf(P,:)=conj(fPPy);
#    for m=2:P/2
#        mLex=mod(m*Lex,P);
#        qf(m+1,:)=fPPy(r*mLex+1);
#        qf(P+1-m,:)=conj(qf(m+1,:));
#    end
#end

    
    entropies = np.full((prime**n_sources),0,dtype=np.float64)

     
    #OBS: the ifft is reliable for little numbers, as probabilities.
    # it just gives problems when concerning about big numbers.
    # when doing the ifft, we may get little, but negative complex part values, which
    # gives us a problem when evaluating to the entropies calculation.
    # So we converted the whole number to float.


    # Here we collect the vT*X entropies, for each v.
    # Note the above notation is naive. T is the transpose, X would be the 
    # samples; we do not use it directly, because of FFT manipulation.
    # This type of showing appears in a later article after the original one.
    for n in possible_n_values:
        
        #this a.real.astype is necessary to not raise ComplexWarning.
        y_empirical_prob[n] = ifft(y_char_tensor[n]).real.astype(np.float64)
        
        if(function == "entropy"):
            entropies[n] = \
            entropyRandomVariable(y_empirical_prob[n],log_base=prime)
        elif(function == "o_gini"):
            entropies[n] = \
            o_gini(y_empirical_prob[n])


    
    # n must be different from zero. So mark as used.
    entropies[0] = np.inf
    
    # just for checking:
#    print("The entropies")
#    for index,_ in enumerate(entropies):
#        print(index,_,index_vector[index])
    ####
    
    estimated_sep_matrix = []
        

    for k in np.arange(1,n_sources+1):
        min_entropy_index = np.argmin(entropies)
        
        #just for checking
        #print("min_entropy_index",min_entropy_index)

        #Mark as used
        entropies[min_entropy_index] = np.inf

        if( k == 1):
            estimated_sep_matrix.append(index_vector[min_entropy_index])        
        
        estimated_sep_matrix = np.array(estimated_sep_matrix,ndmin=2)

        if (k > 1):
            index_vec_as_array = np.array(index_vector[min_entropy_index],ndmin=2)

            estimated_sep_matrix = np.vstack((estimated_sep_matrix,index_vec_as_array))
        

        
        length_k_index_vector = bss.createReversedIndexVectorSize(k,prime)
        #Mark any combination of the rows as used
        for scalars in length_k_index_vector:
            the_vector = np.array(scalars).dot(estimated_sep_matrix)
            the_vector = the_vector % prime
            
            entropies = entropies.reshape(occur_tensor.shape)
    
            the_vector = list(the_vector)
    
            the_vector.reverse()
            # one may wonder why reversing? It's because the values in entropies
            # are already in reversed_index organization. So if we want to 
            # to access the value related to the tuple (3,2,1)
            # we may access (1,2,3) in entropies, while it's in tensor shape
    
            the_vector = tuple(the_vector)
            
            entropies[the_vector] = np.inf
        
        # entropies comes back to array 1-D form
        entropies = entropies.ravel()
        # just for checking:
#        print("The entropies")
#        for index,_ in enumerate(entropies):
#            print(index,_,index_vector[index])
        ####
    estimated_sep_matrix = np.array(estimated_sep_matrix)

    
    return np.array(estimated_sep_matrix,dtype=np.int64)



def execute_AMERICA_online(mean_entropies_list,prime,n_sources):
    K_index_vector = bss.createReversedIndexVectorSize(n_sources,prime)
    entropies = mean_entropies_list
    entropies = np.array(entropies).ravel()
    
    # n must be different from zero. So mark as used.
    entropies[0] = np.inf
    min_entropy_index = np.argmin(entropies)

    #Mark as used
#    print("entropies",entropies)
#    print("index vector", index_vector)
    entropies[min_entropy_index] = np.inf

    estimated_sep_matrix = []
    estimated_sep_matrix.append(K_index_vector[min_entropy_index])
    
    the_shape = bss.generateOccurrenceTensorShape(prime,n_sources)

#    print(entropies)
    length_1_index_vector = bss.createReversedIndexVectorSize(1,prime)
    #Mark any combination of the found tuple as used
    for scalar in length_1_index_vector:
        the_vector = np.array(scalar).dot(estimated_sep_matrix)
        the_vector = the_vector % prime
        
        entropies = entropies.reshape(the_shape)

        the_vector = list(the_vector)
        the_vector.reverse()
        the_vector = tuple(the_vector)
        
        entropies[the_vector] = np.inf

    entropies = entropies.ravel()
#    print(entropies)
#    print("the estimated sep",estimated_sep_matrix)

    for k in np.arange(2,n_sources+1):

        min_entropy_index = np.argmin(entropies)

        #Mark as used
        entropies[min_entropy_index] = np.inf
        
        
        test_matrix = np.array(estimated_sep_matrix,ndmin=2)
        index_vec_as_array = np.array(K_index_vector[min_entropy_index],ndmin=2)

        test_matrix = np.vstack((test_matrix,index_vec_as_array))
        
#        print("the test_matrix",test_matrix)
        
        length_k_index_vector = bss.createReversedIndexVectorSize(k,prime)
        #Mark any combination of the rows as used
        for scalars in length_k_index_vector:
            the_vector = np.array(scalars).dot(test_matrix)
            the_vector = the_vector % prime
            
            entropies = entropies.reshape(the_shape)
    
            the_vector = list(the_vector)
    
            the_vector.reverse()
    
            the_vector = tuple(the_vector)
            
            entropies[the_vector] = np.inf
        
        entropies = entropies.ravel()


        estimated_sep_matrix = test_matrix

    estimated_sep_matrix = np.array(estimated_sep_matrix)

    return estimated_sep_matrix





def trial_AMERICA_EqualPMFs_OnlineEntropy(n_samples,n_sources,src_0_pmf,mix_matrix_type='random',r=100):
    import time
    start_time = time.time()

    

    n_samples = n_samples
    n_sources = n_sources
    finite_field = len(src_0_pmf)
    prime = finite_field

    pmf_list = []
    for _ in range(n_sources):
        pmf_list.append(src_0_pmf)
        
    if(mix_matrix_type == 'random'):
        the_mixing_matrix = bss.generateMixingMatrix(n_sources,prime = finite_field)
    elif(mix_matrix_type == 'identity'):
        the_mixing_matrix = np.eye(n_sources,dtype=int)  
        
    
    assert(bss.round_det_modP(the_mixing_matrix,finite_field) != 0)
    
    the_mixing_matrix  = getGFarray(the_mixing_matrix,finite_field)
    
    the_samples  = []
    
    for _ in range(n_samples):
        generateObservations_AppendToSamples(the_mixing_matrix,the_samples,pmf_list)
        
    empirical_entropies_list = [ [] for _ in range(finite_field ** n_sources)]
    mean_entropies_list = [ [] for _ in range(finite_field ** n_sources)]
    
    
#    print(the_samples)
    K_index_vector = bss.createReversedIndexVectorSize(n_sources,finite_field)
    for index,the_tuple in enumerate(K_index_vector):
        samples_transformed = np.array(the_tuple).dot(np.array(the_samples).T) % prime
#        print("samples transformed")
#        print(list(samples_transformed))
        empirical_entropies_list[index] = online_Init(samples_transformed,r,log_base=finite_field)
        mean_entropies_list[index] = np.asarray(empirical_entropies_list[index],dtype=np.float64).mean()
        
    
#    print("the empirical entropies list")        
#    print(empirical_entropies_list)
#    print("the entropies")
#    print(mean_entropies_list)
#    
#    print(np.array(empirical_entropies_list).shape)


## THIS PART is just the second part of AMERICA ALGORITHM
## TO-DO, transform it into a function
    


#    print(entropies)
#    print(estimated_sep_matrix)
    sep_matrix = execute_AMERICA_online(mean_entropies_list,prime,n_sources)
    print("on-line AMERICA")
    print("Also: W*A")
    W_x_A = getINTarray(sep_matrix.dot(the_mixing_matrix))
    print(W_x_A)
    
    print("the number of unique hits:")
    n_unique_sep_sources = bss.numberOfUniqueHitsInSeparation(W_x_A)
    print(n_unique_sep_sources)
    
    if(n_unique_sep_sources == n_sources):
        print("Fully separated!")
    else:
        print("Not Totally separated...")
    
    
    print("The mix matrix")
    print(getINTarray(the_mixing_matrix))
    
    
    print("The sep_matrix")
    print(sep_matrix)
    print("Must be full rank: det != 0")
    
    
    
    
    assert(bss.round_det_modP(sep_matrix,finite_field) != 0)
    
    print("det of sep_matrix\n",bss.round_det_modP(sep_matrix,finite_field))

    print("run n_samples = ",n_samples," in: ","--- %s seconds ---" % (time.time() - start_time))

#    return (time.time() - start_time)
    return n_unique_sep_sources == n_sources
    
## Test Part:



