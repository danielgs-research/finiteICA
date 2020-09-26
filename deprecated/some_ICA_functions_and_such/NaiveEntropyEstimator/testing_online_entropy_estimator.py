# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 16:54:30 2018

@author: MateusMarcuzzo
Inspired in original Article from:
    Cédric Lauradoux et al
    "Online Entropy Estimation for Non-Binary Sources and Applications on IPhone"
    
Here we have the transcription in python, from the MatLab code of 
    Prof. Daniel Guerreiro trying to test the entropy estimator

later we have a mean of values collected, which converges to the true entropy
    in multiple trials.
"""



import numpy as np
from numpy.random import choice
from entropy_estimators import entropyRandomVariable
from entropy_estimators import entropyRandomVariableFromArray as entropy_on_array
from entropy_estimators import onlineEntropyLauradouxEstimator_Initialization as online_entropy_array
from entropy_estimators import onlineEntropyLauradouxEstimator_LatestSample as online_entropy_update
import matplotlib.pyplot as plt



def calculateNthHarmonicalNumber(n):
    the_sum = 0
    for i in np.arange(1,n+1):
        the_sum += 1/i
        
    return the_sum

#for i in range(100):
#    print(calculateNthHarmonicalNumber(i))


n_samples = 1000
r  = 100

assert(r < n_samples)



log_2 = np.log(2)

# r==1 is a special case. It's descripted in the paper
if(r == 1):
    adjust = 2
else:
    adjust = 1/log_2

entropy_mean_values_among_runs = []

possible_events  = ['head','tail']
probs = [0.55,0.45]
samples = choice(possible_events, size = n_samples, p = probs)
samples = list(samples)

traditional_entropy = entropy_on_array(samples,possible_events,log_base = len(possible_events))

#uncomment for uniform case
#true_entropy = entropyRandomVariable(len(possible_events)*list([1]),log_base=len(possible_events))

#a custom case using probs
true_entropy = entropyRandomVariable(probs,log_base=len(possible_events))


## Esquema do Daniel.

# adjustable
lambda_value = 0.001

#the 'r' in the paper
## "window size" ?


entropy_previous  = 0
entropy_empirical = np.full(n_samples-r,0,dtype=np.float64)
for i in np.arange(r,n_samples):
    # don't have a name at all
    w = 0
    
    # don't have a name too
    j = 1
    
    while(j <= r and samples[i] != samples[i-j]):
        w = w + 1/j
        j = j + 1
        
    entropy_empirical[i-r] = (1-lambda_value)*entropy_previous + lambda_value * w * adjust
    
    entropy_previous = entropy_empirical[i-r]
        
#mudança de base
entropy_empirical = entropy_empirical/np.log2(len(possible_events))
#    
#    print("difference in entropies")
#    print(entropy_empirical[-1] - traditional_entropy)
##    
#
#    plt.plot(np.arange(r+1,n_samples+1),entropy_empirical)
#    plt.plot(np.arange(r+1,n_samples+1),traditional_entropy*np.full(n_samples-r,1,dtype=np.float64))
#    
#    plt.show()
#    plt.clf()

######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################

## Ensaio dois, mudando a média móvel:


#entropy_empirical = np.full(n_samples-r,0,dtype=np.float64)
# Actually, it does'nt need to be an array.
#entropy_empirical = [0] * (n_samples-r)

entropy_empirical = [ 0 for i in range(n_samples-r)]

for i in np.arange(r,n_samples):
    # don't have a name at all
    w = 0
    
    # don't have a name too
    j = 1  
    while(j <= (r-1) and samples[i] != samples[i-j]):
        w = w + 1/j
        j = j + 1
        
    # H_P(X) = H_2(X)/log(P)_base(2)
    entropy_empirical[i-r] =  w *  adjust /np.log2(len(possible_events))
    
#    print(entropy_empirical)
#    plt.plot(np.arange(len(entropy_empirical)),entropy_empirical)
#    plt.show()
#    plt.clf()
# Quero fazer uma média que seja incremental também, senão vai ter que sempre
## somar todo mundo, dividir pelo total
## Isso também não dura pra sempre...
#print(entropy_empirical)
entropy_empirical_mean = np.asarray(entropy_empirical,dtype=np.float64).mean()
entropy_mean_values_among_runs.append(entropy_empirical_mean)

print("the true entropy")
print(true_entropy)

print("the Naive entropy")
print(traditional_entropy)

print("The online entropy")
print(entropy_empirical_mean)


print("(true - Naive)")
print(true_entropy - traditional_entropy)

print("(Naive - online)")
print( traditional_entropy - entropy_empirical_mean )

#    
#    print("(true - online)")
#    print(true_entropy-entropy_empirical_mean)

#print("the true entropy")
#print(true_entropy)
#
#print("entropy mean values among runs")
#print(sum(entropy_mean_values_among_runs)/len(entropy_mean_values_among_runs))

plt.plot(np.arange(len(entropy_mean_values_among_runs)),entropy_mean_values_among_runs,label='entropy mean values among runs')
plt.show()
plt.clf()

plt.plot(np.arange(r+1,n_samples+1),traditional_entropy*np.full(n_samples-r,1,dtype=np.float64),label="traditional")

plt.plot(np.arange(r+1,n_samples+1),entropy_empirical_mean*np.full(n_samples-r,1,dtype=np.float64),label="online")
plt.legend()


plt.show()
plt.clf()

############################################
# the following method just does what the above code does.
print("the entropy empirical mean")
print(entropy_empirical_mean)
print("the value using a method")
print(online_entropy_array(samples,r,log_base=len(possible_events)) . mean())


## The code Below illustrates the difference in speed doing mean with all the
## dataset and with just doing a recursive mean calculation.
#old_mean = entropy_empirical_mean
#counter = 0
#while(True):
#    samples.append('tail')
#    new_entropy = online_entropy_update(samples,r,log_base=len(possible_events))
#    entropy_empirical.append(new_entropy)
#    
#    
#    ## doing mean on all entropies
#    #print(online_entropy_array(samples,r,log_base=len(possible_events)).mean())
#    
#    
#    # recursive mean
#    ## Using the Method Suggest by Prof. Daniel Guerreiro
#    
#    ## M_(n+1) = n/(n+1) * M_n + H_(n+1)/(n+1)
#    
#    ## M_n designa a n-ésima média, H_n a n-ésima estimativa de entropia.
#    ## Note que esta entropia é calculada baseada no r, que depende das
#    ## das amostras (tuplas) que chegam no processo de BSS
#    new_mean = (len(entropy_empirical)-1)/len(entropy_empirical) * old_mean + new_entropy/len(entropy_empirical)
#    
#    old_mean = new_mean
#    
#    print(new_mean)
#    counter+=1
#    print("counter ", counter)
    
