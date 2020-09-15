# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 16:51:21 2019

@author: MateusMarcuzzo
"""

from AMERICA import *

algorithm = executeAMERICAAlgorithm
sample_generator = generateObservations_UpdateTensor

america_trial = trial(algorithm,sample_generator)

#    np.random.seed(1)
n_samples = 20000
prime = 2
n_sources = 11
pmf_list = []
for _ in range(n_sources):
    pmf_list.append(bss.generateRandomProbabilityMassFunction(prime))
    
#print("Distribuições")
#for _ in pmf_list:
#    print(_)
    
mix_matrix = bss.generateMixingMatrix(n_sources,prime)


#    o_gini_trials = 0
#    entropy_trials = 0
#    for k in range(10):
#        trials = 0
#        while(america_trial.execute_trial(n_samples,pmf_list,mix_matrix,"o_gini",verbosity=True)):
#           trials+=1
#            #print(trials)
#        
#        print("trials until fail: ", trials)
#        o_gini_trials+=trials
#        
#        
#        trials=0
#        while(america_trial.execute_trial(n_samples,pmf_list,mix_matrix,"entropy",verbosity=False)):
#            trials+=1
#           # print(trials)
#        
#        print("trials until fail: ", trials)
#        entropy_trials+=trials
#        print("#############")
#              
#    print("o_gini trials until fail ",o_gini_trials)
#    print("entropy trials until fail ",entropy_trials)

america_trial.execute_trial(n_samples,pmf_list,mix_matrix,"entropy",verbosity=True)