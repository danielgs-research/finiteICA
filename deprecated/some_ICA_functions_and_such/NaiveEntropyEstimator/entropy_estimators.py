# -*- coding: utf-8 -*-
"""
Transfered the content here to BSS_functions in 26/02/2019
It just remained some tests.
"""
import numpy as np



############################################################################
    


## This should be put separately
    ##Needs both bss and entropy_estimators
    ## genPMF = generateRandomProbabilityMassFunction
    ## entropy = entropyRandomVariable
    
#def testing_o_gini_vs_entropy():
#    n_sources_ = np.arange(3,30)
#for n_sources in n_sources_:
#    print('n_sources = ',n_sources)
#    n_samples = 10**5
#    
#    prob_ = []
#    
#    for _ in range(n_samples):
#        prob_.append(genPMF(n_sources))
#    
#    before = time.time()
#    for prob in prob_:
#        o_gini(prob)
#    
#    print(time.time()-before)
#
#    before = time.time()
#    for prob in prob_:            
#        entropy(prob,len(prob))
#    
#    print(time.time() - before,"\n")
   
    
    
############################################################################
    
    
    
    
## This is also very important
    # We are measuring how much the o_gini is less-than-equal to entropy.
    # so we have a probabilistically monotonic dependance of gini and entropy!
    
#    def combine(integer):
#    prob =  genPMF(integer)
#    return o_gini(prob) <= entropy(prob,len(prob))
    
    
## ALso
#     counter_true = 0
#    for _ in range (10000):
#        if(combine(2)):
#            counter_true+=1
#    print(counter_true/10000)
        
    
    
    
############################################################################
    
    
    
## Muito importante também! Teste completo. Gráfico em 25/05/2019
    ## EDIT IN 27/02/2019
    # n_sources tem o tamanho dos corpos! Não a quantidade de fontes!!!
#start_time = time.time()
#pmf_sizes  = np.arange(2,30)
#
#n_sources_ = pmf_sizes
#
#total_time_o_gini = []
#total_time_entropy = []
#combine_success_rate = []
#for n_sources in n_sources_:
#    n_samples = 10**5
#    
#    prob_ = []
#    
#    counter_true = 0
#    for _ in range(n_samples):
#        prob_.append(genPMF(n_sources))
#        if(combine(n_sources)):
#            counter_true+=1
#    
#    combine_success_rate.append(counter_true/n_samples)
#    
#
#    before = time.time()
#    for prob in prob_:
#        o_gini(prob)
#        
#    total_time_o_gini.append(time.time()-before)
#    
#    before = time.time()
#    for prob in prob_:
#        entropy(prob,len(prob))
#        
#    total_time_entropy.append(time.time() - before)    
#
#plt.plot(n_sources_,total_time_o_gini,label='o_gini_total')
#plt.plot(n_sources_,total_time_entropy,label='entropy_total')
#
#plt.legend()
#plt.show()
#
#plt.plot(n_sources_,combine_success_rate,label='combine_o_gini_rate')
#plt.legend()
#plt.show()
#
#plt.clf()
#
#print('elapsed total time ', time.time()-start_time)