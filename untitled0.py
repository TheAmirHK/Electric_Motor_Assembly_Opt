# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 11:23:03 2023

@author: amirh
"""

import numpy as np
import random
import statistics 
import matplotlib.pyplot as plt
from tqdm import tqdm 

# In[] Random assembly
def Random_assembly(matrix_form, iterations,admissble_quality):
    
    random.seed(100)    
    
    IDs = np.arange(0, len(matrix_form), 1, dtype=int).tolist()    
    numbers = []
    perf = []
    
    for _ in tqdm(range(iterations)):
        rand_house = random.sample(IDs, len(IDs)) # generate random order of labels 
        rand_body = random.sample(IDs, len(IDs)) # generate random order of labels 
        rand_shaft = random.sample(IDs, len(IDs)) # generate random order of labels 

        t_list = []

        for i in range(len(rand_house)):
            t_list.append(matrix_form[rand_house[i]][rand_body[i]][rand_shaft[i]]) # This line returns assembly response 
               
        n = 0
        for element in t_list:
            
            if element > admissble_quality:
                n = n+1 # Number of conformed assemblies at each iteration
        numbers.append(n)

        new_perf = np.mean(t_list)
        
        perf.append(new_perf)
    evaluation = np.c_[perf,numbers] ### [The average value of quality at each iteration,  Numbner of conformed assemblies at each iteration]         
    print("\nAvergave number of assemblies = %.0f" %np.mean(numbers))                     
    print(["Min:", np.min(perf), ", Mean:", np.mean(perf), "Max:", np.max(perf)])
    
    ### plot histogram graph of the random assemblies
    var = statistics.stdev(perf)
    mean = np.mean(perf)
    bins = 200    
    n, bins, patches = plt.hist(perf, bins=bins, color = [(51/255,161/255,201/255)], density = 1)
    y = ((1 / (np.sqrt(2 * np.pi) * var)) *np.exp(-0.5 * (1 / var * (bins - mean))**2))    
    plt.plot(bins, y, '--',color = "black" )
    plt.axvline(mean, color='k', linestyle='dashed', linewidth=1)    
    plt.xlabel("Quality conformity rate")
    plt.ylabel("Frequency")
    
    return np.mean(evaluation[:,0]), evaluation