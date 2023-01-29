# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 11:41:02 2022

"""

from pyomo.environ import * 
from pyomo.opt import SolverFactory
import numpy as np
import matplotlib.pyplot as plt
import statistics 

# In[] Individual assembly

def Individual_assembly (matrix_form, admissble_quality):   
    population_size = len (matrix_form)
    rang = population_size
    model = ConcreteModel(name="(pairing)")
    A = [(i,j,k) for i in range(rang) for j in range(rang) for k in range(rang) ]
    vel = [i for i in range(rang)]
    model.x = Var(A, within = Binary)
    model.y = Var(A, within = Binary)
    model.number = Var(initialize=0)
    
    def objective_rule(model):
        QD = sum((matrix_form[i,j,k]*model.x[i,j,k] - 1000 * model.y[i,j,k]) for i in range (rang) for j in range (rang) for k in range(rang))
        
        return QD
    model.objective = Objective(rule=objective_rule, sense=maximize)                                  
#---------Define constriants---------
    def Each_column(model, j):
        return (sum(model.x[i,j,k] + model.y[i,j,k]  for i in range(rang ) for k in range(rang ))==1)
    model.constraint = Constraint(vel, rule=Each_column)
    
    def Each_raw(model, i):
        return (sum(model.x[i,j,k] + model.y[i,j,k]  for j in range(rang) for k in range(rang )) ==1)
    model.constraint2 = Constraint(vel, rule=Each_raw)
    
    def Each_diameter(model, k):
        return (sum(model.x[i,j,k] + model.y[i,j,k]  for j in range(rang) for i in range(rang )) ==1)
    model.constraint3 = Constraint(vel, rule=Each_diameter)
                     
    opt = SolverFactory('glpk')
    results = opt.solve(model, tee=True)   
    m = 0
    n = 0

    house_ = []
    body_ = []
    shaft_ = []
    pairs = np.zeros ((population_size,4))
    
    for i in range (rang):
        for j in range(rang):
            for k in range(rang):
            
                if model.x[i,j,k]==1:
                    m= m+1
                    pairs[i] = [i,j,k,matrix_form[i,j,k]]
                    print (i,j,k,matrix_form[i,j,k])                
        if pairs[i,3] < admissble_quality:
            n = n + 1
            house_.append(pairs[i,0])
            body_.append(pairs[i,1])
            shaft_.append(pairs[i,2])
                        
    ### plot histogram graph of the random assemblies
    perf = pairs[:,3]
    var = statistics.stdev(perf)
    mean = np.mean(perf)
    bins = 100    
    n, bins, patches = plt.hist(perf, bins=bins, color = [(51/255,161/255,201/255)], density = 1)
    y = ((1 / (np.sqrt(2 * np.pi) * var)) *np.exp(-0.5 * (1 / var * (bins - mean))**2))    
    plt.plot(bins, y, '--',color = "black" )
    plt.axvline(mean, color='k', linestyle='dashed', linewidth=1)    
    plt.xlabel("Quality conformity rate")
    plt.ylabel("Frequency")
                        
    print ("Residual house",  np.asarray(house_))
    print("Residual body",np.asarray(body_))
    print("Residual shaft",np.asarray(shaft_))       
    print ("Number of residuals %0.f"%n)
    print("number of pairs %0.f"%(m-n))

    print ( "Mean KTE vale = %0.2f"%(np.mean(pairs[:,3])))              
    return pairs

Individual_assembly (matrix_form,0.95)