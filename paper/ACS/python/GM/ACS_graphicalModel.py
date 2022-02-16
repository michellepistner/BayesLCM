# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 13:58:35 2019

@author: m_pis
"""
import os
os.chdir("C:\\Users\\m_pis\\Box Sync\\Python Examples (ACS)\\Graphical Models Example\\private-pgm-master\\src")


import numpy as np
from scipy import sparse
import pandas as pd
import matplotlib
import networkx
from mbi import Dataset, Domain, FactoredInference
import csv


##Now for a different example
# load ACS data set
os.chdir("C:\\Users\\m_pis\\Box Sync\\Python Examples (ACS)\\Graphical Models Example")

E = [0.25,0.50,1.0]

ep_names = ["0_25","0_50","1_00"]
k = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100"]

# spend half of privacy budget to measure all 1 way marginals
np.random.seed(0)

for q in range(0, 3):
    for kk in range(0,100):
        file_name = 'ACSsubset_' + k[kk] + '.csv'
        data = Dataset.load(file_name, 'acs-domain.json')
        domain = data.domain
        total = data.df.shape[0]
        epsilon = E[q]/2
        sigma = 2.0 / (epsilon/len(data.domain))
        
        measurements = []
        for col in data.domain:
            x = data.project(col).datavector()
            y = x + np.random.laplace(loc=0, scale=sigma, size=x.size)
            I = sparse.eye(x.size)
            measurements.append( (I, y, sigma, (col,)) )

# spend half of privacy budget to measure some more 2 and 3 way marginals

        cliques = [('CIT', 'AGEP'), 
                    ('CIT', 'RACWHT'), 
                    ('CIT', 'SEX'),
                    ('CIT', 'PINCP'),
                    ('AGEP', 'RACWHT'),
                    ('AGEP', 'SEX'),
                    ('AGEP', 'PINCP'),
                    ('RACWHT', 'SEX'),
                    ('RACWHT', 'PINCP'),
                    ('SEX', 'PINCP')]
        
        sigma = 2.0 / (epsilon/len(cliques))
        
        for cl in cliques:
            x = data.project(cl).datavector()
            y = x + np.random.laplace(loc=0, scale=sigma, size=x.size)
            I = sparse.eye(x.size)
            measurements.append( (I, y, sigma, cl) )

        # now perform inference to estimate the data distribution
        
        engine = FactoredInference(domain, log=True, iters=10000)
        model = engine.estimate(measurements, total=total, engine='RDA')

        # now answer new queries
        
        privTable = model.project(('CIT', 'AGEP', 'RACWHT','SEX','PINCP')).datavector()
        Table = data.project(('CIT', 'AGEP', 'RACWHT','SEX','PINCP')).datavector()
        
        file_name = "GM_"+ ep_names[q] + "_" + k[kk] + ".csv"
        
        resultFile = open(file_name,'w')
        wr = csv.writer(resultFile, delimiter=",")
        wr.writerows([privTable])


