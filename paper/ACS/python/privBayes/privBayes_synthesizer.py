# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:42:29 2019

@author: m_pis
"""

import os, sys
os.chdir("C:\\Users\\m_pis\\Box Sync\\Python Examples (ACS)\\PrivBayes\\DataSynthesizer\\DataSynthesizer")

from DataDescriber import DataDescriber
from DataGenerator import DataGenerator
from ModelInspector import ModelInspector
from lib.utils import read_json_file, display_bayesian_network

import pandas as pd
import numpy
import scipy
import dateutil
%matplotlib inline

os.chdir("C:\\Users\\m_pis\\Box Sync\\Python Examples (ACS)\\PrivBayes\\DataSynthesizer")
# input dataset
input_data = './data/ACSsubset.csv'
# location of two output files
mode = 'correlated_attribute_mode'
description_file = f'./out/{mode}/1000descrip.json'

# An attribute is categorical if its domain size is less than this threshold.
# Here modify the threshold to adapt to the domain size of "education" (which is 14 in input dataset).
threshold_value = 10

# specify categorical attributes
categorical_attributes = {'cit': True}##{'CIT': True, 'AGEP': True, 'SEX': True, 'RACWHT': True, 'PINCP': True}

# specify which attributes are candidate keys of input dataset.
candidate_keys = {'id': True}

# A parameter in Differential Privacy. It roughly means that removing a row in the input dataset will not 
# change the probability of getting the same output more than a multiplicative difference of exp(epsilon).
# Increase epsilon value to reduce the injected noises. Set epsilon=0 to turn off differential privacy.
E = [0.25,0.50,1.0]

ep_names = ["0_25","0_50","1_00"]
k = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100"]


# The maximum number of parents in Bayesian network, i.e., the maximum number of incoming edges.
degree_of_bayesian_network = 2

# Number of tuples generated in synthetic dataset.
num_tuples_to_generate = 10000

os.chdir("C:\\Users\\m_pis\\Box Sync\\Python Examples (ACS)\\PrivBayes\\DataSynthesizer")

for q in range(0, 3):
     for kk in range(0,100):
        
        synthetic_data = f'./out/{mode}/privBayes_E{q}_{kk}.csv'
 
        epsilon = E[q]

        describer = DataDescriber(category_threshold=threshold_value)
        describer.describe_dataset_in_correlated_attribute_mode(dataset_file=input_data, 
                                                                epsilon=epsilon, 
                                                                k=degree_of_bayesian_network,
                                                                attribute_to_is_categorical=categorical_attributes,
                                                                attribute_to_is_candidate_key=candidate_keys)
        describer.save_dataset_description_to_file(description_file)
        
        display_bayesian_network(describer.bayesian_network)
        
        generator = DataGenerator()
        generator.generate_dataset_in_correlated_attribute_mode(num_tuples_to_generate, description_file)
        generator.save_synthetic_data(synthetic_data)
        
        #file_name = "privBayes_"+ ep_names[q] + "_" + k[kk] + ".csv"
        
        #resultFile = open(file_name,'w')
        #wr = csv.writer(resultFile, delimiter=",")
        #wr.writerows([privTable])