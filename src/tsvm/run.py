#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 10:02:36 2018
@author: Mir, A.

In this module, TSVM classifier can be trained with GridSearchCV and final model
will be saved on disk. Also, results of GridSearch will be saved in log file.
"""

from twinsvm import OVO_TSVM
from dataproc import read_data
from sklearn.model_selection import GridSearchCV
from joblib import dump
import sys
old_stdout = sys.stdout

# Open a file to save logs
log_file = open('../report/tsvm_dd_acc.txt', 'w')
sys.stdout = log_file

# Reading dataset file
train_data, labels, data_name = read_data('../dataset/dd_ACC_num.csv')

# Determine range of parameters
kernel = 'RBF'

c_range = {'C_1': [float(2**i) for i in range(-8, -6)],
             'C_2': [float(2**i) for i in range(-8, -6)]}

gamma_range = {'gamma': [float(2**i) for i in range(-10, -8)]} if kernel == 'RBF' else {}

param_range = {**c_range, **gamma_range}

# Initialize TSVM model
tsvm_model = OVO_TSVM(kernel=kernel)

# Arguments for grid search
cv_fold = 5
n_workers = 1 # Number of CPU threads

result = GridSearchCV(tsvm_model, param_range, cv=cv_fold, n_jobs=n_workers, refit=True,
                      verbose=10)

result.fit(train_data, labels)

# Save best estimator on disk
dump(result.best_estimator_, './tsvm_dd_acc_best_model.joblib')

sys.stdout = old_stdout
log_file.close()

