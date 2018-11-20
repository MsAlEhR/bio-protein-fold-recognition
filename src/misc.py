# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 18:32:47 2018

@author: Mir, A.
"""

import itertools
import matplotlib.pyplot as plt
import numpy as np



def plt_confusion_matrix(cm, classes, title='Confusion matrix',
                         cmap=plt.cm.Blues):
    
    """
    This function plots the confusion matrix
    """
    
    # Normalizing confusion matrix
    cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=90)
    plt.yticks(tick_marks, classes, rotation=90)
    
    fmt = '.2f'
    thresh = cm.max() / 2.0
    
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        
        plt.text(j, i, format(cm[i, j], fmt), horizontalalignment="center", \
                 color="white" if cm[i, j] > thresh else "black")
        
    plt.ylabel("True label")
    plt.xlabel("Predicted label")
    plt.tight_layout()
    
