import numpy as np
import os
import shutil
import posixpath

import matplotlib.pyplot as plt

import wfdb


for i in range (11,21):

    a = 'a' + str(i) 

    f = open(a + '.txt', "w")

    record = wfdb.rdrecord('C:\\Users\\Denis\\Documents\\Uni\\Apnea sleep detector\\apnea-ecg-database-1.0.0\\apnea-ecg-database-1.0.0\\' + a)

    signal,sampling_freq = wfdb.plot_wfdb(record=record, title='')

    #print(signal[0][3])


    for i in range(0, signal[0].shape[0]-1):
     f.write(str(signal[0][i]) + '\n')

    print(sampling_freq) #un campione ogni 100ms 

    f.close()

    





