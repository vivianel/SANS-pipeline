# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 14:50:44 2022

@author: lutzbueno_v
"""
path_dir = 'C:/Users/lutzbueno_v/Documents/GitHub/pipeline_SANS1/Data_Gaia/'

# sample to consider for the creation of the transmission mask
mask_measurement = 'EB'

import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from tabulate import tabulate
from contextlib import redirect_stdout
import json

plt.close('all')
path_dir_an = os.path.join(path_dir, 'analysis/')
path_transmission = os.path.join(path_dir_an, 'transmission/')

class_trans = open(os.path.join(path_transmission, 'trans_files.json'))
class_trans = eval(class_trans.read())

#%% find the reference value
for ii in range(0, len(class_trans['name'])):
    if class_trans['name'][ii] == mask_measurement:
        path_hdf_raw = os.path.join(path_transmission, 'hdf_raw/')
        name_hdf = path_hdf_raw + '/'+ class_trans['name_hdf'][ii]
        file_hdf = h5py.File(name_hdf, 'r')
        img = np.array(file_hdf['entry1/data1/counts'])
        # 1e-4 is chosen as a small value
        img1 = np.log(np.where(img==0, 1e-4, img))
        # calculation of the cuttoff value for the mask
        cutoff = img1[img1>0].mean()
        plt.figure()
        plt.imshow(img1, clim=[0, round(cutoff)], cmap='jet', origin='lower')
        plt.colorbar(orientation = 'vertical', shrink = 0.5).set_label('log(Intensity)')
        im_title = 'Empty_beam'
        plt.title(im_title)
        im_title = str(path_transmission + im_title + '.jpg')
        plt.savefig(im_title)
        img2 = np.where(img1<cutoff, 0, img1)
        mask = np.where(img2>=cutoff, 1, img2)
        plt.figure()
        plt.imshow(mask, cmap='gray', origin='lower')
        plt.colorbar(orientation = 'vertical', shrink = 0.5, ticks=[0, 1]).set_label('Binary')
        im_title = 'Transmission_Mask'
        EB_ref = int(np.sum(np.multiply(img,mask)))
        print('###########################################################')
        print('This is the EB transmission:' + str(EB_ref))
        print('###########################################################')
        plt.title(im_title +  ', Total counts = ' + str(EB_ref))
        im_title = str(path_transmission + im_title + '.jpg')
        plt.savefig(im_title)
        
        
#%% calculate all transmissions
list_trans = []
list_counts = []
path_hdf_raw = os.path.join(path_transmission, 'hdf_raw/')
for ii in range(0, len(class_trans['name'])):
        name_hdf = path_hdf_raw + '/'+ class_trans['name_hdf'][ii]
        file_hdf = h5py.File(name_hdf, 'r')
        img = np.array(file_hdf['entry1/data1/counts'])
        counts = int(np.sum(np.multiply(img,mask)))
        list_counts.append(counts)
        print('This is sample '+ class_trans['name'][ii] + ' transmission:' + str(counts))
        trans = np.divide(counts,EB_ref)
        trans = round(trans, 3)
        list_trans.append(trans)
        print(trans)
class_trans['transmission'] = list_trans
class_trans['counts'] = list_counts
# print the updated list of the transmission files
df_trans = pd.DataFrame(class_trans)
data = tabulate(df_trans, headers='keys', tablefmt='psql')
print(data)
save_file = os.path.join(path_transmission, 'trans_files.txt')
with open(save_file, 'w') as f:
    with redirect_stdout(f):
        print(data)
save_file = os.path.join(path_transmission, 'trans_files.json')
a_file = open(save_file, "w")
json.dump(class_trans, a_file)
a_file.close()