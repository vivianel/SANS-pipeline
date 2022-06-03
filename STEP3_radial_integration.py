# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 17:20:35 2022

@author: lutzbueno_v
"""

path_dir = 'C:/Users/lutzbueno_v/Documents/GitHub/pipeline_SANS1/Data_Gaia/'
file_beam_center = 'EB' #usually we use the empty beam



#load the data
import pyFAI
import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
from calculate_beam_center import calculate_center


plt.close('all')
# find all files in the folder
files = []
path_dir_an = os.path.join(path_dir, 'analysis/')
list_dir = list(os.listdir(path_dir_an))

#print the list
#for f in files:
#    print(f)


for kk in list_dir:
    if kk[0:3] == 'det':
        path_det = os.path.join(path_dir_an, str(kk))
        file_name = os.path.join(path_det, 'det_files_'+ kk[4:] + 'm.json')
        class_file = open(file_name)
        class_file = eval(class_file.read())

        # %% create poni and masks
        path_rad_int = os.path.join(path_det, 'integration/')
        if not os.path.exists(path_rad_int):
            os.mkdir(path_rad_int)
        path_rad_int_fig = os.path.join(path_det, 'figures/')
        if not os.path.exists(path_rad_int_fig):
            os.mkdir(path_rad_int_fig)
        
        if file_beam_center in class_file['name']:
            #open the file for the beam ceneter and radial integrator
            idx = class_file['name'].index(file_beam_center)
        else:
            print('###########################################################')
            print('There is no Empty beam measurement for this configuration!')
            print('###########################################################')
            idx = 1
        path_hdf_raw = os.path.join(path_det, 'hdf_raw/')
        name_hdf = path_hdf_raw + '/'+ class_file['name_hdf'][idx]
    
        file_hdf = h5py.File(name_hdf, 'r')
        
        # get values from the hdf file
        #ai = pyFAI.load('C:/Users/lutzbueno_v/Documents/GitHub/pipeline_SANS1/CH1_SAXS.poni')#azimuthalIntegrator().AzimuthalIntegrator()
        dist = (float(np.asarray(file_hdf['/entry1/SANS/detector/x_position']))/1000)
        det_disp = (float(np.asarray(file_hdf['/entry1/SANS/detector/y_position']))/1000)
        pixel1 = 7.5e-3 # detector size is 7.5 x 7.5 mm²
        pixel2 = 7.5e-3 # detector size is 7.5 x 7.5 mm²
        wavelength = (float(np.asarray(file_hdf['/entry1/data1/lambda'])))/1e8 # in m from nm
        
        # calculate the beam center
        img = np.array(file_hdf['entry1/data1/counts'])
        img1 = np.where(img==0, 1e-4, img)
        det_img = np.log(img1)
        if file_beam_center in class_file['name']:
            bc_x, bc_y = calculate_center(det_img)
            file_name = path_rad_int_fig + 'beam_center_' + class_file['scan'][idx] + '_' + class_file['name'][idx] + '_' +kk[4:] + 'm.jpeg' 
            plt.savefig(file_name)
            plt.close('all')
        else:
            bc_x = round(float(np.asarray(file_hdf['/entry1/SANS/detector/beam_center_x'])), 2)
            bc_y = round(float(np.asarray(file_hdf['/entry1/SANS/detector/beam_center_y'])), 2)
        
        poni2 = bc_x*pixel1;
        poni1 = bc_y*pixel2;
        
        
        #create the radial integrator
        ai = pyFAI.AzimuthalIntegrator(dist=dist, poni1=poni1, poni2=poni2,
         rot1=0, rot2=0, rot3=0,
         pixel1=pixel1, pixel2=pixel2,
         splineFile=None, detector=None, wavelength=wavelength)
        print("\nIntegrator: \n", ai)
        
        # create a mask
        mask = np.zeros([128, 128])
        mask[:,0] = 1
        mask[:,-1] = 1
        mask[0,:] = 1
        mask[-1,:] = 1
        corner = 2
        mask[0:corner,0:corner] = 1
        mask[-corner:-1,0:corner] = 1
        mask[-corner:-1,-corner:-1] = 1
        mask[0:corner,-corner:-1] = 1
        beam_stopper = 6
        mask[round(bc_y)-beam_stopper: round(bc_y)+beam_stopper, round(bc_x)-beam_stopper: round(bc_x)+beam_stopper] = 1
        mask_inv =  (mask == 0).astype(int)
        #plt.imshow(mask, cmap='gray', origin='lower')
        
        # execute the radial inrtegration for all
        
        
        for ii in range(0, len(class_file['name'])):
            name_hdf = path_hdf_raw + '/'+ class_file['name_hdf'][ii]
            file_hdf = h5py.File(name_hdf, 'r')
            # radial integration
            img = np.array(file_hdf['entry1/data1/counts'])
            img1 = np.where(img<=0, 1e-4, img)
            #img1 = ndimage.gaussian_filter(img1, sigma=(0.25), mode = 'nearest')
            
            #plt.imshow(np.multiply(img2,mask_inv), clim=(0, 8), origin='lower', cmap = 'jet')
            #plt.colorbar()
            file_name = path_rad_int + 'integ_' + class_file['scan'][ii] + '_' + class_file['name'][ii] + '_' +kk[4:] + 'm.dat'           
            q, I, sigma = ai.integrate1d(img1,  120, correctSolidAngle=True, mask=mask,
                                   method = 'nosplit_csr', unit = 'q_nm^-1', safe=True, error_model="azimuthal",
                                   filename=file_name)
            
            # plots
            # Convert pixel coordinates starting at the beam center to coordinates in the inverse space (unit: nm ^ -1)
            def x2q(x, wv, dist, pixelsize): 
                return 4*np.pi/wavelength*np.sin(np.arctan(pixelsize*x/dist)/2)
            
            qx = x2q(np.arange(img.shape[1])-bc_x, wavelength, dist, pixel1)
            qy = x2q(np.arange(img.shape[0])-bc_y, wavelength, dist, pixel2)
            extent = [qx.min(), qx.max(), qy.min(), qy.max()]
            
            
            Int = 2
            img2 = img1#np.log(img1)
            if img2[img2>0].mean() > 0:
                clim1 = img2[img2>0].mean()-Int*img2[img2>0].std()
                if clim1 > 0:
                    clim1 = (clim1, img2[img2>0].mean()+Int*img2[img2>0].std())
                else:
                    clim1 = (0, img2[img2>0].mean()+Int*img2[img2>0].std())
            else:
                clim1 = [0, 1]
            fig, axs = plt.subplots(1, 2, num = 1,  figsize=(10, 6), gridspec_kw={'width_ratios':[2, 2]})
            im1 = axs[0].imshow(img2, origin='lower', clim =clim1, cmap = 'jet', extent = np.divide(extent,1e9)) # to have in nm
            fig.colorbar(im1, ax = axs[0], orientation = 'horizontal', shrink = 0.75).set_label('log(Intensity)')
            axs[0].grid(color = 'white', linestyle = '--', linewidth = 0.25)
            axs[0].set(ylabel = r'q$_{y}$ [$nm$$^{-1}$]', xlabel = r'q$_{x}$ [$nm$$^{-1}$]')
                   
            axs[1].plot(q[0:-3], I[0:-3], '*k')
            axs[1].set(xlabel = r'Scattering vector q [$nm^{-1}$]', ylabel = 'Intensity I [a.u.]', xscale = 'log', 
                       yscale = 'log', title = 'Sample: '+  class_file['name'][ii], xlim = [0.001, 0.1])
            axs[1].grid(color = 'gray', linestyle = '--', linewidth = 0.5)
            axs[1].errorbar(q[0:-3], I[0:-3], yerr = sigma[0:-3], color = 'black', lw = 1)
            
            file_name = path_rad_int_fig + 'integ_' + class_file['scan'][ii] + '_' + class_file['name'][ii] + '_' +kk[4:] + 'm.jpeg' 
            plt.savefig(file_name)
            plt.clf()




