# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 10:52:00 2022

@author: lutzbueno_v
"""

def calculate_center(det_img):
    import skimage.io
    import numpy as np
    import matplotlib.pyplot as plt
    
    
    # Read in the image
    # Note - intensities are floating point from [0,1]
    cutoff = np.max(det_img)-det_img[det_img>0].mean()
    im = np.where(det_img<cutoff, 0, det_img)
    im = np.where(im>=cutoff, 1, im)
    
    
    # Threshold the image first then clear the border
    im_clear = skimage.segmentation.clear_border(im > (200.0/255.0))
    
      
    # Show image in figure and hold to place dots in
    plt.figure()
    plt.imshow(np.dstack([im,im,im]))
    
  #  # For each image...
  #  for idx in range(5):
    
  #    # Extract sub image
    img = im_clear#[:,idx*split_point:(idx+1)*split_point]
    
    # Find coordinates of thresholded image
    y,x = np.nonzero(img)
  
    # Find average
    xmean = x.mean()
    ymean = y.mean()
  
    # Plot on figure
    plt.plot(xmean, ymean, 'r+', markersize=5)
  
    bc_x = xmean
    bc_y =  ymean
    plt.title('x_center = ' + str(round(bc_x, 2)) + ', y_center =' + str(round(bc_y, 2)) + ' pixels')
    
    # Show image and make sure axis is removed
    plt.axis('off')
    plt.show()
    
    return bc_x, bc_y