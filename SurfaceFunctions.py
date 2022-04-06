#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 23:02:49 2021

@author: mrr
"""

def mrr_percolation(field,p0,pf,nos,ifplot):
    """
    * field is the matrix where percolating network will be searched for \n
    * p0 and pf are the start and end points of p array \n
    * nos is the number of elements between p0 and pf \n
    * If ifplot is set to 1, figures are plotted at every p value, if set to 0
    figures are plotted only when a network is found \n
    * returned variables:
        1. Pcr
        2. [p, csizeAvg, csizeMax, num]
        3. array of cluster size/area
    [Ref.] Percolation with Python by Prof. Anders Malthe-Sorenssen, Univ, of Oslo
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.ndimage import measurements
    
    pvect =  np.linspace(p0,pf,nos)
    ln = len(pvect)
    pCr = np.nan
    csizeAvg = np.nan
    csizeMax = np.nan
    
    pstats = np.nan*np.ones(ln*4)
    pstats = pstats.reshape(ln,4)

    count = 0    
    for k in range(len(pvect)):
        p = pvect[k]
        m = field < p
        ######## NEXT-NEAREST NEIGHBOUR #######
        s = [[1,1,1],
             [1,1,1],
             [1,1,1]]
        lw, num = measurements.label(m,structure=s)

        
        ###########
        
        ####### NEAREST NEIGHBOUR #############
        
        #lw, num = measurements.label(m)
        
        ###########
        perc_xV = np.intersect1d (lw[0,:],lw[-1,:])
        perc_xH = np.intersect1d (lw[:,0],lw[:,-1])
        percV = perc_xV[(perc_xV>0)]    # vertical percolation
        percH = perc_xH[(perc_xH>0)] 
        
        # in lw, the labels are mere sequential. the following line, instead, labels the clusters
        # based on the number of sites in each cluster. So if 5 sites with labels 1 form a cluster, 
        # it will rename all these 5 clusters as 5
        labelList = np.arange(lw.max()+1)
        area = measurements.sum(m, lw, labelList)
        areaIm = area[lw]
        
        # average and maximum cluster size
        csizeAvg = area.mean() 
        csizeMax = area.max()
        
        pstats[count,:] = [p, csizeAvg, csizeMax, num]
        count = count + 1
        
        if ifplot == 1: # show images before critical p
            plt.figure(figsize=(10,10))
            plt.imshow(areaIm,cmap = 'RdYlBu_r') #cmap='RdYlBu_r') 
            plt.axis('off')
            tit = 'p = {}' 
            tit = tit.format(p)
            #plt.title(tit)
            plt.savefig('testperc.png',bbox_inches=None, pad_inches=0)
            plt.show()
    
        if len(percH)>0 and len(percV)>0:
            text = 'Percolation detected at p = {}'
            pstop = round(p,5)
            print(text.format(pstop))
            pCr = p 
            plt.figure(figsize=(10,10))
            plt.imshow(areaIm,cmap='RdYlBu_r') 
            plt.axis('off')
            plt.savefig('perc_both_dir.tiff',bbox_inches=None, pad_inches=0)
            
            break  
        
    return pCr, pstats, area
        

def mrr_surftenGlobal (fdir, startrec,endrec):
    import sys
    import numpy as np
    
    ppdir = '/Users/muhammadrizwanurrahman/Desktop/MainDrive/My_Storage/ICLResearch/PythonFIles/pyDataViewData/'
    sys.path.append(ppdir)
    
    import postproclib as ppl
    
    normal = 0 
    
    #Get Post Proc Object
    PPObj = ppl.All_PostProc(fdir) # has all the data files
    #Get plotting object
    PObj = PPObj.plotlist # Dictionary contains only the data for all fnames: PPObj.plotlist.keys() gives all the field names

    f_psurf = PObj["psurface"]
    a_psurf = f_psurf.read(startrec,endrec)
    a_psurfxx = a_psurf[:,:,:,:,0]
    a_psurfyy = a_psurf[:,:,:,:,4]
    a_psurfzz = a_psurf[:,:,:,:,8]
    locx, y = f_psurf.profile(axis=normal,startrec=startrec,endrec=endrec) # locx is location

    f_vflux = PObj["vflux"]
    a_vflux = f_vflux.read(startrec,endrec)
    a_vfluxxx = a_vflux[:,:,:,:,0]
    a_vfluxyy = a_vflux[:,:,:,:,4]
    a_vfluxzz = a_vflux[:,:,:,:,8]

    f_dsurfV = PObj["dsurf_vflux"]
    a_dsurfV = f_dsurfV.read(startrec,endrec)
    a_dsurfVxx = a_dsurfV[:,:,:,:,0]

    # instantaneous time
    time_axis = 3
    timeTot = a_psurf.shape[time_axis]
    timeVect = np.arange(0,timeTot)
    Gamma_inst = np.nan*np.ones(timeTot)
    
    res = np.nan*np.ones(len(timeVect)*2)
    res = res.reshape(len(timeVect),2)

    timecount = 0
    for time in timeVect:
        
        a_psurfxx_t = a_psurfxx[:,:,:,time] 
        a_psurfyy_t = a_psurfyy[:,:,:,time]
        a_psurfzz_t = a_psurfzz[:,:,:,time]

        a_vfluxxx_t = a_vfluxxx[:,:,:,time]
        a_vfluxyy_t = a_vfluxyy[:,:,:,time]
        a_vfluxzz_t = a_vfluxzz[:,:,:,time]

        a_dsurfVxx_t = a_dsurfVxx[:,:,:,time]


        # Global surface tension - Only spatial averaged (in y and z planes)
        dims = np.linspace(0,799,800) # 800 bins along x axis

        PxxG = np.nan*np.ones(800)
        PyyG = np.nan*np.ones(800)
        PzzG = np.nan*np.ones(800)
        
        ## averaging over y-z plane (mean of the y-z slice) at every bin along x
        for pidx in dims:
            PxxG[int(pidx)] = a_psurfxx_t[int(pidx)].mean() 
        for pidy in dims:
            PyyG[int(pidy)] = a_psurfyy_t[int(pidy)].mean()
        for pidz in dims:
            PzzG[int(pidz)] = a_psurfzz_t[int(pidz)].mean()
    
        vFluxxxG = np.nan*np.ones(800) 
        vFluxyyG = np.nan*np.ones(800) 
        vFluxzzG = np.nan*np.ones(800)
        
        for vidx in dims:
            vFluxxxG[int(vidx)] = a_vfluxxx_t[int(vidx)].mean()
        for vidy in dims:
            vFluxyyG[int(vidy)] = a_vfluxyy_t[int(vidy)].mean()
        for vidz in dims:
            vFluxzzG[int(vidz)] = a_vfluxzz_t[int(vidz)].mean()
 
        dsurfVxxG = np.nan*np.ones(800)
        for didx in dims:
            dsurfVxxG[int(didx)] = a_dsurfVxx_t[int(didx)].mean()

        pN = PxxG + vFluxxxG + dsurfVxxG 
        pTy = PyyG + vFluxyyG
        pTz = PzzG + vFluxzzG
        pT = 0.5 * (pTy + pTz)

        pNminusT = pN - pT 
        
        loLim = 380
        hiLim = 420
        Gamma = np.trapz(pNminusT[loLim:hiLim],locx[loLim:hiLim])

        # storing Gamma at each timesnap
        Gamma_inst = round(Gamma,5)
        
        res[timecount,:] = [time, Gamma_inst]

        msg = 'Global surface tension is: {} \n'
        print(msg.format(Gamma_inst))
        
        timecount = timecount + 1
        
    return res

    
def mrr_conver2binary (file,thresh):

    import matplotlib.pyplot as plt
    from PIL import Image

    image_file_name = file # should be a tif file, but could be anything
    col = Image.open(image_file_name)
    plt.figure()
    plt.imshow(col)
    plt.title('Original image')
    plt.show()
    
    # cast to gray scale
    gray = col.convert('L')

    conv_thres = thresh
    bw = gray.point(lambda x: 0 if x<conv_thres else 255, '1')
    plt.figure()
    plt.imshow(bw)
    plt.title('Converted image')
    plt.axis('off')
    plt.show()  
    
    return bw
    

def mrr_fractals (filename,thresh, imshow='False', background = 'Light'):
    
    """
    Returns fractal dimension (Hausdorff dimension) of an rgb image. Please DO a trial and error
    to set the threshold value. The value that can successfully identify the intended surface features
    should be used for final analysis. For images, .tif format is recommended.
    * optional argument: imshow if set to 'True', images (original, converted gray & binary) are shown
    * optional argument: background if set to 'Dark', image is inverted
    *** The fractal dimension is that of the black pixels. Therefore, it is essential to compare the original
    and binary images to understand whether the background needs be set to Dark or Light
    It is important to 
    
    Mazority of the fractal analysis algorith is adapted from: 
    https://www.physics.utoronto.ca/apl/fvf/python_code/box_count.py
    """

    import matplotlib.pyplot as plt
    from PIL import Image 
    import numpy as np
    import pylab as pl

    image =  Image.open(filename)  # create PIL image object
    image_grayscale = image.convert('L')   # cast the image into grey scale 0-255 (just in case)
    
    image_matrix = np.asmatrix(image_grayscale).copy()  # copy the greyscale values into numpy matrix format
    # binarize the image (it may be thresholded already)
    # 0 is black 255 is white
    if background == 'Dark':
        for i in range(image.size[1]):
            for j in range(image.size[0]):
                if image_matrix[i,j] > thresh:
                    image_matrix[i,j] = 0    # black
                else: 
                    image_matrix[i,j] = 255   # white
    else:    
        for i in range(image.size[1]):
            for j in range(image.size[0]):
                if image_matrix[i,j] > thresh:
                    image_matrix[i,j] = 255    # white
                else: 
                    image_matrix[i,j] = 0   # black

    # Make a list of black pixel coordinates

    pixels=[]
    for i in range(image.size[1]): 
        for j in range(image.size[0]):
            if image_matrix[i,j] == 0:  #  pixel is black 
                pixels.append((i,j))    #  count it
 
    Lx=image_grayscale.size[0]
    Ly=image_grayscale.size[1]
    
    pixels=pl.array(pixels)   # turn into a pylab array
 
    # computing the fractal dimension
    # considering only scales in a logarithmic list
    scales=np.logspace(1, 8, num=20, endpoint=False, base=2)
    Ns=[]
    # looping over several scales
    for scale in scales:
        #print ('....... Scale :', scale)
        # computing the histogram
        H, edges=np.histogramdd(pixels, bins=(np.arange(0,Lx,scale),np.arange(0,Ly,scale)))
        Ns.append(np.sum(H > 0))
 
    # linear fit, polynomial of degree 1 --> a line
    coeffs=np.polyfit(np.log(scales), np.log(Ns), 1)

    D = - coeffs[0]  #the fractal dimension is the negative of the slope of the line
 
    if imshow =='True':        
    
        plt.figure(frameon=False)
        plt.imshow(image) # an image file, not data file
        plt.title('Original Image')
        plt.axis('off')
        plt.show()
        
        plt.figure(frameon=False)
        plt.imshow(image_grayscale,cmap='gray', vmin = 0, vmax = 255)
        plt.colorbar()
        plt.title('Grayscale conversion')
        plt.axis('off') 
        plt.show()
    
        plt.figure(frameon=False)
        img = Image.fromarray(image_matrix)      	# save the binarizd image
        plt.imshow(img, cmap = 'gray', vmin=0, vmax=255)
        plt.colorbar()
        plt.axis('off') 
        tit1 = 'Binary conversion, threshold: {}'
        plt.title(tit1.format(thresh))
        #img.save(image_file_name + '_binarized.tif')  # as a check against the original image 
    
           
        plt.rcParams["figure.figsize"] = [9,10]
        plt.rcParams["figure.autolayout"] = True
        plt.rcParams['font.size'] = '34'
        plt.rcParams['font.family'] = 'Times New Roman' #'sans-serif'

        plt.figure(frameon=False)
        plt.plot(np.log(scales),np.log(Ns), 'o', mfc='none', label='data')
        plt.plot(np.log(scales), np.polyval(coeffs,np.log(scales)),label='fit')
        plt.xlabel('$\log \epsilon$')
        plt.ylabel('$\log N$')
        
        B = [np.log(scales),np.log(Ns)]
        #plt.xlim([0,5])
        #plt.ylim([4,15])
        #pl.title(image_file_name +' D = '+str(D))
        #pl.legend()
        #pl.savefig(image_file_name+'_powerlaw.pdf')
 
        print('The Hausdorff dimension is ', D )
        #np.savetxt(image_file_name+'_data.txt', (scales,Ns))
    
    return  D 

    
    
    
    
    
    