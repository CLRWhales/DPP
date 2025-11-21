#this script contains the helper functions developed by calder robinson
#%%
import numpy as np
from concurrent.futures import ThreadPoolExecutor
import glob
import dpp.h5pydict as h5pydict
import os
from numpy.lib.stride_tricks import as_strided
import scipy.ndimage as NDI
from scipy.stats import entropy


def faststack(X,n, wind = 1):
     """ fast adjacent channel stack through time
    Parameters
    ----------
    X : 2d np.array 
        2d data array with dimension (T,X).
    n : int
        number of channels to stack together
    wind : 1d np.array
        optional cross channel windowshape, default is square

    Returns
    -------
    stacked : 2d np.array
        data array containing the mean of n adjacent channels 
    
    Note: when using a window other than 1, if n does not divide into the number of X, the final group  average will be returend without the window applied, 
    and it contains fewerthan requested channels.
     """
     rows,cols = X.shape
     trimval = np.mod(cols,n)

     if wind == 1:
        wind = np.ones(n)
     elif len(wind) != n:
         raise ValueError("window not same size as stack number.")
         

     if trimval!=0:
          trimmedmean = X[:,-trimval:].mean(axis = 1)
          X = X[:,:-trimval]
          stacked = X.reshape(-1,n)
          stacked = np.average(stacked, axis = 1, weights = wind).reshape(rows,-1)
          stacked = np.c_[stacked,trimmedmean]
     else:
          stacked = X.reshape(-1,n)
          stacked = np.average(stacked, axis = 1, weights=wind).reshape(rows,-1)

     return stacked

def loadFTX(directory,nworkers = 1,start_f=None, stop_f=None, start_cidx=None, stop_cidx=None): 
    """
    This is a function that can be used to quickly load FTX files in a directory. if there are too many files/too large, can overload the ram. 

    inputs:
    X: directory of files, without trailing slash
    nworkers: how many workers do you want to deploy?
    startf: index of frequency start for slicing
    stopf: index of frequency stop for slicing
    start_cidx: index of channel start for slicing
    stop_cidx: index of channel stop for slicing


    returns:
    data: concatenated np array of data in the directory of dimension FTX
    freqs: np. array of freqs for the data
    times: np array of relative time since the start of the array, ignores gaps
    channels: np array of channel distances of the array

    TODO:
    build in gap handelling in the time index by parsing the file names?

    """

    tmpdir = directory
    for i in range(4):
        clist = glob.glob(os.path.join(tmpdir, 'Dim_Channel.txt'))
        if clist:
            break
        else:
            tmpdir = os.path.split(tmpdir)[0]

    if i == 4:
        raise ValueError("could not find channel, frequency, or time files within 4 levels, check directory")

    channels = np.loadtxt(os.path.join(tmpdir, 'Dim_Channel.txt'))[start_cidx:stop_cidx]
    times =  np.loadtxt(os.path.join(tmpdir, 'Dim_Time.txt'))
    freqs = np.loadtxt(os.path.join(tmpdir, 'Dim_Frequency.txt'))[start_f:stop_f]

    filepaths = sorted( glob.glob(os.path.join(directory, '*.npy')))#[0:10]
    if len(filepaths) == 0:
        raise ValueError("no files found at directory.")
    
    lt = len(times)
    bt = np.arange(len(times),dtype= 'int')

    data = np.empty((len(freqs),len(times)*len(filepaths)+1,len(channels)))

    with ThreadPoolExecutor(max_workers=nworkers) as exe:
        futures = exe.map(np.load,filepaths)
        for i, ff in enumerate(futures):
            tidx = i*lt + bt
            data[:,tidx,:] = ff[start_f:stop_f,:,start_cidx:stop_cidx]

    times = np.arange(start = 0, stop = data.shape[1])*0.25


    return data, freqs,times,channels

def load_meta(filename,metaDetail = 1):
    """
    Extracted from ASN load das file, see for details
    """
    with h5pydict.DictFile(filename,'r') as f:
        # Load metedata (all file contents except data field)
        m = f.load_dict(skipFields=['data']) 
        if metaDetail==1:
            ds=m['demodSpec']
            mon=m['monitoring']
            meta=dict(fileVersion = m['fileVersion'],
                      header     = m['header'],
                      timing     = m['timing'],
                      cableSpec  = m['cableSpec'],
                      monitoring = dict(Gps = mon['Gps'],
                                        Laser=dict(itu=mon['Laser']['itu'])),
                      demodSpec  = dict(roiStart = ds['roiStart'],
                                        roiEnd   = ds['roiEnd'],
                                        roiDec   = ds['roiDec'],
                                        nDiffTau = ds['nDiffTau'],
                                        nAvgTau  = ds['nAvgTau'],
                                        dTau     = ds['dTau']))
        else:
            meta = m
    return meta

def window_rms(a, windowsize):
    a2 = np.square(a)
    window = np.ones(windowsize)/float(windowsize)
    return np.sqrt(np.convolve(a2, window, 'valid'))

def compute_entropy(arr):
    '''
    computes the local timewise spectral entropy of FTX after prewhitening with the mean
    inputs
        arr: FTX np.array
    outputs:
        entropy: 2d np array
        spectral entropy through time, relative to the timewise means of the FTX block, same dimension of TX
    '''
    arr = abs(arr)
    arr /= np.mean(arr, axis = 1)[:,None,:] #prewhiten
    arr /=np.sum(arr,axis = 0) #normalize to psd
    denom = np.log2(arr.shape[0]-1)
    num = np.sum(arr * np.log2(arr), axis = 0)
    entropy = -num/denom
    return entropy

def compute_FK_speed(k, f):
    """
    Computes angle and slope from a reference cell to every other cell in a 2D grid.

    Parameters:
        k (2D np.array): array of  of x positions.
        f (2D np.array): Grid of y positions.

    Returns:
        angle_grid (2D np.array): Angles in radians.
        slope_grid (2D np.array): Slopes (np.inf for vertical lines).
    """
    x,y = np.meshgrid(k,f)
    angle_grid = np.arctan2(y, x)

    with np.errstate(divide='ignore', invalid='ignore'):
        slope_grid = np.true_divide(y, abs(x))
        slope_grid[x == 0] = np.inf  # handle vertical lines

    return  slope_grid,angle_grid


def average_values_by_slope_bins(slope_grid, value_grid, slope_bins):
    """
    Computes average of values from value_grid based on slope bins.

    Parameters:
        slope_grid (2D np.array): Array of slope values.
        value_grid (2D np.array): Array of values to average.
        slope_bins (1D np.array): Array defining slope bin edges.

    Returns:
        averages (1D np.array): Average of values in each slope bin.
        counts (1D np.array): Number of values in each bin.
    """
    slopes_flat = slope_grid.flatten()
    values_flat = value_grid.flatten()

    bin_indices = np.digitize(slopes_flat, slope_bins) - 1  # 0-based bins

    # Mask for valid values (finite slopes and values)
    mask = np.isfinite(slopes_flat) & np.isfinite(values_flat)
    slopes_flat = slopes_flat[mask]
    values_flat = values_flat[mask]
    bin_indices = bin_indices[mask]

    # Sum of values and count per bin
    bin_sums = np.bincount(bin_indices, weights=values_flat, minlength=len(slope_bins) - 1)
    bin_counts = np.bincount(bin_indices, minlength=len(slope_bins) - 1)

    # Avoid divide by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        averages = np.true_divide(bin_sums, bin_counts)
        averages[bin_counts == 0] = np.nan  # Or set to 0, or leave as NaN

    return averages

def foldFK(fk,K):
    negs = np.where(K<0)[0]
    pos = np.where(K>0)[0]
    negV = fk[:,negs]
    posV = fk[:,pos]
    negV = np.fliplr(negV)
    posV = np.append(posV,negV[:,-1:],axis = 1)
    return np.stack([posV,negV],axis = -1)

def KL_div(imglist):
    mi = np.min(imglist)
    ma = np.max(imglist)
    total,_ = np.histogram(np.concat(imglist).flatten(),bins = 128,range=(mi,ma),density = True)
    total+=1e-15
    total /= np.sum(total)
    hists= [np.histogram(im.flatten(),bins = 128,range = (mi,ma), density = True) for im in imglist]
    hists = [(h[0]+1e-15)/np.sum((h[0]+1e-15)) for h in hists]
    ents = [np.round(entropy(fk,total),6)for fk in hists]
    return ents

# def sliding_window_FK(arr, window_shape,overlap = 2, rescale = False, fold = True):
#     step_y, step_x = window_shape[0] // overlap, window_shape[1] // overlap
#     shape = (
#         (arr.shape[0] - window_shape[0]) // step_y + 1,
#         (arr.shape[1] - window_shape[1]) // step_x + 1,
#         window_shape[0],
#         window_shape[1]
#     )
#     strides = (
#         arr.strides[0] * step_y,
#         arr.strides[1] * step_x,
#         arr.strides[0],
#         arr.strides[1]
#     )
    
#     windows = as_strided(arr, shape=shape, strides=strides)
#     weight = np.outer(np.hanning(window_shape[0]),np.hanning(window_shape[1]))

#     velweights = np.hanning(7)
#     velweights = velweights/np.sum(velweights)

#     ks = np.fft.fftshift(np.fft.fftfreq(512,12))
#     freqs = np.fft.rfftfreq(512,1/256)[1:]
#     slope,_ = compute_FK_speed(ks,freqs)
#     slope[slope > 10000] = 10000

#     slopebins = np.linspace(0,10000,num = 101)
#     centers = 0.5* (slopebins[:-1] + slopebins[1:])
#     slopebins[0] = -np.inf
#     slopebins[-1] = np.inf

#     #acumulators
#     results = []
#     pos = []
#     maxs = []
#     Ls = []
#     vels = []
#     entfull = []
#     flags = []

#     # c_min = 1000
#     # c_max = 5000
#     # f = np.fft.rfftfreq(512, d = 1/256)[1:]
#     # k1 = np.fft.fftshift(np.fft.fftfreq(512, d = 12))
#     # kk,ff = np.meshgrid(k1,f)
#     # slow = ff<np.abs(kk*c_min)
#     # fast = ff > np.abs(kk*c_max)
#     # full = slow+fast

#     for j in range(shape[1]): #range
#         map = np.zeros(shape = (int(window_shape[0]/2),window_shape[1]))
#         N = 0
#         intermediate = []
#         for i in range(shape[0]): #time
#             win = windows[i, j]*weight
#             fft_result = np.fft.fftshift(np.abs((np.fft.rfft2(win,axes=(-1,-2),norm = 'forward'))),axes=1)[1:,:]#this removes time DC offset
#             pos.append((i * step_y, j * step_x))
#             map = map+fft_result
#             intermediate.append(fft_result)
#             N = N+1
        
#         mean_img = map/N
#         stdev = np.std(intermediate)

#         for f in intermediate:
#             tmp= NDI.gaussian_filter((f-mean_img)/stdev,sigma=(1))
#             peak_lock = np.unravel_index(np.argmax(tmp), tmp.shape)
#             fmax=peak_lock[0]
#             kmax = np.round(np.abs(peak_lock[1]-window_shape[1]/2))
#             L = int((peak_lock[1]-window_shape[1]/2)<0)
#             Ls.append(L)
#             maxs.append((fmax,kmax))
#             vels.append(np.round(slope[peak_lock]))

#         # tmp = [fk[~full].flatten() for fk in intermediate]
#         # ents = KL_div(20*np.log10(tmp))
#         # thresh = np.median(ents)+2.5*np.std(ents)
#         # flag = np.zeros_like(ents)
#         # flag[ents>thresh] = 1
#         # entfull.extend(ents)
#         # flags.extend(flag)
        
#         results.extend(20*np.log10(intermediate))

#     if rescale:
#         vals = np.stack(results,axis=0)[:,128:,:]
#         #print(vals.shape)
#         low = np.floor(np.percentile(vals,1)) #file wise
#         high = np.ceil(np.percentile(vals,99)) #filewise
#         del vals
#         results = [(255*((r-low)/(high-low))).clip(0, 255).astype(np.uint8) for r in results]

#     if fold:
#         results = [foldFK(r,ks) for r in results] #this folds the fk so pos and neg ks are in separate image channels
    
#     outputs = [(k[0],k[1],k[2],k[3],k[4]) for k in zip(results,pos,maxs,vels,Ls)]

#     return outputs

def sliding_window_FK(arr, window_shape,overlap = 2, rescale = False, fold = True):
    step_y, step_x = window_shape[0] // overlap, window_shape[1] // overlap
    shape = (
        (arr.shape[0] - window_shape[0]) // step_y + 1,
        (arr.shape[1] - window_shape[1]) // step_x + 1,
        window_shape[0],
        window_shape[1]
    )
    #print(shape)
    strides = (
        arr.strides[0] * step_y,
        arr.strides[1] * step_x,
        arr.strides[0],
        arr.strides[1]
    )
    #print(strides)
    windows = as_strided(arr, shape=shape, strides=strides)
    weight = np.outer(np.hanning(window_shape[0]),np.hanning(window_shape[1]))

    ks = np.fft.fftshift(np.fft.fftfreq(512,12))
    freqs = np.fft.rfftfreq(512,1/256)[1:]
    slope,_ = compute_FK_speed(ks,freqs)
    slope[slope > 10000] = 10000

    #acumulators
    results = []
    pos = []
    maxs = []
    Ls = []
    vels = []


    for j in range(shape[1]): #range
        map = np.zeros(shape = (int(window_shape[0]/2),window_shape[1]))
        N = 0
        intermediate = []
        for i in range(shape[0]): #time
            win = windows[i, j]*weight
            fft_result = np.fft.fftshift(np.abs((np.fft.rfft2(win,axes=(-1,-2),norm = 'forward'))),axes=1)[1:,:]#this removes time DC offset
            pos.append((i * step_y, j * step_x))
            #map = map+fft_result
            intermediate.append(fft_result)
            #N = N+1
        mean_img = np.mean(intermediate, axis = 0)
        #mean_img = map/N
        stdev = np.std(intermediate)

        for f in intermediate:
            tmp= NDI.gaussian_filter(f-mean_img/stdev,sigma=(1))
            peak_lock = np.unravel_index(np.argmax(tmp), tmp.shape)
            fmax=peak_lock[0]
            kmax = np.round(np.abs(peak_lock[1]-window_shape[1]/2))
            L = int((peak_lock[1]-window_shape[1]/2)<0)
            Ls.append(L)
            maxs.append((fmax,kmax))
            vels.append(np.round(slope[peak_lock]))

        intermediate = 20*np.log10(intermediate)
        mintermediate = np.mean(intermediate, axis = 0)
        results.extend((intermediate))#-mintermediate))

    if rescale:
        vals = np.stack(results,axis=0)[:,128:,:]
        #vals[vals<0] = 0
        #print(vals.shape)
        low = np.floor(np.percentile(vals,1)) #file wise
        high = np.ceil(np.percentile(vals,99)) #filewise
        del vals
        results = [(255*((r-low)/(high-low))).clip(0, 255).astype(np.uint8) for r in results]

    if fold:
        results = [foldFK(r,ks) for r in results] #this folds the fk so pos and neg ks are in separate image channels
    
    outputs = [(k[0],k[1],k[2],k[3],k[4]) for k in zip(results,pos,maxs,vels,Ls)]

    return outputs


