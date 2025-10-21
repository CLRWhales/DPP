# this function seeks to process DAS FFT dat a s quickly as possible. GPU from cupy should be a drop in replacement
import numpy as np; import time
#from numba import njit

#build in a flag that toggles the zero padding on and off.
def reshape_array_with_overlap(arr, N, N_overlap):
    """ reshapes a typical das array with overlap so that fft can be applied efficiently.

    Parameters
    ----------
    arr : 2d np.array
        2d data array.
    N : int
        number samples in a section
    N_overlap : int
        number of samples to overlap.

    Returns
    -------
    reshaped_arr : 2d np.array
        data array sliced to be able to apply windowing and fft very quickly along 0th axis.
    """
    rows, cols = arr.shape
    
    # Validate that overlap doesn't exceed block size
    if N_overlap >= N:
        raise ValueError("Overlap cannot be greater than or equal to the block size.")
    
    # Split the array into blocks with the specified overlap
    reshaped = []
    start_idx = 0
    while start_idx + N <= rows:
        reshaped.append(arr[start_idx:start_idx + N, :])
        start_idx += N - N_overlap  # Move the starting index to create overlap
    
    # Check if there's a remaining block that is less than N rows
    if start_idx < rows:
        # Calculate how much padding is needed
        remaining_rows = rows - start_idx
        padded_block = np.vstack([
            arr[start_idx:rows, :],  # Remaining rows
            np.zeros((N - remaining_rows, cols))  # Pad the block to size N
        ])
        reshaped.append(padded_block)


    # Concatenate the blocks along the columns axis (axis=1)
    reshaped_arr = np.concatenate(reshaped, axis=1)
    
    return reshaped_arr


def sneakyfft(X,N_samp,N_overlap,N_fft, window,fs):
    """ computes multi channel spectrogram very fast.

    Parameters
    ----------
    X : 2d np.array
        2d data array. (time, channel)
    N_samp : int
        number samples in a spectrogram frame
    N_overlap : int
        number of samples to overlap for each frame.
    N_fft : int
        fft size (samples) to zero pad to.
    Window : np.array
        the output of np window function or equivalent, must be same lenth as N_samp for broadcasting
    fs : int
        sample rate (Hz) of the das data series'

    Returns
    -------
    spec : 3d np.array, complex
        stack of spectrograms for each channel included in the input X, dimension of freq, time, channel 
    f : np.array
        vector of frequencies of the spectrogram array
    t: np.arry
        vector of relative times of the time bins of the spectrogram.
    """
    if N_overlap >= N_samp:
        raise ValueError("Overlap cannot be greater than or equal to the block size.")
    
    if len(window) != N_samp:
        raise ValueError("Window and block size must be the same.")
    
    reshaped = reshape_array_with_overlap(X,N_samp,N_overlap)
    reshaped = reshaped * window[:,None]
    fft_out = np.fft.rfft(reshaped, n = N_fft, axis =0)
   
    nt_slices = fft_out.shape[1]//X.shape[1]
    spec = fft_out.reshape(fft_out.shape[0],nt_slices,X.shape[1])
    #f = fs/N_FFT*np.arange(fs+1)
    f = np.fft.rfftfreq(N_fft,1/fs)
    t = np.arange(nt_slices)*(N_samp-N_overlap)/fs


    return spec, f, t




