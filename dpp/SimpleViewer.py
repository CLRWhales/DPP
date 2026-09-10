#this is a simplt hdf5 das viewer function:
from dpp.simpleDASreader4 import load_DAS_file
from scipy.signal import detrend, resample, butter, sosfiltfilt, decimate, resample_poly
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import os
import argparse

def simpleViewer():
    parser = argparse.ArgumentParser(description="Process a filename or directory.")
    parser.add_argument("path", nargs="?", type=str, help="The name of the file to process")
    args = parser.parse_args()

    signal,meta = load_DAS_file(args.path)

    sos = butter(N = int(6),
                        Wn = 10,
                        btype = 'highpass',
                        fs = 1/meta['header']['dt'],
                        output = 'sos')

    data = sosfiltfilt(sos = sos,
                    x = signal,
                    axis = 0)

    max = np.percentile(np.abs(data),q=75)
    norm = colors.TwoSlopeNorm(vmin=-max, vcenter=0, vmax=max)
    plt.figure()
    plt.imshow(data,aspect = 'auto',cmap = 'seismic',norm = norm)
    plt.colorbar()
    plt.xlabel('channel #')
    plt.ylabel('time sample')
    plt.title(os.path.basename(args.path))
    plt.show()
def main():
  simpleViewer()
if __name__ == '__main__':
  main()