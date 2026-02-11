#%% this script writes a default processing ini file with the config parser function that can be modified to process based on needs
import configparser

config = configparser.ConfigParser()
config['DataInfo'] = {'Directory':'your/directory/here', #directory to where the data lives
                      'n_files' : '6', #how many files do you want concatenated and processed per batch
                      'n_workers': '4', #how many batches do you want to do simultaneously, note if too high, can overload available memory
                      'auto_optimize':'true' #auto pick the number of workers. will take the minimum between n_workers and estimated value
                      }
config['ProcessingInfo'] = {'n_synthetic':'auto', # number of synthetic receiver positions to be used along the fiber, auto means fill it up acording to synthetic spacing, meters means auto, but spacing is in units of meters and auto fills the full extent of the fiber.
                            'synthetic_spacing':'250', #number of channels between the start of each synthetic receiver
                            'c_start':'0', #where to start the channel idx along the fiber? idx or meters as described above
                            'c_end':'',
                            'n_stack':'5', #how many chanels to use within one synthetic receiver
                            'stack':'true', #do you want to stack these channels
                            'stackaxis':'1', #along which axis do you want to stack, keep to 1 to avoid bugs
                            'fs_target' : 'auto', #resample the data in time to this, auto moves it to the highest power of 2
                            'starttime' : '', #optional, in the format HHMMSS if you only want to process some of the data
                            'stoptime': '', # optional, in the format HHMMSS, if you only want to process some of the data
                            'verbose': 'false', # diagnositic information during the run, troubleshooting only
                            'cmn_filt': 'false', #this is a place holder for common mode filtering toggle, it does nothing.
                            'unwr':'false', #do you want to phase unwrap the signal through space 
                            'integrate':'true' #do you want to time intrgrate the signal
                            }
config['FFTInfo'] = {'input_type':'time', #unit for fft parameter definition (time, point)
                     'n_fft':'0.5', #if time, Hz resoltuion of spectrogram, if points, n points to include, will add zero padding if necessary
                     'n_overlap':'0.5', #if time, seconds of window overlap, if points, n points to overlap
                     'n_samp':'0.5', #if time, how may seconds of data to include, if points, how many points to use
                     'do_fk':'false' #do you want to do 1D FFT or 'FK'
                     }
config['SaveInfo'] = {'directory':'your/directory/here', #general directory to save the data within
                      'run_name':'test', # run name to save within the above directory
                      'data_type':'cleaning' # what type of data do you want output? magnitude is a 3d FTX spectrum, complex is the complex FTX, and cleaning is data used for cleaning in the cleaner app
                      }
config['FilterInfo'] = {'type':'highpass', #can be highpass, lowpass, bandpass, bandstop, none
                        'lowcut':'0.5', #defines the lower cut of the filter (Hz), if lowpass, highpass are used, only this value is necessary
                        'highcut':'-1', #defines the upper cut of the fuilter (Hz), if lowpass, highpass, this is ignored
                        'order':'6', #filter order
                        'zerophase':'0', #filter zerophase
                        'alphataper':'0.1' #curently doesnt do anything
                        }
config['FKInfo'] = {'nfft_time':'512', #parameters in units of samples for the fk transform (how many time points) (applied after resampling)
                    'nfft_space':'512', #how many channel to include, both should be a power of 2 for speed (applied after any stacking)
                    'rescale':'true',
                    'overlap':'2',
                    'fold':'true', #do you want to fold your Fks over around the 0 wavenumber
                    'vmin':'1470', #dod oyu want to velocity filter your saved fks? set to -1 for no minimum vel
                    'vmax':'3500', #do you want to velocity filter your saved fks? set to -1 for no max vel
                    'fmin':'0', #do you want to freq filter your saved FKS? set to -1 for no min frequency 
                    'fmax':'80', #do you want to freq filter your saved FKS? set to -1 for no max frequency 
                    'thresh': '5',
                    'sample_method':'none', #do you want to sample other non signals? none for no, same for a random selection of the same number, in files with no water band, randomly pull n images for more diverse training material
                    'n':'50'
                    }

def main():
  with open('example.ini', 'w') as configfile:
    config.write(configfile)
if __name__ == '__main__':
  main()