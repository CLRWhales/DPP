#this script is the first steps into a consistent processing pipeline for bulk data runs based on an ini file io
import os
import glob,time
import numpy as np 
from scipy.signal import detrend, resample, butter, sosfiltfilt
from simpleDASreader4 import load_DAS_file, unwrap, combine_units #nned this if the other functions are uncommented
from DASFFT import sneakyfft
import configparser
import argparse
import Calder_utils as Calder_utils
import math
import datetime
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from functools import partial
import imageio
#handling lack of tk on some linux distros. shoddy, need to fix in the future
try:
    import tkinter as tk
except ImportError:
    available = False
else:
    available = True
    from tkinter import filedialog

def load_INI():
    """
    Load a .ini from argument or file browser.
    Allows selecting a single file or a folder containing multiple .ini files.
    """
    parser = argparse.ArgumentParser(description="Process a filename.")
    parser.add_argument("filename", nargs="?", type=str, help="The name of the file to process")
    args = parser.parse_args()

    if args.filename:
        config = configparser.ConfigParser()
        config.read(args.filename)
        return config

    root = tk.Tk()
    root.title("Select .ini file or folder")
    result = []

    def select_file():
        path = filedialog.askopenfilename(
            title="Select a .ini",
            filetypes=[("INI files", "*.ini")]
        )
        if path:
            result.append(path)
            root.quit()

    def select_folder():
        folder = filedialog.askdirectory(title="Select folder with .ini")
        if folder:
            ini_files = [
                os.path.join(folder, f)
                for f in os.listdir(folder)
                if f.lower().endswith(".ini")
            ]
            ini_files.sort()
            if ini_files:
                result.extend(ini_files)
            root.quit()

    btn_file = tk.Button(root, text="Select a .ini", command=select_file)
    btn_file.pack(pady=10)
    btn_folder = tk.Button(root, text="Select a folder", command=select_folder)
    btn_folder.pack(pady=10)

    root.mainloop()
    root.destroy()
    return result if result else None



def load_file(channels, verbose, filepath):
    """
    load a single DAS file -> from kevinG
    """
    if verbose: 
        fid = os.path.basename(filepath);  print(f"/nloading file {fid}")
    data, meta = load_DAS_file(filepath, chIndex=channels, roiIndex=None, samples=None,
                      integrate=False, unwr=False, metaDetail=1, useSensitivity=False,
                      spikeThr=None)
    return data, meta


def load_files(path_data, channels, verbose, fileIDs):
    """
    distribute multiple file loading over threads -> from kevinG
    """
    # create a thread pool
    with ThreadPoolExecutor(len(fileIDs)) as exe:
        
        file_paths = [os.path.join(path_data, str(fid) + '.hdf5') for fid in fileIDs]
        # load files
        results = exe.map(partial(load_file,channels, verbose), file_paths)
        # collect data
        list_data=[]; list_meta=[];
        for listx in results: 
            list_data.append(listx[0]); list_meta.append(listx[1])

        return (list_data,list_meta, fileIDs)



def preprocess_DAS(data, list_meta, unwr=True, integrate=True, useSensitivity=True, spikeThr=None): 
    """
    preprocess loaded DAS data, for details see SimpleDASreader4  -> from kevinG
    """
    meta = list_meta[0]
    # pre-process raw data
    if unwr or spikeThr or integrate:
        if meta['header']['dataType']<3 or meta['demodSpec']['nDiffTau']==0:
            raise ValueError('Options unwr, spikeThr or integrate can only be\
                             used with time differentiated phase data')
    if unwr and meta['header']['spatialUnwrRange']:
        data=unwrap(data,meta['header']['spatialUnwrRange'],axis=1)
    unit=meta['appended']['unit']
    #print(f"DEBUG: unit={unit}")
    if spikeThr is not None:
        data[np.abs(data)>spikeThr] = 0
    if useSensitivity:
        if 'sensitivity' in meta["header"].keys():
            data/=meta['header']['sensitivity']
            sensitivity_unit=meta['header']['sensitivityUnit']
        elif 'sensitivities' in meta["header"].keys(): 
            data/=meta['header']['sensitivities'].item()
            sensitivity_unit=meta['header']['sensitivityUnits'].item()
        else: 
            raise KeyError("!! no sensitivity keys in meta dict")
        unit=combine_units([unit, sensitivity_unit],'/')
        #print(f"DEBUG: sensitivity_unit={sensitivity_unit}, unit={unit}")
    if integrate:
        data=np.cumsum(data,axis=0)*meta['header']['dt']
        unit=combine_units([unit, unit.split('/')[-1]])
        #print(f"DEBUG: unit={unit}")
    meta['appended']['unit']=unit
   
    #update all metas
    for metax in list_meta: 
        metax['appended']['unit']=unit
        
    return (data, list_meta)


def LPS_block(path_data,channels,verbose,config, fileIDs):
    """
    Load, Process, and Save, a single block of das data files for further analysis and visualization 
    """
    #load into list
    do_fk = config['FFTInfo'].getboolean('do_fk')

    # if do_fk:
    #     FKchans = channels
    #     channels = None

    list_data, list_meta, _ = load_files(path_data = path_data,
                                      channels = channels,
                                      verbose = verbose,
                                      fileIDs= fileIDs)
    data =  np.concatenate(list_data, axis=0)
    data, list_meta = preprocess_DAS(data, list_meta)

    # if do_fk:
    #     data = data[:,FKchans]

    data /= 1E-9 #scaling into units of strain is handled, this moves it to nano strain? 

    if config['ProcessingInfo'].getboolean('cmn_filt'):
        data-=np.median(data,axis = 1)[:,None]
    #do stacking
    n_stack = int(config['ProcessingInfo']['n_stack'])
    chans = list_meta[0]['appended']['channels']

    if config['ProcessingInfo'].getboolean('stack'):
        data = Calder_utils.faststack(data,n_stack)
        chIDX = np.arange(start = 0,stop = len(chans), step = n_stack)
    else:
        chIDX = np.arange(start = 0,stop = len(chans))

    #do resampling
    dt = list_meta[0]['header']['dt']
    dx = list_meta[0]['appended']['channels'][chIDX][1] - list_meta[0]['appended']['channels'][chIDX][0]
    if config['ProcessingInfo']['fs_target'] == 'auto':
        fs_target = 2**math.floor(math.log(1/dt,2))
    else:
        fs_target = int(config['ProcessingInfo']['fs_target'])

    num = data.shape[0]/(1/fs_target)*dt
    num = num.astype(int)

    #perfomring some lowpass AA filtering
    sos1 = butter(N = 6,Wn = fs_target/2, btype='low', fs = 1/dt, output= 'sos')
    data = sosfiltfilt(sos1,data,axis = 0)
    data = resample(data,num ,axis = 0);
    data=detrend(data, axis=0, type='linear');
    dt_new = 1/fs_target

    #filtering
    cuts = [float(config['FilterInfo']['lowcut']),float(config['FilterInfo']['highcut'])]
    dofilt = True
    

    match config['FilterInfo']['type']:
        case 'lowpass':
            cuts = cuts[0]
        case 'highpass':
            cuts = cuts[0]
        case 'bandpass':
            cuts = cuts
        case 'bandstop':
            cuts = cuts
        case 'none':
            dofilt = False
        case _:
            TypeError('input must be either "lowpass", "highpass","bandpass","bandstop", or "none')
    
    
    
    if dofilt:
        sos = butter(N = int(config['FilterInfo']['order']),
                    Wn = cuts,
                    btype = config['FilterInfo']['type'],
                    fs = fs_target,
                    output = 'sos')

        data = sosfiltfilt(sos = sos,
                        x = data,
                        axis = 0)
        
    
    if do_fk:
        times = np.arange(data.shape[0])*dt_new
        windowshape = (int(config['FKInfo']['nfft_time']),int(config['FKInfo']['nfft_space']))
        overlap = int(config['FKInfo']['overlap'])
        rsbool = config['FKInfo'].getboolean('rescale')
        fold= config['FKInfo'].getboolean('fold')
        fks = Calder_utils.sliding_window_FK(data,windowshape,overlap,rsbool,fold) 
        del data
        freqs = np.fft.rfftfreq(n=windowshape[0],d=dt_new)
        WN = np.fft.fftshift(np.fft.fftfreq(n=windowshape[1],d=dx))
        savetype = 'FK'
        # chIDX = FKchans
    else:
    # STFT 
        match config['FFTInfo']['input_type']:
            case 'point':
                N_samp= int(config['FFTInfo']['n_samp'])
                N_overlap = int(config['FFTInfo']['n_overlap'])
                N_fft = int(config['FFTInfo']['n_fft'])
                
            case 'time':
                N_samp= int(fs_target*float(config['FFTInfo']['n_samp']))
                N_overlap = int(N_samp*float(config['FFTInfo']['n_overlap']))
                N_fft = int(fs_target/float(config['FFTInfo']['n_fft']))
                
            case _:
                raise TypeError('input must be either "point" or "time", if time, make sure divisions yield power of 2 for speed')
            
        window = np.hanning(N_samp)
        spec, freqs, times = sneakyfft(data,N_samp,N_overlap,N_fft, window,fs_target)
        savetype = config['SaveInfo']['data_type']
        del data


    #saving info
    date = list_meta[0]['header']['time']
    fdate = datetime.datetime.fromtimestamp(int(date),tz = datetime.timezone.utc).strftime('%Y%m%dT%H%M%S')
    odir = config['Append']['outputdir']

    if fileIDs[0] == config['Append']['first']:
        freqname = os.path.join( odir, 'Dim_Frequency.txt')
        np.savetxt(freqname,freqs)

        channeldistance = list_meta[0]['appended']['channels'][chIDX]
        channelname= os.path.join(odir ,'Dim_Channel.txt')
        np.savetxt(channelname,channeldistance)

        timename = os.path.join(odir, 'Dim_Time.txt')
        np.savetxt(timename,times)

        if savetype=='FK':
            wnName = os.path.join(odir , 'Dim_Wavenumber.txt')
            pname = os.path.join(odir,'Dim_max.txt')
            posmax = fks[-1][1] + (fs_target,)
            np.savetxt(pname,posmax)
            np.savetxt(wnName,WN)

        cfgname = os.path.join(odir , 'config.ini')
        with open(cfgname, 'w') as configfile:
            config.write(configfile)

        
    match savetype: 
        case 'magnitude':
            magdir = os.path.join(odir , 'Magnitude')
            #Path(magdir).mkdir()
            os.makedirs(magdir, exist_ok=True)
            np.abs(spec,out = spec)
            np.log10(spec,out=spec)
            spec*=10
            #spec = 10*np.log10(abs(spec))
            fname = 'FTX' + str(fs_target) + '_' + fdate +'Z'
            data_name = os.path.join(magdir,fname)
            np.save(data_name,np.real(spec))

            
        case 'complex':
            compdir = os.path.join(odir , 'Complex')
            #Path(compdir).mkdir(exist_ok=True)
            os.makedirs(compdir, exist_ok=True)
            #spec = 10*np.log10(spec)
            fname = 'FTX' + str(fs_target) + '_' + fdate +'Z'
            data_name = os.path.join(compdir,fname)
            np.save(data_name,spec)

        case 'cleaning':
            #these thresholds are based on Robins values.
            mfreq = np.max(freqs)
            cuts = [0.5,5,30,130,mfreq+1] 
            for (l,h) in zip(cuts[:-1],cuts[1:]):
                if l > mfreq:
                    break
                f_idx = np.where(np.logical_and(freqs >= l,freqs <= h))
                TX = 10*np.log10(np.mean(abs(spec[f_idx[0],:,:]),axis = 0))
                cleandir = os.path.join(odir , str(l) + 'Hz_' + str(h) + 'Hz')
                #Path(cleandir).mkdir(exist_ok=True)
                os.makedirs(cleandir, exist_ok= True)
                fname = 'TX' + str(fs_target) + '_' + fdate +'Z'
                fout = os.path.join(cleandir, fname)
                np.save(fout,TX)

        case 'LTSA':
            LTSAdir = os.path.join(odir , 'LTSA')
            os.makedirs(LTSAdir, exist_ok=True)
            LTSA = np.mean(abs(spec),axis = 1)
            fname = 'FX_LTSA' + str(fs_target) + '_' + fdate +'Z'
            data_name = os.path.join(LTSAdir,fname)
            np.save(data_name,LTSA)
        
        case 'entropy':
            ENTdir = os.path.join(odir , 'Entropy')
            os.makedirs(ENTdir,exist_ok=True)
            np.abs(spec,out = spec)
            entropy = Calder_utils.compute_entropy(spec)
            fname = 'ENT' + str(fs_target) + '_' + fdate + 'Z'
            data_name = os.path.join(ENTdir,fname)
            np.save(data_name,entropy)

        case 'FK':
            FKDir = os.path.join(odir , 'FK',fdate + 'Z')
            os.makedirs(FKDir,exist_ok=True)
            for fk in fks:
                fname = 'FS'+str(fs_target)+'_T'+ str(fk[1][0]) + '_X' + str(fk[1][1]) + '_F' + str(fk[2][0]) + '_K' +str(fk[2][1]) +'_V'+ str(fk[3])+ '_L'+str(fk[4]) + '_' + fdate + 'Z.png'
                data_name = os.path.join(FKDir,fname)
                #print(fk[1].dtype)
                imageio.imwrite(data_name,fk[0])
                
            
        case _:
            raise TypeError('input must be either "magnitude", "complex","cleaning","LTSA","Entropy", or doFk must be set to true')
        

def main(config_path=None):
    # Load config
    import configparser
    config = configparser.ConfigParser()
    config.read(config_path)

    if config_path:
        config.read(config_path)
    else:
        # Tkinter selection fallback
        cfg = load_INI()
        if not cfg:
            raise ValueError("could not find config")
        return cfg

    #setup 
    filepaths = sorted( glob.glob(os.path.join(config['DataInfo']['directory'], '*.hdf5')))
    files = [os.path.basename(f) for f in filepaths] 
    fileIDs = [int(item.split('.')[0]) for item in files]

    if len(fileIDs)== 0:
        print(os.path.join(config['DataInfo']['directory'], '*.hdf5'))
        #raise ValueError("Data files cannot be found, check path ends with slash")
    
    if len(config['ProcessingInfo']['starttime'])>0:
        starttime = int(config['ProcessingInfo']['starttime'])
        stoptime = int(config['ProcessingInfo']['stoptime'])
        fileIDs = [i for i in fileIDs if i >= starttime and i <= stoptime]
        if len(fileIDs) == 0:
            raise ValueError("time snippet requested does not exist in File list, check if correct")

    if type(fileIDs[0]) == int: 
            fileIDs = ['{:06d}'.format(fid) for fid in fileIDs]  
    fileIDs_int = np.array(fileIDs, dtype=np.int32)
    assert fileIDs_int[0] == fileIDs_int.min() and fileIDs_int[-1] == fileIDs_int.max()

    n_files = int(config['DataInfo']['n_files'])
    list_fids = [fileIDs[x:x+n_files] for x in range(0, len(fileIDs), n_files)]

    #find channels
    firstfile = os.path.join(config['DataInfo']['directory'] , fileIDs[0] + '.hdf5')
    channels = []
    meta = Calder_utils.load_meta(firstfile)
    chans = meta['header']['channels']

    n_synth = config['ProcessingInfo']['n_synthetic']
    
    match n_synth:
        case '-1':
            channels = None

        case 'auto':
            spacing = int(config['ProcessingInfo']['synthetic_spacing'])
            nstack = int(config['ProcessingInfo']['n_stack'])
            c_start = int(config['ProcessingInfo']['c_start'])
            n_synthetic = np.floor(len(chans)/spacing)
            n_synthetic = int(n_synthetic)
            for i in range(n_synthetic):
                channels.extend([x+(i*spacing)+c_start for x in range(0,nstack)])
            channels[:] = [x for x in channels if x <= len(chans)-1]

        case 'meter':
            dx = int(chans[1]-chans[0])
            spacing = int(np.rint(int(config['ProcessingInfo']['synthetic_spacing'])/dx))
            nstack = int(config['ProcessingInfo']['n_stack'])
            c_start = int(np.rint(int(config['ProcessingInfo']['c_start'])/dx))
            n_synthetic = np.floor(len(chans)/spacing)
            n_synthetic = int(n_synthetic)

            if nstack > spacing:
                nstack = nstack - (nstack - spacing - 1) #makes the nstacks fit nicely within spacing
                print('reducing stack size to fit in spacing.')
            
            for i in range(n_synthetic):
                channels.extend([x+(i*spacing)+c_start for x in range(0,nstack)])
            channels[:] = [x for x in channels if x <= len(chans)-1]

        case _:
            c_start = int(config['ProcessingInfo']['c_start'])
            for i in range(int(config['ProcessingInfo']['n_synthetic'])):
                channels.extend([x+(i*int(config['ProcessingInfo']['synthetic_spacing']))+c_start for x in range(0,int(config['ProcessingInfo']['n_stack']))])
            channels[:] = [x for x in channels if x <= len(chans)-1]
    
    n_workers = int(config['DataInfo']['n_workers'])
    verbose = config['ProcessingInfo'].getboolean('verbose')
    path_data = config['DataInfo']['directory'] 

    #making the output directory
    tnow = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
    outputdir = os.path.join(config['SaveInfo']['directory'],config['SaveInfo']['run_name']+tnow)
    os.makedirs(outputdir)
    
    config['Append'] = {'first':fileIDs[0],
                        'outputdir':outputdir}




    if verbose:
        print(path_data)
        print(outputdir)
        print(channels) 


    # with ProcessPoolExecutor(max_workers= n_workers) as executor:
    #     executor.map(partial(LPS_block, path_data,channels,verbose, config), list_fids)
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(partial(LPS_block, path_data,channels,verbose, config), lf)for lf in list_fids]
        for future in as_completed(futures):
            future.result()
        
if __name__ == '__main__':
    t_ex_start = time.perf_counter()

    ini_list = load_INI()
    if not ini_list:
        raise ValueError("No .ini selected.")

    print(f"\n=== {len(ini_list)} file(s) .ini selected ===")

    for ini_file in ini_list:
        print(f"\n--- Start of {ini_file} ---")
        main(config_path=ini_file)
        print(f"--- End pf {ini_file} ---")

    t_ex_end = time.perf_counter()
    print(f"\n=== Duration: {t_ex_end - t_ex_start:.2f}s ===")




# %%
