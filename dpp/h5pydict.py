#!/home/camaro/anaconda3/envs/python3/bin/python
# -*- coding: utf-8 -*-

from h5py import File,Group,Dataset

import numpy as np
import warnings

class DictFile(File):
    '''
    Loads the whole hdf5-file into a dictionary tree with numpy-array data.
    Usage:
    Reading:
        with h5py.DictFile('test.hdf5','r') as f:
            data=f.load_dict()
    equivalent to:
        data = h5py.load('test.hdf5')
    Reading only a spesific field/group:
        with h5py.DictFile('test.hdf5','r') as f:
            data=f.load_dict(field='big')
    equivalent to:
        data = h5py.load('test.hdf5',field = 'big')
    Reading all data, except some spesific fields:
        with h5py.DictFile('test.hdf5','r') as f:
            data=f.load_dict(skipFields=['big'])
    equivalent to:
        data = h5py.load('test.hdf5',skipFields = ['big'])
    Saving:
        with h5py.DictFile(filename,'w') as f:
            f.save_dict(data)
    equivalent to:
        h5py.save('test.hdf5',data)
            
    Appending to an existing file:
        with h5py.DictFile(filename,'a') as f:
            f.save_dict(data)

    Only load the keys:
        with h5py.DictFile('test.hdf5','r') as f:
            keys=f.load_keys()
    equivalent to:
        h5py.load_keys('test.hdf5')

    '''
    def __init__(self,filename,mode,verbose=False,compress_ndarray_size=-1,compression_type='gzip',compression_opts=4,**kwargs):
        
        '''
        Paramters:
        ------------        
        filename: str
             name of dict file 
        mode: char
            'r','w','a' etc            
        verbose: bool
        compress_ndarray_size: int
            compress ndarray with gzip when size of array is above compress_ndarray_size
        compression_opts: int
            set compression level
        '''
        self.verbose=verbose
        self.compress_ndarray_size=compress_ndarray_size
        self.compression_opts=compression_opts
        self.compression_type=compression_type
        if verbose:
            print ("Open file: " + filename)
        #import pdb;pdb.set_trace()            
        super(DictFile,self).__init__(filename,mode,**kwargs)


    def save_dict(self,data):
        """Save entire dictionary tree in hdf5 format"""
        #import pdb; pdb.set_trace()
        if self.verbose: print("Save dictinary tree to: " + str(self.filename))
        self.__recursion_save(data,self,0)

    def load_dict(self,field=None,skipFields=[],getData=True):
        """Load entire dictionary tree in hdf5 format.
        Parameters:
        ----------
        field : hdf5 group or str
            Only load a spesific field
            If None, the whole file is loaded.
        skipFields: list 
            list of fields not to load
        getData : list or bool
            is a list keys(at final level) that will be loaded or True to load all.
        """
        if self.verbose: print("Load dictinary from tree: " + str(self.filename))
        data = {}
        self.skipFields = skipFields
        if isinstance(field,str): field = self[field]
        if isinstance(field,Dataset):
            data = field[()]
        else:
            self.__recursion_load(self if field is None else field,data,0,getData)
        return data

    def load_keys(self):
        """Get the keys of the dictionary"""
        if self.verbose: print("Loading tree structure: " + str(self.filename))
        data = {}
        self.__recursion_load(self,data,0,False)
        return data


    def __recursion_save(self,tree,parentgroup, depth = 0):
        if tree == None or len(tree) == 0:
            if self.verbose:  print("\t" * depth, "-")
        else:
            for key, val in tree.items():
#                if key=='phaseDemodPar':
#                    import pdb; pdb.set_trace()
                
                if isinstance(val,dict):
#                    if key in str(parentgroup):
#                        import pdb;pdb.set_trace()
                    if self.verbose:  print("\t" * depth + str(key))
                    newgroup=parentgroup.create_group(str(key))
                    try:
                        self.__recursion_save(val,newgroup, depth+1)
                    except RuntimeError as er:
                        print(val,newgroup, depth+1)
                        
                        raise(er)
                        
                        
                    
                else:
                    compression = None
                    if not isinstance(val,np.ndarray): val = np.array(val)
                    elif (val.size >= self.compress_ndarray_size) and (self.compress_ndarray_size >=0):
                        compression = self.compression_type
                        #print("compress")
                    try:
                        if compression == 'gzip':
                            parentgroup.create_dataset(str(key),data=val,compression=compression,compression_opts=self.compression_opts)
                        else:
                            parentgroup.create_dataset(str(key),data=val,compression=compression)
                    except TypeError as e:
                        #print(key,val,e)
                        # Fix unicode error python3
                        if isinstance(val,np.ndarray):
                            val2=np.array(val,dtype='S')
                            parentgroup.create_dataset(str(key),data=val2)
                        else:
                            print(e)
                            #import pdb;pdb.set_trace()
                                
                        
                if self.verbose: self.__print_info(str(key),val,depth)


    def __recursion_load(self,tree,parentdict, depth = 0, getData=True):
        
        if tree == None or not isinstance(tree,Group):  # siste nivå
           
            if self.verbose:  print("\t" * depth, "-")
        else:
            
            for key, val in tree.items():
                if key in self.skipFields:
                    continue
                if isinstance(key,str) and key.isdigit():
                    key = int(key)
                else:
                    key = str(key)
                if isinstance(val,Group):
                    if self.verbose:  print("\t" * depth + key)
                    parentdict[key]={}
                    self.__recursion_load(val,parentdict[key],depth+1,getData)
                else:
                    #if key=='nSkip':import pdb;pdb.set_trace()
                    if (getData is True) or (isinstance(getData,list) and getData.count(key)>0):
                        if val.dtype.char in 'SO':
                            if val.shape==():
                                if val[()]==b'None':
                                    parentdict[key]=None
                                else:
                                    try:
                                        parentdict[key]=val[()].decode("utf-8")
                                    except:
                                        parentdict[key]=val[()]
                            else:       
                                
                                parentdict[key]=np.array(val).astype(str)
                        elif 'int' in str(val.dtype):
                            if val.shape==():
                                
                                parentdict[key]=int(np.array(val).astype(int))
                            else:
                                parentdict[key]=np.array(val).astype(int)
                            
                        else:
                            
                            parentdict[key]=np.array(val[()])
                            with warnings.catch_warnings():
                                # Cause all warnings to always be triggered.
                                warnings.filterwarnings("error")
                                try:
                                    if parentdict[key]=='None':
                                        parentdict[key]=None
                                except:
                                    pass
                        
                            
                    else:
                        
                        parentdict[key] = val.shape
                    if self.verbose:  self.__print_info(key,val,depth)

    def __print_info(self,key,val,depth):
        
        if (len(val.shape)==0  or (len(val.shape)<2 and len(val)<10)):
            print("\t" * depth + ("%s (%s):" %(key,str(val.dtype))).ljust(30) + str(val[()]))
        else:
            print("\t" * depth + ("%s (%s):" %(key,str(val.dtype))).ljust(30) + str(val.shape))
#        else:
#            print("\t" * depth + ("%s (%s):" %(key,type(val).__name__)).ljust(30)+ str(val))
#        try:
#            print("\t " * depth + ("%s (%s):" %(key,str(val.dtype))).ljust(30) + str(val.shape))
#        except:
#            print("\t " * depth + ("%s (%s):" %(key,type(val).__name__)).ljust(30)\
#                    + ('('+str(len(val))+')' if hasattr(val,'__len__') else str(val)))

def save(filename,data,**kwargs):
    '''save dictionary data to hdf5-file

        Parameters:
        -----------
        filename: str
             name of dict file 
        data: dict
            dictionary with the data to be saved in hdf5-file            
        verbose: bool
        compress_ndarray_size: int
            compress ndarray with gzip when size of array is above compress_ndarray_size
        compression_opts: int
            set compression level
    '''
    #f=DictFile(filename,'w',**kwargs)
    with DictFile(filename,'w',**kwargs) as f: #a better way of handling files that does not end up with a open file on exceptions
        f.save_dict(data)
    #f.close()

def load(filename,field=None,skipFields=[],**kwargs):
    '''read dictionary data from hdf5-file
    
        Parameters:
        -----------
        filename: str
             name of dict file 
        field : str
            Only load a spesific field
            If None, the whole file is loaded.
        skipFields: list 
            list of fields not to load

        verbose: bool
        compress_ndarray_size: int
            compress ndarray with gzip when size of array is above compress_ndarray_size
        compression_opts: int
            set compression level
            
        Returns:
        ----------
        data: dict
            dictionary with the fields in hdf5-file
 
    '''
    with DictFile(filename,'r',**kwargs) as f:
        data =  f.load_dict(field=field,skipFields=skipFields)
    #import pdb;pdb.set_trace()
    
    return data

def load_keys(filename,**kwargs):
    '''read dictionary data from hdf5-file'''
    with DictFile(filename,'r',**kwargs) as f:
        data =  f.load_keys()
    #import pdb;pdb.set_trace()
    
    return data

def rename_field(filename,old_name,new_name):
    ''' Rename a hdf5 field name 
    
    Parameters:
    -----------
    filename str
        filename to load and change field name in
    old name str
        old field name split to subdirionaty using /
        example a/b
    '''
    data=load(filename)
#    
#    
#    if "/" in old_name:
#        okeys=old_name.split('/')
#        nkeys=new_name.split('/')
#        for old_key,new_key in enumerate(len([::-1]:
#            field_val=data.pop(old_name)
#    else:
    field_val=data.pop(old_name)
    data[new_name]=field_val
    save(filename,data)

if __name__=='__main__':

    import numpy as np
    switch=6
    if switch==0:
        #print dict
       from pp_dict import print_tree as pprint;           
       #data=load('/data/fsi/param/sweeps/bw80MHz_sw20us_p100us_s300MHzScale-35/bw80MHz_sw20us_p100us_s300MHzScale-35_005.hdf5')
       data=load('/raid1/fsi/adc/rayxtalk01/00000016.hdf5')
       pprint(data)
       
    elif switch==1:
        #Test rename
        import os
        ffids=range(16,20)
        folder="/raid1/fsi/adc/rayxtalk01/"
        for ffid in ffids:
            print (ffid)
            fname="%08d.hdf5" % ffid
            fname=os.path.join(folder,fname)
             
            rename_field(fname,'chDemuxPar','chDemux')
    elif switch==2:
        # Update field
        import os
        ffids=range(1,26)
        for ffid in ffids:
            try:
                folder="/raid1/fsi/adc/rayxtalk01/"
                print(ffid)
                fname="%08d.hdf5" % ffid
                fname=os.path.join(folder,fname  )          
                data=load(fname)
                data['receiver']={'att':4}
                save(fname,data)
            except IOError:
                pass
    elif switch==3: 
       # Update field in h5py file
       # One need to delete then open as append and then append new datafield
       from pp_dict import print_tree as pprint;            
       import os
       import h5py

       folder='/raid1/fsi/dem/20171013/';ffids=range(10,523)
       #folder='/raid1/fsi/dem/20171014/';ffids=range(334,486)
       #folder='/raid1/fsi/dem/20171011/';ffids=range(10,523):
       #folder='/raid1/fsi/dem/20171012/';ffids=range(364,436)
       #for 
       for ffid in ffids:
           print(ffid)
           fname=os.path.join(folder,'%08d.hdf5' %ffid)
           
           # Read
           f=h5py.File(fname,'r')
           print(f['channelHeader']['sensitivity'][()])
           ds=f['channelHeader']['sensitivity']
           newData=ds[()]*0+10**(-20.6/20)
           f.close()
           if 1:
               #write
               f=h5py.File(fname,'a')
               del f['channelHeader']['sensitivity']
               dset = f.create_dataset('channelHeader/sensitivity', data=newData)
               f.close()
    elif switch==4:        
       import os
       folder='/raid1/fsi/S2/test2/'
       #fname='145323.hdf5' # rlvl
       fname='152719.hdf5' # phase
       
       keys=load_keys(os.path.join(folder,fname))
       if "S2" in keys:
           print('S2',fname)
       elif "phi" in keys:
           print("phi",fname)
       elif "DECxI" in keys:
           print("DECxI",fname)
           
           
       if 0:
           data=load(os.path.join(folder,fname),verbose=True)
           
           folder='/raid1/fsi/dem/26kmFBG_0k05_02/'
           fname='%08d.hdf5'%12 # phase
           data2=load(os.path.join(folder,fname),verbose=True)
           

    elif switch==5:
         import os,sys
         thispath=os.path.dirname(os.path.realpath(__file__));
         sys.path.append(thispath+'/../FSI/')
         ffidpath='/raid1/fsi/exps/20190426/dem/'
         ffid=75120 # Forward sweep
         ffid=85045 # Reverse sweep
         fname='%06d.hdf5'%ffid # phase
         ffidpath=os.path.join(ffidpath,fname)
         from  fsi_common2 import load_ffid
         from logger import logit;
         from fsi_demod2 import Fsi;

         cfsi=Fsi(ffidpath)
         cfsi.load_ffid()

    elif switch==6:
        ffidfile = '/opt/data/fsi/exps/ALMA2_self_noise/20200401/proc/095810.hdf5'
        d = load(ffidfile)#,field='monitoring')
        print(d)
