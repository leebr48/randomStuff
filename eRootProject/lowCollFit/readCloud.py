# This script is meant to read the *_opt* files generated by the STELLOPT GADE optimizer.
# You can run the optimizer for 1 generation to generate a 'cloud' of random configurations,
# then use this script to pull data for those configurations.

# User options #FIXME this should probably be made into argparse stuff

loadPickle = 'allData.pkl' # None to load new data, path of pickle file to load old data
directory = '/cobra/u/lebra/data/w7x/reactor/twentySecondObj12b' # Should be able to handle absolute and relative paths
ext = 'W7X_REACTOR_woptim_forSfincs' #file extension, input.ext
#FIXME I guess you could put filters for the surfaces you look at and such in here?... How best to do it, I'm not sure.

####################################################################################################

# Import necessary modules

import warnings
import os
import numpy as np
import pickle as pkl
from glob import glob
from scipy.io import netcdf_file

# Define handy functions

warnings.filterwarnings("error") # Necessary so the UserWarning saying 'this file is empty' is caught in the function below
def loadData(file):
    
    try:
        data = np.loadtxt(file)
    except UserWarning:
        data = None

    return data

def get_eps_eff(outFile):
    
    rawData = loadData(outFile)

    if rawData is None:
        eps_eff = None
    else:
        eps_eff = rawData[:,1] ** (2/3)
    
    return eps_eff

def get_DKES_data(optFile):

    rawData = loadData(optFile)

    if rawData is None:
        DKESdata = None
    else:
        DKESdata = rawData

    return DKESdata

def sort_DKES_data(optFile):

    DKESdata = get_DKES_data(optFile)

    if DKESdata is None:
        sortedDKESdata = None
    else: # This uses notation similar to the dkesout.ext_opt#_s* and stellopt.ext files
        sortedDKESdata = {
                          'L11p':    DKESdata[0,0],
                          'L11m':    DKESdata[1,0],
                          'L33p':    DKESdata[0,1],
                          'L33m':    DKESdata[1,1],
                          'L31p':    DKESdata[0,2],
                          'L31m':    DKESdata[1,2],
                          'SCALE11': DKESdata[2,0],
                          'SCALE33': DKESdata[2,1],
                          'SCALE31': DKESdata[2,2]
                         }

    return sortedDKESdata

def filter_DKES_data(optFile):

    sortedDKESdata = sort_DKES_data(optFile)

    if sortedDKESdata is None:
        filteredDKESdata = None
    else:
        filteredDKESdata = np.average((sortedDKESdata['L11p'], sortedDKESdata['L11m']))

    return filteredDKESdata

# Get data

if loadPickle is None:

    # Identify all configurations in the cloud

    basewoutFileName = os.path.join(directory, 'wout_' + ext)
    matchingwoutFiles = sorted(glob(basewoutFileName + '_opt*.nc'))

    # Loop through files, extract desired data, and pair with the configuration

    allData = {}

    for woutFile in matchingwoutFiles:
        
        # Instantiate a dictionary that will hold all the information forthis configuration

        configData = {}

        # Isolate the opt number to keep configurations straight

        optNum = int(woutFile.split('_opt')[-1].strip().split('.')[0])

        if optNum == 0: # There is an 'extra' wout file, so we skip it
            continue
       
        # Get configuration data

        wout = netcdf_file(woutFile, mode='r', mmap=False)

        vmecData = {
                    'aspect':    wout.variables['aspect'][()],
                    'betatotal': wout.variables['betatotal'][()],
                    'b0':        wout.variables['b0'][()],
                    'volavgB':   wout.variables['volavgB'][()],
                    'Aminor_p':  wout.variables['Aminor_p'][()],
                    'Rmajor_p':  wout.variables['Rmajor_p'][()],
                    'volume_p':  wout.variables['volume_p'][()],
                    'iota':      wout.variables['iotaf'][()],
                    'rmnc':      wout.variables['rmnc'][()],
                    'zmns':      wout.variables['zmns'][()]
                   }

        configData['vmecData'] = vmecData

        # Get eps_eff data
        
        epseffFileName = os.path.join(directory, 'neo_out.' + ext + '_opt' + str(optNum))
        eps_eff = get_eps_eff(epseffFileName)
        configData['eps_eff'] = eps_eff

        # Get DKES data
        
        baseDKESFileName = os.path.join(directory, 'opt_dkes.' + ext + '_opt' + str(optNum))
        matchingDKESFileNames = sorted(glob(baseDKESFileName + '_s*'))

        DKESdata = np.array([])
        
        for DKESFile in matchingDKESFileNames:
            DKESdata = np.append(DKESdata, filter_DKES_data(DKESFile))
        
        configData['L11'] = DKESdata

        # Append configuration data to the allData dictionary

        allData[optNum] = configData

    with open('allData.pkl', 'wb') as f:
        pkl.dump(allData, f)

else:
    
    with open(loadPickle, 'rb') as f:
        allData = pkl.load(f)

# FIXME you're going to need to write your own stuff to normalize to D_11^* and nu^*, then use Hakan's fit function.
