# This is a custodial script that allows you to pull profiles out of profs files and print them to the screen.
# This is useful, for instance, if you wish to set up an NTSS run that begins where another ended.
# Actually, this may be a waste of time...

fileToRead = 'profs_1b'

# Load modules
import numpy as np
np.set_printoptions(suppress=True)

# Indices of relevant quantities - use the values listed in the file, they will be converted to Python indices automatically
inds = {
'r':1, # In m
'ne':2, # In 10^20 m^-3
'nHe':5, # In 10^20 m^-3
'Te':6, # In keV
'TD':7, # In keV
'Er':9, # In kV/m
}

# Handy functions
def fixInd(ind):
    return ind - 1

def loadData(filePath):
    skiprows = 1
    return  np.loadtxt(filePath, skiprows=skiprows)

def writeText(ar):
    return np.array2string(ar, separator='    ', precision=7, max_line_width=np.inf).replace('[','').replace(']','')

# Print stuff
data = loadData(fileToRead)

neStr = 'xNe\n\n' + 'r   ' + writeText(data[:, fixInd(inds['r'])]) + '\n' + 'xNe   ' + writeText(data[:, fixInd(inds['ne'])]) + '\n\n' + '**********\n\n'
nHeStr = 'xNHe\n\n' + 'r   ' + writeText(data[:, fixInd(inds['r'])]) + '\n' + 'xNHe   ' + writeText(data[:, fixInd(inds['nHe'])]) + '\n\n' + '**********\n\n'
TeStr = 'xTe\n\n' + 'r   ' + writeText(data[:, fixInd(inds['r'])]) + '\n' + 'xTe   ' + writeText(data[:, fixInd(inds['Te'])]) + '\n\n' + '**********\n\n'
TDStr = 'xTD\n\n' + 'r   ' + writeText(data[:, fixInd(inds['r'])]) + '\n' + 'xTD   ' + writeText(data[:, fixInd(inds['TD'])]) + '\n\n' + '**********\n\n'
ErStr = 'xEr\n\n' + 'r   ' + writeText(data[:, fixInd(inds['r'])]) + '\n' + 'xEr   ' + writeText(data[:, fixInd(inds['Er'])]) + '\n\n' + '**********\n\n'

print(neStr)
print(nHeStr)
print(TeStr)
print(TDStr)
print(ErStr)
