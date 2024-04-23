# Settings
profsFilePaths = ['profs_1b', 'profs_2_eroot', 'profs_Lee_3'] # Path to the main file. Auxiliary files (such as *_Dij) will also be loaded by the program.
figNamePrefix = 'Ers'
rhoMin = 0
rhoMax = 1 # Especially useful if the edge is broken... # Should be 0.85 or 1
axisFontSize = 24
legendFontSize = 14
xSizeInches = 8
ySizeInches = 6
useRho = True
showTempScreenThresh = False
savePlots = False # FIXME
showPlots = True
fileExt = 'pdf'
dpi = 600

# Indices of relevant quantities - use the values listed in the file, they will be converted to Python indices automatically
mainInds = {
'r':1, # In m
'ne':2, # In 10^20 m^-3
'nD':3, # In 10^20 m^-3
'nT':4, # In 10^20 m^-3
'nHe':5, # In 10^20 m^-3
'Te':6, # In keV
'TD':7, # In keV
'TT':8, # In keV
'Er':9, # In kV/m
'p':11, # In Pa
'Ibs':63, # In A
'vaciota':69, # Uses susceptance matrices to exclude bootstrap current from the calculation, like vaciota in STELLOPT
'iota':70 # Uses susceptance matrices to include bootstrap current in the calculation
}

LInds = {
'r':1, # In m
'L11e':2, # In m^2/s
'L12e':3, # In m^2/s
'L11D':8, # In m^2/s
'L12D':9, # In m^2/s
'L11T':14, # In m^2/s
'L12T':15, # In m^2/s
'L11He':20, # In m^2/s
'L12He':21 # In m^2/s
}

################################################################
# Import necessary modules
import numpy as np
import matplotlib.pyplot as plt

plt.rc('font', size=axisFontSize)
plt.rc('legend', fontsize=legendFontSize)

# Handy functions
def fixInd(ind):
    return ind - 1

def loadData(filePath, fileType='main'):
    global a
    if fileType == 'main':
        skiprows = 1
        inds = mainInds
    elif fileType == 'L':
        skiprows = 2
        inds = LInds
    else:
        raise IOError('Unkown fileType.')
    data = np.loadtxt(filePath, skiprows=skiprows)
    rInd = fixInd(inds['r'])
    a = np.max(data[:, rInd])
    return data[(data[:, rInd]/a >= rhoMin) & (data[:, rInd]/a <= rhoMax)]

def loadVec(ind, fileType='main'):
    if fileType == 'main':
        return mainFilteredData[:, fixInd(ind)]
    elif fileType == 'L':
        return LFilteredData[:, fixInd(ind)]
    else:
        raise IOError('Unkown fileType.')

def makePlot(xdata, ydata, ylabel, figName, leg=None, linestyles=None, fileExt=fileExt, yticks=None, ymin=None, ymax=None):
    plt.subplots(figsize=(xSizeInches, ySizeInches))
    if linestyles is None:
        plt.plot(xdata, ydata)
    else:
        for X, Y, L in zip(xdata.T, ydata.T, linestyles):
            plt.plot(X, Y, linestyle=L)
    if ymin is not None:
        plt.ylim(ymin=ymin)
    if ymax is not None:
        plt.ylim(ymax=ymax)
    if yticks is not None:
        plt.yticks(np.arange(np.floor(np.min(ydata)), np.ceil(np.max(ydata))+yticks, yticks))
    plt.xlabel(xlab)
    plt.ylabel(ylabel)
    if leg is not None:
        plt.legend(leg, loc='best')
    if savePlots:
        plt.savefig(figName+'.'+fileExt, bbox_inches='tight', dpi=dpi)

def multiPlot(xdata, ydataList):
    y = np.column_stack(ydataList)
    x = np.tile(xdata, (y.shape[1],1)).T
    return x, y

def calcDelta12(species):
    L11 = 'L11' + species
    L12 = 'L12' + species
    return vecs[L12] / vecs[L11] #NOTE: I think NTSS uses the old Maassberg definition (without the 3/2)... if you include the 3/2 you get negative numbers in nonsensical places

radVecs = []
ErVecs = []
L11RatVecs = []
IbsVecs = []
iotaVecs = []

for profsFilePath in profsFilePaths:

    # Load data
    mainFilteredData = loadData(profsFilePath, fileType='main')
    LFilteredData = loadData(profsFilePath + '_Dij', fileType='L')

    # Grab relevant variables
    vecs = {}
    for variable, index in mainInds.items():
        vecs[variable] = loadVec(index, fileType='main')
    for variable, index in LInds.items():
        vecs[variable] = loadVec(index, fileType='L')

    # Calculate temperature screening threshold using Beidler's 2022 Simons talk
    d12e = calcDelta12('e')
    d12D = calcDelta12('D')
    d12T = calcDelta12('T')
    d12He = calcDelta12('He')
    gamma = vecs['nT'] / vecs['nD']
    d12i = (vecs['L11D'] * d12D + gamma * vecs['L11T'] * d12T) / (vecs['L11D'] + gamma * vecs['L11T'])
    ni = vecs['nD'] + vecs['nT']
    ZHe = 2
    tempScreenThresh = (ZHe * d12i - d12He) / (ZHe * d12e + d12He) / (1 + ZHe * vecs['nHe'] / ni)
    
    # Append info to our arrays that will be plotted
    if useRho:
        radVecs.append(vecs['r'] / a)
        xlab = r'$\rho$'
    else:
        radVec.append(vecs['r'])
        xlab = r'$r$ (m)'
    ErVecs.append(vecs['Er'])
    L11RatVecs.append(vecs['L11e']/(0.5*(vecs['L11D']+vecs['L11T'])))
    IbsVecs.append(vecs['Ibs'])
    iotaVecs.append(vecs['iota'])

# Plot things
makePlot(np.asarray(radVecs).T, np.asarray(ErVecs).T, 'test', 'test')
plt.show()
quit()
makePlot(*multiPlot(xData, [vecs['ne'], vecs['nD'], vecs['nT'], vecs['nHe']]), r'Density ($10^{20}~\mathrm{m^{-3}}$)', figNamePrefix + '_n', leg=stdLeg, linestyles=stdStyle, ymin=0, ymax=nMax)

if showPlots:
    plt.show()
