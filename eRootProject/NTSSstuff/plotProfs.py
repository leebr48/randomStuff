# Settings
profsFilePath = 'profs_1c' # Path to the main file. Auxiliary files (such as *_Dij) will also be loaded by the program. # Should be 'profs_1c', 'profs_2b', 'profs_3b', or 'profs_w7xhm'
figNamePrefix = 'lee1'
rhoMin = 0
rhoMax = 1 # Especially useful if the edge is broken... # Should be 0.85 or 1
nMax = 2.5
TMax = 25.1
axisFontSize = 24
legendFontSize = 14
xSizeInches = 8
ySizeInches = 6
useRho = True
showTempScreenThresh = False
savePlots = True
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
from matplotlib.ticker import FuncFormatter

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

def my_formatter(x, pos):
    if x == 0:
        return str(int(0))
    elif x == 1:
        return str(int(1))
    else:
        return '{0:.1f}'.format(x)

def makePlot(xdata, ydata, ylabel, figName, leg=None, linestyles=None, fileExt=fileExt, yticks=None, ymin=None, ymax=None):
    formatter = FuncFormatter(my_formatter)
    plt.subplots(figsize=(xSizeInches, ySizeInches))
    if linestyles is None:
        plt.plot(xdata, ydata)
    else:
        for X, Y, L in zip(xdata.T, ydata.T, linestyles):
            plt.plot(X, Y, linestyle=L)
    plt.xlim(xmin=np.min(xdata), xmax=np.max(xdata))
    plt.gca().xaxis.set_major_formatter(formatter)
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

# Plot things
if useRho:
    xData = vecs['r'] / a
    xlab = r'$\rho$'
else:
    xData = vecs['r']
    xlab = r'$r$ (m)'

stdLeg = ['e', 'D', 'T', 'He']
stdStyle = ['solid', 'solid', 'dashed', 'dotted']

makePlot(*multiPlot(xData, [vecs['ne'], vecs['nD'], vecs['nT'], vecs['nHe']]), r'Density ($10^{20}~\mathrm{m^{-3}}$)', figNamePrefix + '_n', leg=stdLeg, linestyles=stdStyle, ymin=0, ymax=nMax)
makePlot(*multiPlot(xData, [vecs['Te'], vecs['TD'], vecs['TT'], vecs['TT']]), r'Temperature (keV)', figNamePrefix + '_T', leg=stdLeg, linestyles=stdStyle, ymin=0, ymax=TMax)
makePlot(xData, vecs['Er'], r'Radial Electric Field (kV/m)', figNamePrefix + '_Er')
makePlot(xData, vecs['p'], r'Pressure (Pa)', figNamePrefix + '_p', ymin=0)
makePlot(*multiPlot(xData, [vecs['L11e'], vecs['L11D'], vecs['L11T']]), r'$L_{11}$ ($\mathrm{m^2/s}$)', figNamePrefix + '_L11s', leg=['e', 'D', 'T'], linestyles=stdStyle[:-1], ymin=0)
if showTempScreenThresh:
    makePlot(*multiPlot(xData, [vecs['L11e']/(0.5*(vecs['L11D']+vecs['L11T'])), tempScreenThresh]), r'$ 2 L_{11}^{e} / \left(L_{11}^{D}+L_{11}^{T}\right) $', figNamePrefix + '_L11rat', leg=['Actual', 'He Tmp. Scrn. Thresh.'], yticks=0.5, ymin=0)
else:
    makePlot(xData, vecs['L11e']/(0.5*(vecs['L11D']+vecs['L11T'])), r'$ 2 L_{11}^{e} / \left(L_{11}^{D}+L_{11}^{T}\right) $', figNamePrefix + '_L11rat', yticks=0.5)
makePlot(xData, vecs['Ibs'] / 1000, r'Bootstrap Current (kA)', figNamePrefix + '_Ibs')
makePlot(xData, vecs['vaciota'], r'Vacuum Rotational Transform', figNamePrefix + '_vaciota')
makePlot(xData, vecs['iota'], r'Rotational Transform', figNamePrefix + '_iota')

if showPlots:
    plt.show()
