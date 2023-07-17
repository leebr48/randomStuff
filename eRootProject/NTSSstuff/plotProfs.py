# Settings
profsFilePath = 'profs_1b'
rhoMin = 0
rhoMax = 0.85 # Especially useful if the edge is broken...
useRho = True
savePlots = True
showPlots = True
fileExt = 'png'

# Indices of relevant quantities - use the values in the 'profs' file, they will be converted to Python indices automatically
inds = {
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
'L11e':16, # In m^2/s
'L11D':20, # In m^2/s
'L11T':24, # In m^2/s
'Ibs':63, # In A
'vaciota':69, # Uses susceptance matrices to exclude bootstrap current from the calculation, like vaciota in STELLOPT
'iota':70 # Uses susceptance matrices to include bootstrap current in the calculation
}

################################################################
# Import necessary modules
import numpy as np
import matplotlib.pyplot as plt

# Handy functions
def fixInd(ind):
    return ind - 1

def loadVec(ind):
    return filteredData[:, fixInd(ind)]

def makePlot(xdata, ydata, ylabel, figName, leg=None, fileExt=fileExt, yticks=None, ymin=None):
    plt.figure()
    plt.plot(xdata, ydata)
    if ymin is not None:
        plt.ylim(ymin=ymin)
    if yticks is not None:
        plt.yticks(np.arange(np.floor(np.min(ydata)), np.ceil(np.max(ydata))+yticks, yticks))
    plt.xlabel(xlab)
    plt.ylabel(ylabel)
    if leg is not None:
        plt.legend(leg, loc='best')
    if savePlots:
        plt.savefig(figName+'.'+fileExt, bbox_inches='tight', dpi=400)

def multiPlot(xdata, ydataList):
    y = np.column_stack(ydataList)
    x = np.tile(xdata, (y.shape[1],1)).T
    return x, y

# Load data
data = np.loadtxt(profsFilePath, skiprows=1)
rInd = fixInd(inds['r'])
a = np.max(data[:, rInd])
filteredData = data[(data[:, rInd]/a >= rhoMin) & (data[:, rInd]/a <= rhoMax)]

# Grab relevant variables
vecs = {}
for variable, index in inds.items():
    vecs[variable] = loadVec(index)

# Plot things
if useRho:
    xData = vecs['r'] / a
    xlab = r'$\rho$'
else:
    xData = vecs['r']
    xlab = r'$r$ (m)'

makePlot(*multiPlot(xData, [vecs['ne'], vecs['nD'], vecs['nT'], vecs['nHe']]), r'Density ($10^{20}~\mathrm{m^{-3}}$)', 'n', leg=['e', 'D', 'T', 'He'], ymin=0)
makePlot(*multiPlot(xData, [vecs['Te'], vecs['TD'], vecs['TT'], vecs['TT']]), r'Temperature (keV)', 'T', leg=['e', 'D', 'T', 'He'], ymin=0)
makePlot(xData, vecs['Er'], r'Radial Electric Field (kV/m)', 'Er')
makePlot(xData, vecs['p'], r'Pressure (Pa)', 'p', ymin=0)
makePlot(*multiPlot(xData, [vecs['L11e'], vecs['L11D'], vecs['L11T']]), r'$L_{11}$ ($\mathrm{m^2/s}$)', 'L11s', leg=['e', 'D', 'T'], ymin=0)
makePlot(xData, vecs['L11e']/(0.5*(vecs['L11D']+vecs['L11T'])), r'$ 2 L_{11}^{e} / \left(L_{11}^{D}+L_{11}^{T}\right) $', 'L11rat', yticks=0.5)
makePlot(xData, vecs['Ibs'] / 1000, r'Bootstrap Current (kA)', 'Ibs')
makePlot(xData, vecs['vaciota'], r'Vacuum Rotational Transform', 'vaciota')
makePlot(xData, vecs['iota'], r'Rotational Transform', 'iota')

if showPlots:
    plt.show()
