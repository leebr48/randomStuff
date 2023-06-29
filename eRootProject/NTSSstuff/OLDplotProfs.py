# FIXME probably get the volume of your run spot on 1900 m3 and remake these plots... if that configuration isn't totally broken
# Settings
profsFilePath = 'profs'
rmin = 0 # In m.
rmax = 1.3 # In m. Especially useful if the edge is broken...
useRho = True
savePlots = True
showPlots = True

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
'iota':71
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

def makePlot(xdata, ydata, ylabel, figName, leg=None, fileExt='pdf'):
    plt.figure()
    plt.plot(xdata, ydata)
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
filteredData = data[(data[:, rInd] >= rmin) & (data[:, rInd] <= rmax)]

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

makePlot(*multiPlot(xData, [vecs['ne'], vecs['nD'], vecs['nT'], vecs['nHe']]), r'Density ($10^{20}~\mathrm{m^{-3}}$)', 'n', leg=['e', 'D', 'T', 'He'])
makePlot(*multiPlot(xData, [vecs['Te'], vecs['TD'], vecs['TT'], vecs['TT']]), r'Temperature (keV)', 'T', leg=['e', 'D', 'T', 'He'])
makePlot(xData, vecs['Er'], r'Radial Electric Field (kV/m)', 'Er')
makePlot(xData, vecs['p'], r'Pressure (Pa)', 'p')
makePlot(*multiPlot(xData, [vecs['L11e'], vecs['L11D'], vecs['L11T']]), r'Thermal Transport Coefficients ($\mathrm{m^2/s}$)', 'L11s', leg=['e', 'D', 'T'])
makePlot(xData, vecs['L11e']/(0.5*(vecs['L11D']+vecs['L11T'])), r'$ 2 L_{11}^{e} / \left(L_{11}^{D}+L_{11}^{T}\right) $', 'L11rat')
makePlot(xData, vecs['Ibs'] / 1000, r'Bootstrap Current (kA)', 'Ibs')
makePlot(xData, vecs['iota'], r'Rotational Transform', 'iota')

if showPlots:
    plt.show()
