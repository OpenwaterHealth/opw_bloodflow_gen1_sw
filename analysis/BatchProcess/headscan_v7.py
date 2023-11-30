# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 11:55:05 2020 

@author: soren
"""

from re import sub
from HeadscanData import * #Import all classes
import numpy as np
import os
import batchfile
import proplot as pplt

### file info ###
batch = True
batchname = 'HHC_1B_Schrodinger_ICU'
#'HHC_1B_Schrodinger_Healthy' 'HHC_1B_Schrodinger_ICU'
#'FrontSt_Schrodinger_20210511and12' 'PulseModulationExperiments' 'WaterfallExperiments'
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
#########################################
###### do not modify below ##############
#########################################

calnames, scannames, datapath = eval('batchfile.' + batchname + '()') #batchfile.getscannames()
nsubjects = len(scannames)
datapath = os.path.expanduser(datapath)
savepath = datapath + '/Results'

numScansTot = 0
for scanName in scannames:
    numScansTot = numScansTot+len(scanName)

dataList = [{}] * numScansTot
subjectID = np.zeros(numScansTot)
ind = 0
### loop of subjects and scans
for ii in range(nsubjects):
    if batch:
        nscans = len(scannames[ii])
        #plt.close('all')
    for jj in range(nscans): 
        
        print('subject ' + str(ii+1) + ',  scan ' + str(jj+1) + ' ' + scannames[ii][jj])
    
        ### data and algo parameters ###
        if batch:
            calname = datapath + '/' + calnames[ii]
            scanname = datapath + '/' + scannames[ii][jj]
        
        #Testing new object oriented code
        scanParams  = HeadScanParams(scanname, calname)
        subjectScan = HeadScanData(scanParams)
        subjectScan.LoadScanDataAndPreprocess(cropTimeStart=1,cropTimeEnd=-1)
        subjectScan.RunWaveletFilteringOnSignal()
        if 1:#ii:
            subjectScan.GetPulseFeatures(lowNoise=True)
        else:
            subjectScan.GetPulseFeatures()
        '''for ind,scan in enumerate(subjectScan.scanData):
            fig, ax = pplt.subplots(nrows=1, ncols=1)
            ax.format(suptitle='Patient #'+str(ii+4)+' scan '+str(ind), ylabel='contrast')
            ax[0].plot(scan.contrast)
            pplt.show()'''
        dataList[ind] = subjectScan
        subjectID[ind] = ii+1
        ind = ind+1
asd = 0