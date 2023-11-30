# %%
from HeadscanData import * #Import all classes
import matplotlib.pyplot as plt
#import proplot as pplt
import numpy as np
import pandas as pd
import batchfile

np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
plt.close('all')

##################################
### load data set #1 (Healthy) ###
##################################

### file info ###
scanorder = 'hhc30'#'cuff30'#'hhc30' #'zigzag' #'circle'
batch = True
saveresult = True
offset = 0  #0=none, 1=ax+b, 2=ax+b+c/x
fitalltimes = False
batchname = 'HHC_1B_Schrodinger_Healthy'

calnames, scannames, datapath = eval('batchfile.' + batchname + '()')
nsubjects = len(scannames)
datapath = os.path.expanduser(datapath)
savepath = datapath + '/Results'

numScansTot = 0
for scanName in scannames:
    numScansTot = numScansTot+len(scanName)

dataList = [[]] * numScansTot  #will be list of objects which contain all scan info
subjectID = np.zeros(numScansTot)
ind = 0
### loop of subjects and scans
for ii in range(nsubjects):
    if batch:
        nscans = len(scannames[ii])
    for jj in range(nscans): 
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
        dataList[ind] = subjectScan
        subjectID[ind] = ii+1
        ind = ind+1
        print('Healthy Scans completed = ' + str(ind))

##############################
### load data set #2 (ICU) ###
##############################

### file info ###
scanorder = 'hhc30'#'cuff30'#'hhc30' #'zigzag' #'circle'
batch = True
saveresult = True
offset = 0  #0=none, 1=ax+b, 2=ax+b+c/x
fitalltimes = False
batchname = 'HHC_1B_Schrodinger_ICU'

calnames, scannames, datapath = eval('batchfile.' + batchname + '()')
nsubjects = len(scannames)
datapath = os.path.expanduser(datapath)
savepath = datapath + '/Results'

numScansTot = 0
for scanName in scannames:
    numScansTot = numScansTot+len(scanName)

dataListICU = [[]] * numScansTot #list of objects for ICU data
subjectIDICU = np.zeros(numScansTot)
ind = 0
### loop of subjects and scans
for ii in range(nsubjects):
    if batch:
        nscans = len(scannames[ii])
    for jj in range(nscans): 
        ### data and algo parameters ###
        if batch:
            calname = datapath + '/' + calnames[ii]
            scanname = datapath + '/' + scannames[ii][jj]
        
        #Testing new object oriented code
        scanParams  = HeadScanParams(scanname, calname)
        subjectScan = HeadScanData(scanParams)
        subjectScan.LoadScanDataAndPreprocess(cropTimeStart=1,cropTimeEnd=99)
        subjectScan.RunWaveletFilteringOnSignal()
        subjectScan.GetPulseFeatures(lowNoise=True)
        dataListICU[ind] = subjectScan
        subjectIDICU[ind] = ii+1
        ind = ind+1
        print('ICU Scans completed = ' + str(ind))

# %%
        
def medianParameter(parameter, goodScans, npositions, nrepeats, camera):
    npatients = len(goodScans)
    medianValues = np.zeros((npatients, npositions))        
    for patient in range(npatients):  #for each patient
        print('++++ patient ' + str(patient))
        values = np.nan*np.zeros((3,npositions))
        for position in range(npositions):    #for each position
            for repeat in range(nrepeats):    #for each scan
                print(position)
                if (position+1) in goodScans[patient][repeat]:   
                    #values[repeat, position]= dataList[3*patient+repeat].amplitude[position, camera]
                    mystring = 'dataList[3*patient+repeat].' + parameter + '[position, camera]'
                    if parameter=='intensityMeansAtPos' or parameter=='contrastAvgsAtPos': 
                        mystring = 'dataList[3*patient+repeat].scanData[position].' + parameter + '[camera]'
                    values[repeat, position]= eval(mystring)
        medianValues[patient, :] = np.nanmedian(values, axis=0)
    return medianValues


Healthy = dict()
ICU = dict()
goodScansHealthy = batchfile.HHC_1B_Schrodinger_Healthy_Validation()
goodScansICU   = batchfile.HHC_1B_Schrodinger_ICU_Validation()
npositions = 30
nrepeats = 3
cameras = [0, 3]
parameters = ['intensityMeansAtPos', 'contrastAvgsAtPos', 'amplitude', 'areaUnderCurve', 'kurtosis', 'modulationDepth', 'pulseCanopy', 'pulseOnset', 'pulseOnsetProp', 'pulseRate', 'secondMoment', 'skewness']

for parameter in parameters:
    for camera in cameras: 
        Healthy[parameter+str(camera)] = medianParameter(parameter, goodScansHealthy, npositions, nrepeats, camera)
        Healthy[parameter+str(camera)+'LR'] = Healthy[parameter+str(camera)][:, :15] - Healthy[parameter+str(camera)][:, 15:]
        ICU[parameter+str(camera)] = medianParameter(parameter, goodScansICU, npositions, nrepeats, camera)
        ICU[parameter+str(camera)+'LR'] = ICU[parameter+str(camera)][:, :15] - ICU[parameter+str(camera)][:, 15:]

#%%
def numberplots(Healthy, ICU, ids_healthy, ids_icu, cameras, npairs):
        
    for parameter in parameters:
        for camera in cameras:
            title = parameter+str(camera)
            if npairs==15:
                title += 'LR'
            hdata = Healthy[title]
            idata = ICU[title]
            plt.figure()
            plt.title(title)
            linelist = np.arange(-0.5, npairs, 3)
            miny = np.minimum(np.nanmin(hdata), np.nanmin(idata))
            maxy = np.maximum(np.nanmax(hdata), np.nanmax(idata))
            plt.plot([14.4, 14.4], [miny, maxy], ':', color=[0.5, 0.5, 0.5])
            plt.plot([14.6, 14.6], [miny, maxy], ':', color=[0.5, 0.5, 0.5])
            for i in linelist:
                plt.plot([i, i], [miny, maxy], ':', color=[0.5, 0.5, 0.5])
            for pair in range(npairs):
                ys = hdata[:, pair]
                ys = ys[np.isfinite(ys)]
                xs = pair*np.ones((len(ys),))
                ind=0
                for x, y in zip(xs, ys):
                    plt.text(x, y, str(ids_healthy[ind]), color="green", fontsize=12, va='center', ha='center')        
                    ind += 1    
                ys = idata[:, pair]
                ys = ys[np.isfinite(ys)]
                xs = pair*np.ones((len(ys),))
                ind=0
                for x, y in zip(xs, ys):
                    if ids_icu[ind]==8: #ICH
                        color = 'yellow'
                    elif ids_icu[ind]<4:
                        color = 'blue'
                    else:
                        color = 'red'
                    plt.text(x, y, str(ids_icu[ind]), color=color, fontsize=12, va='center', ha='center')        
                    ind += 1        
        
        
ids_healthy = [35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 58, 59, 60, 61]
ids_icu = [1, 4, 5, 6, 7, 8]
cameras = [0, 3]
npairs = int(npositions/2) 
numberplots(Healthy, ICU, ids_healthy, ids_icu, cameras, npairs)
#numberplots(Healthy, ICU, ids_healthy, ids_icu, cameras, 30)       

"""
# %%
############################################
### put healthy data in Pandas dataframe ###
############################################        
#plotPath = '/Users/soren/Desktop/plots/'
cameras   = 4
totalGoodPos = 0
bucket = [0, 0, 2, 4, 4, 1, 1, 3, 5, 5] #Map to waterfall experiment locations - we divide by 3 to remain agnostic to measurement position(1 or 2 or 3)
goodScansHealthy = batchfile.HHC_1B_Schrodinger_Healthy_Validation()
for patient in goodScansHealthy:
    for scan in patient:
        totalGoodPos += len(scan)

#Subject, Group, Camera, Location, Oritentation, RelOnset, Canopy, Onset
#Kurtosis, Skewness, ContrastMean, IntensityMean, Amplitude, AUC, ModDepth, SecondMoment
#Group 1-handheld, 2-waterfall@14Hz, 3-waterfall@28Hz, 10-HHC ICU
allFeatureVals = np.zeros( ((totalGoodPos+2)*4*30, 16) )
i = 0
for sid, patient in enumerate(goodScansHealthy):
    for scan in patient:
        for goodPos in scan:
            print(i)
            allFeatureVals[i:i+4,0] = sid #Subject
            allFeatureVals[i:i+4,1] = 1 #Group
            allFeatureVals[i:i+4,2] = np.asarray([0,1,2,3]) #Camera
            allFeatureVals[i:i+4,3] = bucket[int((goodPos-1)/3)] #Location
            allFeatureVals[i:i+4,4] = goodPos #Position number 1-30 hhc30 order
            allFeatureVals[i:i+4,5] = dataList[sid].pulseOnsetProp[goodPos-1,:] #Relative Onset
            allFeatureVals[i:i+4,6] = dataList[sid].pulseCanopy[goodPos-1,:] #Canopy
            allFeatureVals[i:i+4,7] = dataList[sid].pulseOnset[goodPos-1,:] #Onset
            allFeatureVals[i:i+4,8] = dataList[sid].kurtosis[goodPos-1,:] #Kurtosis
            allFeatureVals[i:i+4,9] = dataList[sid].skewness[goodPos-1,:] #Skewness
            allFeatureVals[i:i+4,10] = dataList[sid].scanData[goodPos-1].contrastAvgsAtPos #ContrastMean
            allFeatureVals[i:i+4,11] = dataList[sid].scanData[goodPos-1].intensityMeansAtPos #IntensityMean
            allFeatureVals[i:i+4,12] = dataList[sid].amplitude[goodPos-1,:] #Amplitude
            allFeatureVals[i:i+4,13] = dataList[sid].areaUnderCurve[goodPos-1,:] #AUC
            allFeatureVals[i:i+4,14] = dataList[sid].modulationDepth[goodPos-1,:] #ModDepth
            allFeatureVals[i:i+4,15] = dataList[sid].secondMoment[goodPos-1,:] #SecondMoment
            i += 4
WaterfallFeaturesPD = pd.DataFrame(allFeatureVals, columns=[
    'Subject','Group','Camera','Location','Orientation',
    'RelativeOnset','Canopy','Onset','Kurtosis','Skew',
    'ContrastMean', 'IntensityMean', 'Amplitude', 'AUC', 'ModDepth', 'SecMom'])

# %%
############################################
###   put ICU data in Pandas dataframe   ###
############################################    
goodScansICU   = batchfile.HHC_1B_Schrodinger_ICU_Validation()
#scanLocNOrient = ['Forehead Vertical','Forehead Horizontal','Fissure','Temple Vertical','Temple Horizontal']
totalGoodPos = 0
for patient in goodScansICU:
    for scan in patient:
        totalGoodPos += len(scan)
allFeatureVals = np.zeros( (totalGoodPos*4, 16) ) #Patient, Camera, Location, Oritentation, RelOnset, Canopy, Onset
i = 0
for pid, patient in enumerate(goodScansICU):
    for scan in patient:
        for goodPos in scan:
            allFeatureVals[i:i+4,0] = pid+4 #Patient
            allFeatureVals[i:i+4,1] = 10 #Group
            allFeatureVals[i:i+4,2] = np.asarray([0,1,2,3]) #Camera
            allFeatureVals[i:i+4,3] = bucket[int((goodPos-1)/3)] #Location
            allFeatureVals[i:i+4,4] = goodPos #Position number 1-30 hhc30 order
            allFeatureVals[i:i+4,5] = dataListICU[pid].pulseOnsetProp[goodPos-1,:] #Relative Onset
            allFeatureVals[i:i+4,6] = dataListICU[pid].pulseCanopy[goodPos-1,:] #Canopy
            allFeatureVals[i:i+4,7] = dataListICU[pid].pulseOnset[goodPos-1,:] #Onset
            allFeatureVals[i:i+4,8] = dataListICU[pid].kurtosis[goodPos-1,:] #Kurtosis
            allFeatureVals[i:i+4,9] = dataListICU[pid].skewness[goodPos-1,:] #Skewness
            allFeatureVals[i:i+4,10] = dataListICU[pid].scanData[goodPos-1].contrastAvgsAtPos #ContrastMean
            allFeatureVals[i:i+4,11] = dataListICU[pid].scanData[goodPos-1].intensityMeansAtPos #IntensityMean
            allFeatureVals[i:i+4,12] = dataListICU[pid].amplitude[goodPos-1,:] #Amplitude
            allFeatureVals[i:i+4,13] = dataListICU[pid].areaUnderCurve[goodPos-1,:] #AUC
            allFeatureVals[i:i+4,14] = dataListICU[pid].modulationDepth[goodPos-1,:] #ModDepth
            allFeatureVals[i:i+4,15] = dataListICU[pid].secondMoment[goodPos-1,:] #SecondMoment
            i += 4
ICUFeaturesPD = pd.DataFrame(allFeatureVals, columns=[
    'Patient','Group','Camera','Location','Orientation','RelativeOnset','Canopy',
    'Onset','Kurtosis','Skew','ContrastMean','IntensityMean','Amplitude',
    'AUC','ModDepth','SecMom'])

# %%
print(pid,allFeatureVals.shape, totalGoodPos)
# %%
nsubjects = int(np.amax(WaterfallFeaturesPD['Subject']))
npairs = 15
ncameras = int(np.amax(WaterfallFeaturesPD['Camera']))
for position in range(npairs):
    fig, ax = plt.subplots(nrows=1, ncols=4)
    for subject in range(nsubjects):
        for cam in range(ncameras):
            leftIntensity  = WaterfallFeaturesPD.loc[ (WaterfallFeaturesPD['Camera'] == cam) & (WaterfallFeaturesPD['Orientation'] == position) \
                             & (WaterfallFeaturesPD['Subject']==subject) & (WaterfallFeaturesPD['Group'] != 0) ]['IntensityMean'].to_numpy()
            rightIntensity = WaterfallFeaturesPD.loc[ (WaterfallFeaturesPD['Camera'] == cam) & (WaterfallFeaturesPD['Orientation'] == position+15) \
                             & (WaterfallFeaturesPD['Subject']==subject) & (WaterfallFeaturesPD['Group'] != 0) ]['IntensityMean'].to_numpy()
            print(cam, leftIntensity, rightIntensity)
            print('#########')

fig, ax = plt.subplots(nrows=1, ncols=int(cameras))
#ax.format(suptitle=' ', xlabel='Absolute Canopy Difference(L-R)', ylabel='Mean Intensity Difference(L-R)')
fig.text(0.5,0.04, 'Absolute Canopy Difference(L-R)', ha="center", va="center")
fig.text(0.05,0.5, 'Mean Intensity Difference(L-R)', ha="center", va="center", rotation=90)
markerType = ['+','x']
markerColor = ['r','m','y','k']
for cam in range(3,4):
    lfInt = WaterfallFeaturesPD.loc[ (WaterfallFeaturesPD['Camera'] == cam) & (WaterfallFeaturesPD['Location'] == 4) & 
                            (WaterfallFeaturesPD['Orientation'] == 0) & (WaterfallFeaturesPD['Group'] != 0) ] \
                            ['IntensityMean'].to_numpy()
    rfInt = WaterfallFeaturesPD.loc[ (WaterfallFeaturesPD['Camera'] == cam) & (WaterfallFeaturesPD['Location'] == 5) & 
                            (WaterfallFeaturesPD['Orientation'] == 0) & (WaterfallFeaturesPD['Group'] != 0) ] \
                            ['IntensityMean'].to_numpy()
    lfCan = WaterfallFeaturesPD.loc[ (WaterfallFeaturesPD['Camera'] == 0) & (WaterfallFeaturesPD['Location'] == 0) & 
                            (WaterfallFeaturesPD['Orientation'] == 0) & (WaterfallFeaturesPD['Group'] != 0) ] \
                            ['Canopy'].to_numpy()
    rfCan = WaterfallFeaturesPD.loc[ (WaterfallFeaturesPD['Camera'] == 0) & (WaterfallFeaturesPD['Location'] == 1) & 
                            (WaterfallFeaturesPD['Orientation'] == 0) & (WaterfallFeaturesPD['Group'] != 0) ] \
                            ['Canopy'].to_numpy()
    #aa = pd.DataFrame(np.stack(skewDiffByLoc,axis=1 ), columns=pd.Index(l))
    print(lfCan.shape,rfCan.shape,lfInt.shape,rfInt.shape)
    ax[cam].scatter(np.abs(lfCan-rfCan)[:54],(lfInt[:54]-rfInt),color='blue5')
    titleStr = 'Camera '+str(cam+1)
    ax[cam].format(title=titleStr)
    for loc in range(1):
        if loc%2:
            continue
        for pid in range(len(goodScansICU)):
            for orient in range(2):
                lfInt = ICUFeaturesPD.loc[ (ICUFeaturesPD['Camera'] == cam) & (ICUFeaturesPD['Location'] == loc+4) & 
                                        (ICUFeaturesPD['Orientation'] == orient) & 
                                        (ICUFeaturesPD['Patient'] == pid+4) ] ['IntensityMean'].to_numpy() 
                rfInt = ICUFeaturesPD.loc[ (ICUFeaturesPD['Camera'] == cam) & (ICUFeaturesPD['Location'] == loc+5) & 
                                        (ICUFeaturesPD['Orientation'] == orient) & 
                                        (ICUFeaturesPD['Patient'] == pid+4) ] ['IntensityMean'].to_numpy()
                lfCan = ICUFeaturesPD.loc[ (ICUFeaturesPD['Camera'] == cam) & (ICUFeaturesPD['Location'] == loc) & 
                                        (ICUFeaturesPD['Orientation'] == orient) & 
                                        (ICUFeaturesPD['Patient'] == pid+4) ] ['Canopy'].to_numpy() 
                rfCan = ICUFeaturesPD.loc[ (ICUFeaturesPD['Camera'] == cam) & (ICUFeaturesPD['Location'] == loc+1) & 
                                        (ICUFeaturesPD['Orientation'] == orient) & 
                                        (ICUFeaturesPD['Patient'] == pid+4) ] ['Canopy'].to_numpy()
                                        
                if lfInt.size and rfInt.size:
                    if lfInt.size!=rfInt.size or lfCan.size!=lfInt.size:
                        lfInt = np.mean(lfInt); rfInt = np.mean(rfInt)
                        lfCan = np.mean(lfCan); rfCan = np.mean(rfCan)
                    diffInt = lfInt-rfInt; diffCan = lfCan-rfCan
                    ax[cam].scatter(np.abs(diffCan),diffInt,marker=markerType[orient],color=markerColor[pid],markersize=40)


# %%
print(ICUFeaturesPD.loc[ (ICUFeaturesPD['Camera'] == 0) & (ICUFeaturesPD['Location'] == 0) & 
                                        (ICUFeaturesPD['Orientation'] == 0) & 
                                        (ICUFeaturesPD['Patient'] == pid+4) ] ['Canopy'].to_numpy() )
# %%
"""