# %% This cell reads in all the data
from numpy.lib.function_base import delete
from HeadscanData import * #Import all classes
import matplotlib.pyplot as plt
import proplot as pplt
import numpy as np
import pandas as pd
import batchfile, sys

#%%
### file info ###
batch = True
batchname = 'HHC_1B_Schrodinger_Healthy'
#'HHC_1B_Schrodinger_Healthy' 'HHC_1B_Schrodinger_ICU'
#'FrontSt_Schrodinger_20210511and12' 'PulseModulationExperiments' 'WaterfallExperiments'
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
calnames, scannames, datapath = eval('batchfile.' + batchname + '()')
nsubjects = len(scannames)
datapath = os.path.expanduser(datapath)
savepath = datapath + '/Results'

numScansTot = 0
for scanName in scannames:
    numScansTot = numScansTot+len(scanName)

dataList = [[]] * numScansTot
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

### file info ###
batch = True
batchname = 'HHC_1B_Schrodinger_ICU'

calnames, scannames, datapath = eval('batchfile.' + batchname + '()')
nsubjects = len(scannames)
datapath = os.path.expanduser(datapath)
savepath = datapath + '/Results'

numScansTot = 0
for scanName in scannames:
    numScansTot = numScansTot+len(scanName)

dataListICU = [[]] * numScansTot
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
        subjectScan.LoadScanDataAndPreprocess(cropTimeStart=1,cropTimeEnd=-1)
        subjectScan.RunWaveletFilteringOnSignal()
        subjectScan.GetPulseFeatures(lowNoise=True)
        dataListICU[ind] = subjectScan
        subjectIDICU[ind] = ii+1
        ind = ind+1
# %% This cell creates the pandas frames
plotPath = '/Users/kedar/Desktop/plots/'
cameras   = 4
bucket = [0, 0, 2, 4, 4, 1, 1, 3, 5, 5] #Map to waterfall experiment locations - we divide by 3 to remain agnostic to measurement position(1 or 2 or 3)
scanLocNOrient = ['Forehead Vertical','Forehead Horizontal','Fissure','Temple Vertical','Temple Horizontal']
goodScansHealthy = batchfile.HHC_1B_Schrodinger_Healthy_Validation()

totalGoodPos = 0
for volunteer in goodScansHealthy:
    for scan in volunteer:
        totalGoodPos += len(scan)

#Subject, Group, Camera, Location, Oritentation, RelOnset, Canopy, Onset
#Kurtosis, Skewness, ContrastMean, IntensityMean, Amplitude, AUC, ModDepth, SecondMoment
#Group 1-healthy, 10-HHC ICU
allFeatureVals = np.zeros( ((totalGoodPos)*4, 16) )
i = 0
for sid, volunteer in enumerate(goodScansHealthy):
    for scanNo, allGoodPosInScan in enumerate(volunteer):
        ind = sid*3+scanNo
        for goodPos in allGoodPosInScan:
            allFeatureVals[i:i+4,0] = sid #Volunteer
            allFeatureVals[i:i+4,1] = np.asarray([0,1,2,3]) #Camera
            allFeatureVals[i:i+4,2] = bucket[int((goodPos-1)/3)] #Location
            allFeatureVals[i:i+4,3] = goodPos #Position number 1-30 hhc30 order
            allFeatureVals[i:i+4,4] = 1 #Healthy grop - Green color
            allFeatureVals[i:i+4,5] = dataList[ind].pulseOnsetProp[goodPos-1,:] #Relative Onset
            allFeatureVals[i:i+4,6] = dataList[ind].pulseCanopy[goodPos-1,:] #Canopy
            allFeatureVals[i:i+4,7] = dataList[ind].pulseOnset[goodPos-1,:] #Onset
            allFeatureVals[i:i+4,8] = dataList[ind].kurtosis[goodPos-1,:] #Kurtosis
            allFeatureVals[i:i+4,9] = dataList[ind].skewness[goodPos-1,:] #Skewness
            allFeatureVals[i:i+4,10] = dataList[ind].scanData[goodPos-1].contrastAvgsAtPos #ContrastMean
            allFeatureVals[i:i+4,11] = dataList[ind].scanData[goodPos-1].intensityMeansAtPos #IntensityMean
            allFeatureVals[i:i+4,12] = dataList[ind].amplitude[goodPos-1,:] #Amplitude
            allFeatureVals[i:i+4,13] = dataList[ind].areaUnderCurve[goodPos-1,:] #AUC
            allFeatureVals[i:i+4,14] = dataList[ind].modulationDepth[goodPos-1,:] #ModDepth
            allFeatureVals[i:i+4,15] = dataList[ind].secondMoment[goodPos-1,:] #SecondMoment
            i += 4
HHCHealthyFeaturesPD = pd.DataFrame(allFeatureVals, columns=[
    'Subject','Camera','Location','Position','Group',
    'RelativeOnset','Canopy','Onset','Kurtosis','Skew',
    'ContrastMean', 'IntensityMean', 'Amplitude', 'AUC', 'ModDepth', 'SecMom'])

goodScansICU   = batchfile.HHC_1B_Schrodinger_ICU_Validation()
patientIDs =   [1,3,4,5,6,7,8,9,10,11]
patientGroup = [3,3,2,2,2,2,4,2,4 ,3]
#2 - LVO Red ID#6 and ID#9 were major LVOs
#3 - Non-LVO ischemic stroke blue color
#4 - ICH(Hemorrage) Yellow
totalGoodPos = 0
for patient in goodScansICU:
    for scan in patient:
        totalGoodPos += len(scan)
allFeatureVals = np.zeros( (totalGoodPos*4, 16) ) #Patient, Camera, Location, Oritentation, RelOnset, Canopy, Onset
i = 0
for pid, patient in enumerate(goodScansICU):
    for scanNo,scan in enumerate(patient):
        ind = pid*3+scanNo
        for goodPos in scan:
            allFeatureVals[i:i+4,0] = patientIDs[pid] #Patient
            allFeatureVals[i:i+4,1] = np.asarray([0,1,2,3]) #Camera
            allFeatureVals[i:i+4,2] = bucket[int((goodPos-1)/3)] #Location
            allFeatureVals[i:i+4,3] = goodPos #Position number 1-30 hhc30 order
            allFeatureVals[i:i+4,4] = patientGroup[pid] #Group
            allFeatureVals[i:i+4,5] = dataListICU[ind].pulseOnsetProp[goodPos-1,:] #Relative Onset
            allFeatureVals[i:i+4,6] = dataListICU[ind].pulseCanopy[goodPos-1,:] #Canopy
            allFeatureVals[i:i+4,7] = dataListICU[ind].pulseOnset[goodPos-1,:] #Onset
            allFeatureVals[i:i+4,8] = dataListICU[ind].kurtosis[goodPos-1,:] #Kurtosis
            allFeatureVals[i:i+4,9] = dataListICU[ind].skewness[goodPos-1,:] #Skewness
            allFeatureVals[i:i+4,10] = dataListICU[ind].scanData[goodPos-1].contrastAvgsAtPos #ContrastMean
            allFeatureVals[i:i+4,11] = dataListICU[ind].scanData[goodPos-1].intensityMeansAtPos #IntensityMean
            allFeatureVals[i:i+4,12] = dataListICU[ind].amplitude[goodPos-1,:] #Amplitude
            allFeatureVals[i:i+4,13] = dataListICU[ind].areaUnderCurve[goodPos-1,:] #AUC
            allFeatureVals[i:i+4,14] = dataListICU[ind].modulationDepth[goodPos-1,:] #ModDepth
            allFeatureVals[i:i+4,15] = dataListICU[ind].secondMoment[goodPos-1,:] #SecondMoment
            i += 4
ICUFeaturesPD = pd.DataFrame(allFeatureVals, columns=[
    'Subject','Camera','Location','Position','Group','RelativeOnset','Canopy',
    'Onset','Kurtosis','Skew','ContrastMean','IntensityMean','Amplitude',
    'AUC','ModDepth','SecMom'])

#%% Get median among three scans for each valid position and subtract left-right
ICUFeaturesPDMedian = ICUFeaturesPD.groupby(['Subject','Camera','Position']).median().reset_index()
HHCHealthyFeaturesPDMedian = HHCHealthyFeaturesPD.groupby(['Subject','Camera','Position']).median().reset_index()
print( 'Ideal length:', len(goodScansICU)*30*4, ' Got:', ICUFeaturesPDMedian.shape, ' From:', ICUFeaturesPD.shape )
print( 'Ideal length:', len(goodScansHealthy)*30*4, ' Got:', HHCHealthyFeaturesPDMedian.shape, ' From:', HHCHealthyFeaturesPD.shape )

deleteRows = []
for index, row in ICUFeaturesPDMedian.iterrows():
    if row['Position']<=15:
        rowPair = ICUFeaturesPDMedian.loc[(ICUFeaturesPDMedian['Subject']  == row['Subject']) &
                                         (ICUFeaturesPDMedian['Camera']   == row['Camera']) &
                                         (ICUFeaturesPDMedian['Position'] == row['Position']+15)
                                        ]
        if rowPair.empty:
            deleteRows.append(index)
        else:
            if row['Subject'] == 1 or row['Subject'] == 7:
                row.iloc[5:] = row.iloc[5:]-rowPair.iloc[0,5:]
            else:
                row.iloc[5:] = rowPair.iloc[0,5:]-row.iloc[5:]
ICUFeaturesPDMedian = ICUFeaturesPDMedian.drop(ICUFeaturesPDMedian.index[deleteRows])
ICUFeaturesPDMedian = ICUFeaturesPDMedian[ICUFeaturesPDMedian['Position'] <= 15]

deleteRows = []
for index, row in HHCHealthyFeaturesPDMedian.iterrows():
    if row['Position']<=15:
        rowPair = HHCHealthyFeaturesPDMedian.loc[(HHCHealthyFeaturesPDMedian['Subject']  == row['Subject']) &
                                                 (HHCHealthyFeaturesPDMedian['Camera']   == row['Camera']) &
                                                 (HHCHealthyFeaturesPDMedian['Position'] == row['Position']+15)
                                                ]
        if rowPair.empty:
            deleteRows.append(index)
        else:
            row.iloc[5:] = row.iloc[5:]-rowPair.iloc[0,5:]
HHCHealthyFeaturesPDMedian = HHCHealthyFeaturesPDMedian.drop(HHCHealthyFeaturesPDMedian.index[deleteRows])
HHCHealthyFeaturesPDMedian = HHCHealthyFeaturesPDMedian[HHCHealthyFeaturesPDMedian['Position'] <= 15]

print('Final lengths', HHCHealthyFeaturesPDMedian.shape, ICUFeaturesPDMedian.shape )

#%%
from sklearn.tree import DecisionTreeClassifier
def ClassifyTopAndReturnBestClf(classify,llvo,lhealthy):
    bestTp = 0
    fpCur  = 1
    clf = None
    for _ in range(1000):
        clfCur = DecisionTreeClassifier(max_features=2,max_depth=1, max_leaf_nodes=2)
        clfCur.fit(classify.drop(['Group'], axis=1), classify['Group'])
        tp = sum(llvo['Group']==clfCur.predict(llvo.drop(['Group'],axis=1))) / len(llvo)
        fp = sum(lhealthy['Group']!=clfCur.predict(lhealthy.drop(['Group'],axis=1)))/len(lhealthy)
        if not clf or tp>bestTp or (tp==bestTp and fp<fpCur):
            bestTp = tp
            fpCur = fp
            clf = clfCur
    tp = sum(llvo['Group']==clf.predict(llvo.drop(['Group'],axis=1)))/5#len(llvo)
    fp = sum(lhealthy['Group']!=clf.predict(lhealthy.drop(['Group'],axis=1)))/len(lhealthy)
    return tp,fp,clf

# %% Run a simple decision tree classifier and store the classifier with the best TP rate
for cam in range(0,4,3):
    for pos in range(1,16):
        lhealthy = HHCHealthyFeaturesPDMedian.loc[(HHCHealthyFeaturesPDMedian['Camera']   == cam) &
                                                 (HHCHealthyFeaturesPDMedian['Position'] == pos) ].iloc[:,4:]
        llvo = ICUFeaturesPDMedian.loc[(ICUFeaturesPDMedian['Camera']   == cam) &
                                      (ICUFeaturesPDMedian['Position'] == pos) &
                                      (ICUFeaturesPDMedian['Group'] == 2) ].iloc[:,4:]
        lnonlvo = ICUFeaturesPDMedian.loc[(ICUFeaturesPDMedian['Camera']   == cam) &
                                      (ICUFeaturesPDMedian['Position'] == pos) &
                                      (ICUFeaturesPDMedian['Group'] > 2) ].iloc[:,4:]
        classify = pd.concat([lhealthy,llvo])
        llvo = classify.loc[(classify['Group']==2)]
        lhealthy = classify.loc[(classify['Group']==1)]
        print('Position:', pos,' Cam:', cam)
        printPos = False
        for _ in range(2):
            tp,fp,clf = ClassifyTopAndReturnBestClf(classify,llvo,lhealthy)
            if tp>0.5 and fp<0.15:
                printPos = True
            if printPos:
                print('TP Rate:',tp,clf.predict(llvo.drop(['Group'],axis=1)))
                print('FP Rate:',fp,clf.predict(lhealthy.drop(['Group'],axis=1)))
                print('Non LVO ICU:',clf.predict(lnonlvo.drop(['Group'],axis=1)))
            feats = []
            for featNm,featImp in zip(classify.columns[1:],clf.feature_importances_):
                if featImp!=0:
                    feats.append((featNm,featImp))
            if printPos:
                print('Feats:', feats)
            classify = classify.drop(feats[0][0],axis=1)
            lhealthy = lhealthy.drop(feats[0][0],axis=1)
            llvo = llvo.drop(feats[0][0],axis=1)
            lnonlvo = lnonlvo.drop(feats[0][0],axis=1)
        #print(' ')
# %% Plot a few positions
import plotly.express as px
poss = [15,10,9,7,8,5,3]
cams = [3, 3, 3,3,3,0,3]
for i,pos in enumerate(poss):
    cam = cams[i]
    healthy = HHCHealthyFeaturesPDMedian.loc[(HHCHealthyFeaturesPDMedian['Camera'] == cam) & (HHCHealthyFeaturesPDMedian['Position'] == pos) ]
    lvo = ICUFeaturesPDMedian.loc[(ICUFeaturesPDMedian['Camera'] == cam) & (ICUFeaturesPDMedian['Position'] == pos) ]
    plotpd = pd.concat([healthy,lvo])
    tempStr = 'plotpd' + str(pos)
    setattr(sys.modules[__name__], tempStr, plotpd)
print(sys.modules[__name__])
# %%
#Position: 15  Cam: 3
#TP Rate: 0.8 [2. 2. 2. 1. 2.]
#Feats: [('IntensityMean', 1.0)]
pos = 15; cam = 3
fig = px.scatter(plotpd15, x='IntensityMean', y='ContrastMean', color='Group')
fig.show()

#Position: 10  Cam: 3
#TP Rate: 0.6 [2. 1. 2. 1. 2.]
#FP Rate: 0.037037037037037035 [1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 2. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]
#Non LVO ICU: [2. 1. 1.]
#Feats: [('Amplitude', 1.0)]
fig = px.scatter(plotpd10, x='Amplitude', y='Amplitude', color='Group')
fig.show()

#Position: 5  Cam: 0
#TP Rate: 0.6 [2. 1. 2. 1. 2.]
#FP Rate: 0.0 [1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]
#Non LVO ICU: [1. 1.]
#Feats: [('Amplitude', 1.0)]
fig = px.scatter(plotpd5, x='Amplitude', y='Amplitude', color='Group')
fig.show()

#Position: 9  Cam: 3
#TP Rate: 0.6 [1. 2. 2. 1. 2.]
#FP Rate: 0.06666666666666667 [1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 2. 1. 1. 1. 2. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]
#Non LVO ICU: [1. 2.]
#Feats: [('IntensityMean', 1.0)]
#TP Rate: 0.4 [2. 1. 2. 1. 1.]
#FP Rate: 0.03333333333333333 [1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 2. 1. 1. 1.]
#Non LVO ICU: [1. 1.]
#Feats: [('AUC', 1.0)]
fig = px.scatter(plotpd9, x='IntensityMean', y='AUC', color='Group')
fig.show()


# %%
plotList = []
for _,row in plotpd15.iterrows():
    sub = row['Subject']
    gp  = row['Group']
    #if sub==7 and gp==2: #Exclude patient 7
    #    continue
    rowPair5 = plotpd5.loc[(plotpd5['Subject']  == sub ) & (plotpd5['Group'] == gp)]
    if rowPair5.empty:
        amp5 = np.NaN
    else:
        amp5 = rowPair5['Amplitude'].iloc[0]
    rowPair9 = plotpd9.loc[(plotpd9['Subject']  == sub ) & (plotpd9['Group'] == gp)]
    if rowPair9.empty:
        int9 = np.NaN
        auc9 = np.NaN
    else:
        int9 = rowPair9['IntensityMean'].iloc[0]
        auc9 = rowPair9['AUC'].iloc[0]
    rowPair10 = plotpd10.loc[(plotpd10['Subject']  == sub ) & (plotpd10['Group'] == gp)]
    if rowPair10.empty:
        amp10 = np.NaN
    else:
        amp10 = rowPair10['Amplitude'].iloc[0]
    rowPair3 = plotpd3.loc[(plotpd3['Subject']  == sub ) & (plotpd3['Group'] == gp)]
    rowPair7 = plotpd7.loc[(plotpd7['Subject']  == sub ) & (plotpd7['Group'] == gp)]
    rowPair8 = plotpd8.loc[(plotpd8['Subject']  == sub ) & (plotpd8['Group'] == gp)]
    if rowPair3.empty:
        int3 = np.NaN; amp3 = np.NaN
    else:
        int3 = rowPair3['IntensityMean'].iloc[0]; amp3 = rowPair3['Amplitude'].iloc[0]
    if rowPair7.empty:
        int7 = np.NaN; amp7 = np.NaN
    else:
        int7 = rowPair7['IntensityMean'].iloc[0]; amp7 = rowPair7['Amplitude'].iloc[0]
    if rowPair8.empty:
        int8 = np.NaN; amp8 = np.NaN
    else:
        int8 = rowPair8['IntensityMean'].iloc[0]; amp8 = rowPair8['Amplitude'].iloc[0]
    #print(row['IntensityMean'])
    plotList.append(np.asarray([sub,gp,row['IntensityMean'],amp5,amp10,row['ContrastMean'],int9,auc9,
                                int3,amp3,int7,amp7,int8,amp8]))
plotList = np.asarray(plotList)
plotListPD = pd.DataFrame(plotList, columns=['Subject','Group','IntensityMeanP15C3','AmplitudeP5C0',
                        'AmplitudeP10C3','ContrastMeanP15C3','IntensityMeanP9C3','AUCP9C3',
                        'IntensityMeanP3C3','AmplitudeP3C3','IntensityMeanP7C3','AmplitudeP7C3',
                        'IntensityMeanP8C3','AmplitudeP8C3',])

# %%
colours = []
colorList = ['g','r','b','y','k']
labelList = ['IHC','LVO','MIVO','Hmrg','Unknown']
for num in plotList[:,1].astype(int):
    colours.append(colorList[num-1])
plotListPD['Color'] = colours
#print(plotList)

# %%
fig = px.scatter_3d( x=plotList[:,-3], y=plotList[:,-2], z=plotList[:,-1], color=colours)
fig.show()

# %%
import matplotlib.pyplot as plt
fig,ax = plt.subplots()
ax.scatter(x=plotListPD['IntensityMeanP15C3'], y=plotListPD['AmplitudeP5C0'],color=plotListPD['Color'] )
ax.set_xlabel('L-R Difference Intensity Mean Position 15 Camera 3')
ax.set_ylabel('L-R Difference Amplitude Position 5 Camera 3')
for index, row in plotListPD.iterrows():
    print(row)
    if row['Group']==1:
        continue
    ax.annotate(int(row['Subject']), (row['IntensityMeanP15C3'], row['AmplitudeP5C0']))
plt.show()
plt.scatter(x=plotListPD['IntensityMeanP15C3'], y=plotListPD['AmplitudeP5C0'],color=plotListPD['Color'] )
plt.xlabel('L-R Difference Intensity Mean Position 15 Camera 3')
plt.ylabel('L-R Difference Amplitude Position 5 Camera 3')
plt.show()
plt.scatter(x=plotListPD['IntensityMeanP15C3'], y=plotListPD['ContrastMeanP15C3'],color=plotListPD['Color'] )
plt.xlabel('L-R Difference Intensity Mean Position 15 Camera 3')
plt.ylabel('L-R Difference Contrast Mean Position 15 Camera 3')
plt.show()
plt.scatter(x=plotListPD['IntensityMeanP9C3'], y=plotListPD['AUCP9C3'],color=plotListPD['Color'] )
plt.xlabel('L-R Difference Intensity Mean Position 9 Camera 3')
plt.ylabel('L-R Difference AUC Position 9 Camera 3')
plt.show()
# %%
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
clf = LinearDiscriminantAnalysis()
clf.fit( plotListPD.loc[(plotListPD['Group'] < 3)][['IntensityMeanP15C3','AmplitudeP5C0']],plotListPD.loc[(plotListPD['Group'] < 3)]['Group'])
proj1d = clf.transform(plotListPD[['IntensityMeanP15C3','AmplitudeP5C0']])
ldaPd = pd.DataFrame(proj1d,columns=['LDA'])
ldaPd['Group'] = plotListPD['Group']
ldaPd['Color'] = plotListPD['Color']

fig, ax = plt.subplots(ncols=1)
for gp in ldaPd['Group'].unique():
    ldaPd[ldaPd['Group']==gp].LDA.hist(alpha=0.7, ax=ax, color=colorList[int(gp)-1], label=labelList[int(gp)-1], bins=np.arange(-3.2,4.2,0.2))
plt.show()

clf = LinearDiscriminantAnalysis()
clf.fit( plotListPD.loc[(plotListPD['Group'] < 3)][['IntensityMeanP15C3','ContrastMeanP15C3']],plotListPD.loc[(plotListPD['Group'] < 3)]['Group'])
proj1d = clf.transform(plotListPD[['IntensityMeanP15C3','ContrastMeanP15C3']])
ldaPd = pd.DataFrame(proj1d,columns=['LDA'])
ldaPd['Group'] = plotListPD['Group']
ldaPd['Color'] = plotListPD['Color']

fig, ax = plt.subplots(ncols=1)
for gp in ldaPd['Group'].unique():
    ldaPd[ldaPd['Group']==gp].LDA.hist(alpha=0.7, ax=ax, color=colorList[int(gp)-1], label=labelList[int(gp)-1], bins=np.arange(-3.2,4.2,0.2))

clf = LinearDiscriminantAnalysis()
clf.fit( plotListPD.loc[(plotListPD['Group'] < 3)][['IntensityMeanP9C3','AUCP9C3']],plotListPD.loc[(plotListPD['Group'] < 3)]['Group'])
proj1d = clf.transform(plotListPD[['IntensityMeanP9C3','AUCP9C3']])
ldaPd = pd.DataFrame(proj1d,columns=['LDA'])
ldaPd['Group'] = plotListPD['Group']
ldaPd['Color'] = plotListPD['Color']

fig, ax = plt.subplots(ncols=1)
for gp in ldaPd['Group'].unique():
    ldaPd[ldaPd['Group']==gp].LDA.hist(alpha=0.7, ax=ax, color=colorList[int(gp)-1], label=labelList[int(gp)-1], bins=np.arange(-3.2,4.2,0.2))


# %%
goodScanNums = []
goodScPosVol = []
goodScPosICU = []
for patient in  goodScansICU:
    for patientScan in patient:
        goodScanNums.append(len(patientScan))
        goodScPosICU.extend(patientScan)
for subject in  goodScansHealthy:
    for subjectScan in subject:
        goodScanNums.append(len(subjectScan))
        goodScPosVol.extend(subjectScan)

# %% Headset locations - plot intensity mean and amplitude
fig,ax = plt.subplots()
ax.scatter(x=plotListPD['IntensityMeanP3C3'], y=plotListPD['AmplitudeP3C3'],color=plotListPD['Color'] )
ax.set_xlabel('L-R Difference Intensity Mean Position 3 Camera 3')
ax.set_ylabel('L-R Difference Amplitude Position 3 Camera 3')
for index, row in plotListPD.iterrows():
    if row['Group']==1:
        continue
    ax.annotate(int(row['Subject']), (row['IntensityMeanP3C3'], row['AmplitudeP3C3']))
plt.show()

fig,ax = plt.subplots()
ax.scatter(x=plotListPD['IntensityMeanP7C3'], y=plotListPD['AmplitudeP7C3'],color=plotListPD['Color'] )
ax.set_xlabel('L-R Difference Intensity Mean Position 7 Camera 3')
ax.set_ylabel('L-R Difference Amplitude Position 7 Camera 3')
for index, row in plotListPD.iterrows():
    if row['Group']==1:
        continue
    ax.annotate(int(row['Subject']), (row['IntensityMeanP7C3'], row['AmplitudeP7C3']))
plt.show()

fig,ax = plt.subplots()
ax.scatter(x=plotListPD['IntensityMeanP8C3'], y=plotListPD['AmplitudeP8C3'],color=plotListPD['Color'] )
ax.set_xlabel('L-R Difference Intensity Mean Position 8 Camera 3')
ax.set_ylabel('L-R Difference Amplitude Position 8 Camera 3')
for index, row in plotListPD.iterrows():
    if row['Group']==1:
        continue
    ax.annotate(int(row['Subject']), (row['IntensityMeanP8C3'], row['AmplitudeP8C3']))
plt.show()
# %%
