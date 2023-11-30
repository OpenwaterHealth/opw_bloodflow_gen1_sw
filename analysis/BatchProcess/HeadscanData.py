from speckle_fcns_v7 import amplitude
import numpy as np
import re, csv, os, boto3, json, botocore, copy, pywt, math 
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
from scipy.stats import pearsonr, describe, moment
from scipy import signal
from patsy import dmatrix
import statsmodels.api as sm
import matplotlib.pyplot as plt

class HeadScanParams(object):
    """Class to store the scan parameters - device parameters, physical constant, processing flags"""
    
    def __init__(self, scanPath, calPathIn = None, scanOrder = 'hhc30', offsetIn = 0 ):
        """
            Creates a parameter object of HeadScanParams class
            
            Parameters
            ----------
            scanName : str
                Folder on disk where the scan data is stored
            calPathIn : str
                Folder to the disk where the calibration data is stored
            scanOrder : str
                zigzag, circle, hhc30, cuff30
                Pattern used to scan subject
            offsetIn : int
                0=none, 1=ax+b, 2=ax+b+c/x
                offset correction type used to correct bias and non-linearity in the mean vs contrast graph 
        """
        self.n      = 1.4,               # index of refraction
        self.wv     = np.array([8.5e-4]) # wavelength, mm
        self.mua    = np.array([0.09])   # absorption coefficient, mm^-1
        self.musp   = np.array([1.2])    # reduced scattering coefficient, mm^-1
        self.Db     = 1e-6               # diffusion coefficient of the scatterers, mm^2/s
        self.dV2    = 1                  # mean squared velocity of the scatterers, (mm/s)^2
        self.Tmax   = 1e-3               # maximum exposure time being simulated, s
        self.dtau   = 1e-7               # discretization of time step length, s
        self.zb     = 0.1                # extrapolated boundary
        self.K      = 0.12               # it's K the captial version of k *******
        self.cameras2use = np.array([0, 1, 2, 3])
        self.fitBeta = 0
        self.bet     = 0.25
        self.I0      = 1

        self.offsetMode  = offsetIn      #0=none, 1=ax+b, 2=ax+b+c/x

        self.dataPath    = scanPath
        self.calPath     = calPathIn
        self.scannerName = ''

        self.cameraIDs   = []
        self.gain        = []
        self.rho         = []
        self.exposure    = 0
        self.expectedNumDataPoints = 0
        self.dt          = 0
        self.notes       = None

        self.minbpm = 30  # minimum heartbeats per minute to detect
        self.maxbpm = 200 # maximum heartbeats per minute to detect

        if scanOrder == 'zigzag':
            self.scanPosition = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5])
            self.scanLR       = np.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1])
            self.theta        = 3.14*np.array([-0.85, 0.85, -0.7, 0.7, -0.55, 0.55, -0.42, 0.42, -0.28, 0.28, -0.12, 0.12])
        elif scanOrder == 'circle':
            self.scanPosition = np.array([0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0])
            self.scanLR       = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1])
            self.theta        = 3.14*np.array([-0.85, -0.7, -0.55, -0.42, -0.28, -0.12, 0.12, 0.28, 0.42, 0.55, 0.7, 0.85])
        elif scanOrder == 'hhc30':
            self.scanPosition = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
            self.scanLR       = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
            self.theta        = np.nan
        elif scanOrder == 'cuff30':
            self.scanPosition = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14])
            self.scanLR       = np.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1])
            self.theta        = np.nan

        self.nLocations = len(self.scanPosition)


class SingleScanData(object):
    """Class that stores the data from one location in the scan"""
    def __init__(self):
        """
            Initializes an object for the class that stores the data from one location in the scan
        """
        self.valid         = True
        self.numDataPoints = 0
        #Corresponds to the list of column headers in the CSV file
        self.csvParams     = ['camera', 'frame', 'timeStamp', 'saturated', 'rawMean', 'rawStdDev',
                              'mean', 'stdDev', 'contrast', 'temperature']
        #for param in self.csvParams: #Compact way
        #        setattr(self, param, 0)
        #Adding above list explcitly so that auto complete in IDEs pick up the names
        self.camera = 0; self.frame = 0; self.timeStamp = 0; self.saturated = 0; self.rawMean = 0
        self.rawStdDev = 0; self.mean = 0; self.stdDev = 0; self.contrast = 0; self.temperature = 0
        self.intensityMeansAtPos   = 0 #mean_mean
        self.stdAvgsAtPos          = 0 #std_mean
        self.contrastAvgsAtPos     = 0 #contrast_mean
        self.darkMeansAtPos        = 0 #mean_dark_mean
        self.darkStdAvgsAtPos      = 0 #std_dark_mean
        self.darkContrastAvgsAtPos = 0 #contrast_dark_mean 

        #Can probably remove the following once we have confidence in calibrations
        self.meanDark     = 0
        self.stdDevDark   = 0
        self.contrastDark = 0
        self.allParams    = self.csvParams.copy()
        self.allParams.extend( ['meanDark','stdDevDark','contrastDark'] )
        
    def ReadBloodFlowProtoCSV2(self,fname):
        """
            Reads the blood flow data from a single location scan
            
            Parameters
            ----------
            fname : str
                Path to csv file on disk where the scan data is stored
        """
        # initialize variables
        with open(fname, mode='r') as csvfile:
            readCSV = csv.reader(csvfile, delimiter=',')
            self.numDataPoints = 0
            for row in readCSV:
                self.numDataPoints += 1
        self.numDataPoints = self.numDataPoints - 1
        
        for param in self.csvParams:
                setattr(self, param, np.zeros((self.numDataPoints,)))
                
        # read data
        with open(fname) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=',')

            for i, row in enumerate(readCSV):
                if i>0 and len(row)==10:
                    for index, param in enumerate(self.csvParams):
                        curArray = getattr(self,param)
                        curArray[i-1] = row[index] #row 0 typically has the column headers
                elif i>0:
                    print('CSV Format Has Changed!!!')

    def SortByCameraAndTimeStampProto2(self):
        """
            Sorts and splits the data by camera ID, then sorts by time stamp so that light/dark frames are consistent
        """

        allCameras, numMeasurementsPerCam = np.unique(self.camera, return_counts=True)
        numCameras = np.size(allCameras)
        numMeasurements = np.amax(numMeasurementsPerCam)

        #Collect the indices for the camera
        allCamInds = []
        for i, cam in enumerate(allCameras):
            allCamInds.append( np.argwhere(self.camera == cam) )
        self.numDataPoints = numMeasurements

        #Slice data by camera
        for param in self.csvParams:
            tempArray = np.zeros((numMeasurements, numCameras))
            for i, camInd in enumerate(allCamInds):
                if numMeasurementsPerCam[i] == numMeasurements:
                    tempArray[:,i] = np.squeeze( getattr(self, param)[camInd] )
                else:
                    tempArray[:numMeasurementsPerCam[i],i] = np.squeeze( getattr(self, param)[camInd] )
                    tempArray[numMeasurementsPerCam[i]:,i] = np.nan
            setattr(self, param, tempArray)
        
        #Sort data by time stamp
        sortIndex = np.argsort(self.timeStamp, axis=0)
        for param in self.csvParams:
            setattr(self, param, np.take_along_axis(getattr(self, param), sortIndex, axis=0 ))

    def RemoveDark2(self):
        """
            Subtracts the mean and stdDev between alternate frames (laser is pulsed on/off)
        """
        endIdx = self.rawMean.shape[0]-self.rawMean.shape[0]%2

        self.mean   = np.abs(self.rawMean[1:endIdx:1, :] - self.rawMean[0:endIdx-1:1, :])
        self.mean   = (self.mean[1::2, :] + self.mean[0:-1:2, :])/2
        self.mean   = np.vstack((self.mean, np.squeeze(np.abs(self.rawMean[-2:-1, :] - self.rawMean[-1:, :]))))

        self.stdDev = np.sqrt(np.abs(self.rawStdDev[1:endIdx:1, :]**2 - self.rawStdDev[0:endIdx-1:1, :]**2))
        self.stdDev = (self.stdDev[1::2, :] + self.stdDev[0:-1:2, :])/2
        self.stdDev   = np.vstack((self.stdDev, np.squeeze(np.sqrt(np.abs(self.rawStdDev[-2:-1, :]**2 - self.rawStdDev[-1:, :] ** 2)))))
        
        self.cameraMeansAtLocation = np.mean(self.mean, axis=0)
        
        self.contrast = self.stdDev / self.mean
        
        #Alternate even indices to correspond to bright frames
        self.camera      = self.camera[1:endIdx:2, :]
        self.timeStamp   = self.timeStamp[1:endIdx:2, :]
        self.temperature = self.temperature[1:endIdx:2, :]

        #Can probably remove the following once we have confidence in calibrations
        self.meanDark     = self.rawMean[0:endIdx:2, :]
        self.stdDevDark   = self.rawStdDev[0:endIdx:2, :]
        self.contrastDark = self.stdDevDark/self.meanDark

    def RemoveDark2PulseModulation(self):
        """
            Subtracts the mean and stdDev between four frames - laser is off, normal pulse is second, laser is off for the third and
            a decoherent pulse is the fourth in the repeated cycle
        """
        endIdx = self.rawMean.shape[0]-self.rawMean.shape[0]%4

        #Get index for modulated pulse and dark indices - normal pulse should have the highest std dev
        sortedInds = np.argsort(self.rawStdDev[0:4][:,0])

        self.mean   = np.abs(self.rawMean[sortedInds[-1]:endIdx:4,:] - self.rawMean[sortedInds[0]:endIdx:4,:])
        self.stdDev = np.sqrt(np.abs(self.rawStdDev[sortedInds[-1]:endIdx:4,:]**2 - self.rawStdDev[sortedInds[0]:endIdx:4,:]**2))
        self.stdDevPM = np.sqrt(np.abs(self.rawStdDev[sortedInds[-1]:endIdx:4,:]**2 - self.rawStdDev[sortedInds[-2]:endIdx:4,:]**2))
        
        self.cameraMeansAtLocation = np.mean(self.mean, axis=0)
        
        self.contrastPM = self.stdDevPM / self.mean
        self.contrast   = self.stdDev / self.mean
        
        #Alternate even indices to correspond to bright frames
        self.camera      = self.camera[1:endIdx:4, :]
        self.timeStamp   = self.timeStamp[1:endIdx:4, :]
        self.temperature = self.temperature[1:endIdx:4, :]

        #Can probably remove the following once we have confidence in calibrations
        self.meanDark     = self.rawMean[0:endIdx:4, :]
        self.stdDevDark   = self.rawStdDev[0:endIdx:4, :]
        self.contrastDark = self.stdDevDark/self.meanDark
    
    def AverageTimeSeriesProto(self):
        """
            Computes the moments for a subset of the measurements
        """
        self.intensityMeansAtPos   = np.nanmean(self.mean, axis=0)
        self.stdAvgsAtPos          = np.nanmean(self.stdDev, axis=0)
        self.contrastAvgsAtPos     = np.nanmean(self.contrast, axis=0)
        self.darkMeansAtPos        = np.nanmean(self.meanDark, axis=0)
        self.darkStdAvgsAtPos      = np.nanmean(self.stdDevDark, axis=0)
        self.darkContrastAvgsAtPos = np.nanmean(self.contrastDark, axis=0)

    def CropTimeSeries(self, cropStart = None, cropEnd = None):
        """
            Discards the first and last few measurements .

            Parameters
            ----------
            cropStart : int
                Keeps data starting from this index. Leave empty to keep all data.
            cropEnd : int
                Keeps data up to this index. Leave empty to keep all data.
        """
        if cropEnd and cropStart:
            for param in self.allParams:
                setattr( self, param, getattr(self,param)[cropStart:cropEnd,:] )
                return
        if cropEnd:
            for param in self.allParams:
                setattr( self, param, getattr(self,param)[:cropEnd,:] )
        if cropStart:
            for param in self.allParams:
                setattr( self, param, getattr(self,param)[cropStart:,:] )


class HeadScanData(object):
    """
        Class to handle the data from a scan and the processing functions asscoiated with it
    """
    def __init__(self, scanParamIn:HeadScanParams, displayFitsIn=False):
        """
            Creates a an object to handle a head scan data
            
            Parameters
            ----------
            scanParamIn : HeadScanParams
                Initialize a headscan parameters object with scanner name, data paths for scan
                and calibration path and pass it into this class
            displayFitsIn : bool
                Flag to enable flow fit plots
        """
        self.completedScan = 0
        self.nAttempts     = 0
        self.displayFits   = displayFitsIn
        self.scanParams    = scanParamIn
        self.scanData      = None
        self.goldenPulse   = {}
        self.pulseDeteced  = None
        self.pulseRate     = None
        self.areaUnderCurve = None
        self.amplitude      = None
        self.modulationDepth= None
        self.skewness       = None
        self.kurtosis       = None
        self.pulseCanopy    = None
        self.pulseOnset     = None
        self.pulseOnsetProp = None

    def ReadLogFile(self):
        """
            Reads the log file and checks if the individual location scans have succeded
            along with the number of attempts
        """
        logFile = os.path.join(self.scanParams.dataPath,'log.txt')
        myFile  = open(logFile, "r")
        content = myFile.read()
        myFile.close()
        contentList = content.split("\n")

        location = 0
        self.completedScan = [True]*self.scanParams.nLocations
        self.nAttempts     = [0]*self.scanParams.nLocations

        for i in range(len(contentList)):
            if re.search('Capture at point', contentList[i]):
                location +=1
            if re.search('Capture complete', contentList[i]):
                self.nAttempts[location-1] +=1
            if re.search('WARNING: Location', contentList[i]):
                print(contentList[i])
                self.completedScan[location-1] = False

    def GetJSONInfo(self):
        """
            Reads the JSON info file - number of cameras, images, gain, rho(distance between cams)
            exposure(pulse width)
        """
        ### load json file ###
        with open(os.path.join(self.scanParams.dataPath , 'scan_metadata.json')) as jsonFile:
            jsonDict = json.load(jsonFile)
        numCams = jsonDict.get('cameraParameters').get('numCameras')
        scannerCameras = list(jsonDict.get('cameraParameters').get('cameraInfo').keys())
        self.scanParams.scannerName = list(jsonDict.get('cameraParameters').get('cameraInfo').keys())[0]
        for i in range(len(scannerCameras)):
            if scannerCameras[i] == self.scanParams.scannerName:
                cids = list(jsonDict.get('cameraParameters').get('cameraInfo').keys())[i+1:i+5]
                self.scanParams.cameraIDs = cids
        gain = np.zeros((numCams,))
        rho = np.zeros((numCams,))
        for i in range(numCams):
            gain[i]=jsonDict.get('cameraParameters').get('cameraInfo').get(self.scanParams.cameraIDs[i]).get('gain')
            rho[i]=jsonDict.get('cameraParameters').get('cameraInfo').get(self.scanParams.cameraIDs[i]).get('separation_mm')
        self.scanParams.gain     = gain
        self.scanParams.rho      = rho
        self.scanParams.expectedNumDataPoints = jsonDict.get('cameraParameters').get('numImages')
        self.scanParams.exposure = np.array([jsonDict.get('delayParameters').get('pulseWidth_s')])
        self.scanParams.dt       = 2/jsonDict.get('cameraParameters').get('frameAcquisitionRate_Hz')
        self.scanParams.notes    = jsonDict.get('sampleParameters').get('experimentNotes')

    def S3CheckAndDownloadDataset(self, folderName):
        '''
            Uses credentials stored in the profile - download s3 cli utility and setup for the computer
            Checks the local disk for folder and downloads from S3 if it is not downloaded 
            Assumes the path contains BucketName/DeviceName/DataFolder are in the path passed 
            
            Parameters
            ----------
            folderName : string or os path
                Path should end in BucketName/DeviceName/DataFolder
        '''

        logFile = os.path.join(folderName,'log.txt')
        if os.path.exists(logFile):
            return #Data exists on disk - nothing to do
        
        try:
            deviceFolder, recordingName = os.path.split(folderName)
            bucketFolder, deviceName    = os.path.split(deviceFolder)
            rootFolder,   bucketName    = os.path.split(bucketFolder)
            if not os.path.isdir(rootFolder):
                os.makedirs(rootFolder)
            if not os.path.isdir(bucketFolder):
                os.makedirs(bucketFolder)
            if not os.path.isdir(folderName):
                os.makedirs(folderName)
           
            #Need a resource to download
            s3 = boto3.resource('s3')
            
            #Need a client to paginate(search s3)
            s3client = boto3.client('s3')
            paginator = s3client.get_paginator('list_objects')
            operation_parameters = {'Bucket': bucketName, 'Prefix': deviceName+'/'+recordingName}
            page_iterator = paginator.paginate(**operation_parameters)
            
            for page in page_iterator:
                for fileItem in page['Contents']:
                    fileName  = bucketFolder+'/'+fileItem['Key']
                    writePath = os.path.expanduser(fileName)
                    writeDir, _ = os.path.split(writePath)
                    if not os.path.isdir(writeDir): #Skipping subfolders for now - need to create on disk before writing
                        continue
                    s3.Bucket(bucketName).download_file(fileItem['Key'], writePath)

        except botocore.exceptionsNoCredentialsError as error:
            print(error)
            print('Either fix the aws credentials or download the dataset ' + folderName + 'to run this script.')

    def CheckDataFromCamsAndDeleteDroppedCamParams(self):
        '''
            Checks the dataset for dropped cameras and deletes the camera info from the data structures so that
            they don't factor into subsequent analysis
        '''
        allCameraIDs = [ np.unique(scan.camera[~np.isnan(scan.camera)]) for scan in self.scanData if scan ]
        camerasPresent = np.unique(allCameraIDs)
        droppedCameras = np.setdiff1d( [int(id) for id in self.scanParams.cameraIDs], camerasPresent.astype(int) )
        dropCamInds    = []
        for dropCam in droppedCameras:
            dropCamInds.append(self.scanParams.cameraIDs.index(str(dropCam)))
            self.scanParams.cameraIDs.remove(str(dropCam))
        dropCamInds.sort(reverse=True)
        #JSON may have correct data and gain and rho may not be populated so check sizes before deleting
        for i in dropCamInds:
            if self.scanParams.cameras2use.shape[0] > len(self.scanParams.cameraIDs):
                self.scanParams.cameras2use = np.delete(self.scanParams.cameras2use, i)
            if self.scanParams.gain.shape[0] > len(self.scanParams.cameraIDs):
                self.scanParams.gain        = np.delete(self.scanParams.gain, i)
            if self.scanParams.rho.shape[0]  > len(self.scanParams.cameraIDs):
                self.scanParams.rho         = np.delete(self.scanParams.rho, i)
    
    def LoadScanDataAndPreprocess(self, cropTimeStart = None, cropTimeEnd = None, pulseModulation = False):
        '''
            Main function to load the scan data do the initial processing. Functions included are:
            1. Check if data is on disk and download if needed
            2. Load calibration data
            3. Load individual scans data sorted by camera and time stamp
            4. Subtract dark frames
            5. Check for dropped cameras
            6. Calculate time series averages
            
            Parameters
            ----------
            cropTimeStart(optional) : int
                Crops the loaded time series starting with this index
            
            cropTimeEnd(optional) : int
                Crops the loaded time series and ends with this index

            waveletPreporcess(optiona): bool
        '''
        #Check and dowload data needed for processing
        if self.scanParams.calPath:
            self.S3CheckAndDownloadDataset(self.scanParams.calPath)
        self.S3CheckAndDownloadDataset(self.scanParams.dataPath)

        ### Read json and log files ###
        self.GetJSONInfo()
        self.ReadLogFile()

        ### load calibration data ###
        if self.scanParams.calPath: #I0 is initialized to 1 by default
            calScanData = SingleScanData()
            calScanData.ReadBloodFlowProtoCSV2(self.scanParams.calPath + '/data.csv')
            calScanData.SortByCameraAndTimeStampProto2()
            calScanData.RemoveDark2()
            calScanData.AverageTimeSeriesProto()
            self.scanParams.I0 = calScanData.intensityMeansAtPos
        
        ### load data ###
        self.scanData = [[]] * self.scanParams.nLocations
        csvNames = [ 'location_' + str(i+1) + '.csv' for i in range(self.scanParams.nLocations) ]
        for i,fileName in enumerate(csvNames):
            currentScan = SingleScanData()
            if not os.path.exists(os.path.join(self.scanParams.dataPath,fileName)) or not self.completedScan[i]:
                currentScan.valid     = False
                self.completedScan[i] = False #Set if csv is missing but says it is completed in the log
                continue

            currentScan.ReadBloodFlowProtoCSV2(os.path.join(self.scanParams.dataPath,fileName))
            currentScan.SortByCameraAndTimeStampProto2()
            if pulseModulation:
                currentScan.RemoveDark2PulseModulation()
            else:
                currentScan.RemoveDark2()
            currentScan.AverageTimeSeriesProto()

            if cropTimeStart or cropTimeEnd:
                currentScan.CropTimeSeries(cropTimeStart,cropTimeEnd)

            self.scanData[i] = currentScan

        self.CheckDataFromCamsAndDeleteDroppedCamParams()

    #####################################################
    ### Functions to process acquired data start here ###
    #####################################################
    def LowPassWaveletFilter(self, signal, thresh = 0.1, wavelet="sym2", minBand=3):
        '''
            Filters about the threshold fraction of the energy from the wavelets with support smaller than the band
            passed in with minBand
            
            Parameters
            ----------
            signal : array
                Signal is a 2D array with channels in the order (samples,channels)
            
            thresh : float
                Percentage of the energy to remove from the wavelets with smaller signal support(higher frequencies)
                (default is 0.1 - this is usually the signal with some high frequency noise)

            wavelet : string
                Wavelet type to use in the processing.
                (default is sym2 - correlates best with our pulse shape)
            
            minBand: int
                Wavelet band from which to start thresholding.
            
            Returns
            -------
            reconstructedSignal: array
                Returned in a 2D array the same shape as the input channel data
        '''
        reconstructedSignal = np.empty_like(signal)
        for i in range(signal.shape[1]):
            thresh = thresh*np.nanmax(signal[:,i])
            coeff = pywt.wavedec(signal[:,i], wavelet, mode="per" )
            coeff[minBand:] = (pywt.threshold(i, value=thresh, mode="soft" ) for i in coeff[minBand:])
            reconstructedSignal[:,i] = pywt.waverec(coeff, wavelet, mode="per" )
        return reconstructedSignal

    def HighPassWaveletFilter(self, signal, wavelet='sym2'):
        '''
            Filter removes about the threshold fraction of the energy from the wavelet decomposition for
            wavelets which have support up to 1/4th the length of the signal. The threshold for each channel
            is a factor of the standard deviation of the coefficients in the decomposed bands.
            
            Parameters
            ----------
            signal : array
                Signal is a 2D array with channels in the order (samples,channels)
            
            wavelet : string
                Wavelet type to use in the processing.
                (default is sym2 - correlates best with our pulse shape)

            Returns
            -------
            reconstructedSignal: array
                Returned in a 2D array the same shape as the input channel data
            
            lowPassCoeffStart: int
                The band up to which the coefficients were thresholded
        '''
        coeffsChs = pywt.wavedec (signal, wavelet, axis=0)
        lowPassCoeffStart = 0
        for i in range(len(coeffsChs)):
            if (coeffsChs[i].shape[0]) > math.log2(signal.shape[0])*1.5:
                lowPassCoeffStart = i
                break
        stdDev = np.zeros(signal.shape[1])
        for i in range(lowPassCoeffStart):
            stdDev = np.maximum(stdDev,np.nanstd(coeffsChs[i], axis=0))
        threshCh = np.clip(stdDev*10,0.1,0.95)

        for i in range(lowPassCoeffStart):
            for ch in range(signal.shape[1]):
                coeffsChs [i][:,ch] = pywt.threshold (coeffsChs[i][:,ch], threshCh[ch] * max (coeffsChs[i][:,ch]) )
        reconstructedSignal = pywt.waverec (coeffsChs, wavelet, axis=0)
        return reconstructedSignal, lowPassCoeffStart
    
    def RunWaveletFilteringOnSignal(self, lowPassThreshod = 0.1, attribute = 'contrast'):
        '''
            Runs the wavelet thresholding on the signal. Divided into high and low bands at wavelet support of roughly
            1/4th the length of the signal. High pass thresholding is applied by using the 
            
            Parameters
            ----------
            lowPassThreshod : float
                Percentage of the energy to remove from the wavelets with smaller signal support(higher frequencies)
                (default is 0.1 - this is usually the signal with some high frequency noise)
            
            attribute : string
                Atrribute of scans on which the thresholding should be applied
                (default is the contrast)
        '''
        for scan in self.scanData:
            if not scan:
                continue
            arrayForProcessing = getattr(scan, attribute)
            signalMean = np.nanmean(arrayForProcessing,axis=0)
            filteredSignalHp, highBandStart = self.HighPassWaveletFilter(arrayForProcessing-signalMean)
            filteredSignal = self.LowPassWaveletFilter(filteredSignalHp, lowPassThreshod, minBand=highBandStart)
            setattr(scan,attribute,filteredSignal+signalMean)

    def RemoveShotNoise( self ):
        '''
            Runs the shot noise subtraction model 
        '''
        for i, scan in enumerate(self.scanData):
            if not self.completedScan[i]:
                continue
            scan.stdDev   = np.sqrt( scan.stdDev**2 - self.scanParams.K * self.scanParams.gain * scan.mean )
            scan.contrast = scan.stdDev/scan.mean
            np.nanmean()

    def GetPulseFeatures(self, lowNoise = False, attribute = 'contrast'):
        '''
            Runs pulse detection and averages the pulses together in each camera to get the exemplar pulse in each
            camera

            Parameters
            ----------
            lowNoise : bool
                Assumes low noise and tries to find every pulse in a scan with a Butterworth low pass filter.
                Can lead to incorrect results if the channel inputs are noisy.

            attribute : string
                Atrribute of scans on which the thresholding should be applied
                (default is the contrast)
        '''
        #Initialize structures
        self.goldenPulse    = [[]] * self.scanParams.nLocations
        self.areaUnderCurve = np.zeros((self.scanParams.nLocations, len(self.scanParams.cameraIDs)))
        self.amplitude      = np.zeros((self.scanParams.nLocations, len(self.scanParams.cameraIDs)))
        self.modulationDepth= np.zeros((self.scanParams.nLocations, len(self.scanParams.cameraIDs)))
        self.skewness       = np.zeros((self.scanParams.nLocations, len(self.scanParams.cameraIDs)))
        self.kurtosis       = np.zeros((self.scanParams.nLocations, len(self.scanParams.cameraIDs)))
        self.secondMoment   = np.zeros((self.scanParams.nLocations, len(self.scanParams.cameraIDs)))
        self.pulseDeteced   = np.zeros((self.scanParams.nLocations, len(self.scanParams.cameraIDs)), dtype=bool)
        self.pulseCanopy    = np.zeros((self.scanParams.nLocations, len(self.scanParams.cameraIDs)))
        self.pulseOnset     = np.zeros((self.scanParams.nLocations, len(self.scanParams.cameraIDs)))
        self.pulseOnsetProp = np.zeros((self.scanParams.nLocations, len(self.scanParams.cameraIDs)))
        self.pulseRate      = np.full((self.scanParams.nLocations, len(self.scanParams.cameraIDs)),np.nan)

        noPulseCount = 0
        fewPulseCount = 0
        for loc,scan in enumerate(self.scanData):
            if not scan:
                continue
            y = getattr(scan, attribute)

            allNan = False
            for ch in y.T:
                if np.isnan(ch).all():
                    allNan = True
            if allNan:
                print('No waveform at scan location ', loc+1, ': Skipping bad data')
                self.goldenPulse[loc] = [ np.nan*np.zeros((10,)) for i in range(y.shape[1]) ]
                self.pulseDeteced[loc,:] = periodAll != 0
                continue

            N = int(np.power(2,np.ceil(np.log2(y.shape[0]))))  # nearest power of 2 length for FFT
            
            ### approximate period by finding biggest peak in Fourier spectrum ###
            m = np.mean(y, axis=0)
            Y_hat = np.abs(np.fft.fftshift(np.fft.fft(y-m, n=N, axis=0), axes=0))
            mini = np.argmin(Y_hat, axis=0) #Will be freq = 0 because DC is removed before FFT
            maxi = np.argmax(Y_hat, axis=0)
            periodAll = N/np.abs(mini-maxi)
            bpm = 60/self.scanParams.dt/periodAll
            periodAll[bpm<self.scanParams.minbpm] = 0
            periodAll[bpm>self.scanParams.maxbpm] = 0
            period = np.median(periodAll)
            if period == 0 or min(periodAll) == 0:
                print('No waveform at scan location ', loc+1, ': Skipping bad data')
                self.goldenPulse[loc] = [ np.nan*np.zeros((10,)) for i in range(y.shape[1]) ]
                self.pulseDeteced[loc,:] = periodAll != 0
                continue
            
            self.pulseRate[loc,:] = bpm
            self.pulseDeteced[loc,:] = periodAll != 0
            tempInds = np.where( periodAll == 0 )
            self.pulseRate[loc,tempInds] = np.nan

            ### separate periods ###
            allChFits = []
            if lowNoise:
                fewPulseCountCur , noPulseCountCur, allChFits = self.FindLowNoisePulseAndAverage(y, period, loc)
            else:
                fewPulseCountCur , noPulseCountCur = self.FindNoisyPulseAndAverage(y, period, loc)
            noPulseCount += noPulseCountCur
            fewPulseCount += fewPulseCountCur
            for ch in range(y.shape[1]):
                if np.isnan(self.goldenPulse[loc][ch]).all():
                    continue
                gpDistribution = self.goldenPulse[loc][ch]-np.amin(self.goldenPulse[loc][ch])
                gpDistribution = np.amax(gpDistribution)-gpDistribution
                gpDistribution = gpDistribution / np.sum(gpDistribution)
                self.areaUnderCurve[loc,ch] = self.GetAreaUnderCurve(np.copy(gpDistribution))
                self.amplitude[loc,ch]      = np.amax(self.goldenPulse[loc][ch]) - np.amin(self.goldenPulse[loc][ch])
                self.modulationDepth[loc,ch]= self.amplitude[loc,ch]/np.mean(self.goldenPulse[loc][ch])
                gpDescription  = describe(gpDistribution)
                self.skewness[loc,ch] = gpDescription.skewness
                self.kurtosis[loc,ch] = gpDescription.kurtosis
                if allChFits and allChFits[ch]:
                    canopy, onset, onsetProp, secMoment = \
                        self.ComputeWaveformAttributesFromFit(allChFits[ch],self.goldenPulse[loc][ch].shape[0])
                else:
                    canopy, onset, onsetProp, secMoment = \
                        self.ComputeWaveformAttributes(self.goldenPulse[loc][ch])
                self.pulseCanopy[loc,ch]    = canopy
                self.pulseOnset[loc,ch]     = onset
                self.pulseOnsetProp[loc,ch] = onsetProp
                self.secondMoment[loc,ch] = secMoment

        if fewPulseCount+noPulseCount:
            print('Number of camera locations where few(1 or 2) pulses were recovered is:',fewPulseCount,' no pulses recovered:', noPulseCount)

    def FindLowNoisePulseAndAverage(self, y, period, loc):
        '''
            Runs pulse detection and crops, then fits a spline to the pulses in each camera to get the exemplar pulse in each
            camera

            Parameters
            ----------
            y : 2D array
                First dimension is the time the second is the number of channels/cameras

            period : float
                Period(in number of samples) of pulse detected in the signal

            loc : int
                Index of location where y(the signal being processed) was acquired
            
            Returns
            -------
            fewPulseCount: int
                Count where fewer than three pulses are detected in a channel

            noPulseCount: int
                Count where no pulse is detected in a channel

            fits: list of arrays
                Coefficients for the spline fit
        '''
        noPulseCount  = 0
        fewPulseCount = 0
        ### separate periods ###
        ##Take away all the higher frequency components from the signal to help find the start of systolic
        bpm = 60/self.scanParams.dt/period
        cut   = np.median(bpm)/60
        soslp = signal.butter(2,cut,'lp',fs=1/self.scanParams.dt,output='sos')
        m = np.mean(y, axis=0)
        ySmoothed = signal.sosfilt(soslp,y-m,axis=0)
        minPeriod = int(np.floor(period))-int(period*0.15) #allow for a 15% variation decrease in pulse width
        chPeaksSm = [find_peaks(ySmoothed[:,ch], distance=minPeriod)[0] for ch in range(ySmoothed.shape[1])]

        #Make the arrays the same size so that we can use numpy tools on the array
        lenPeaks  = max(map(len,chPeaksSm))
        for tmpInd1, pkList in enumerate(chPeaksSm):
            if len(pkList) == lenPeaks:
                mxLenPeakInd  = tmpInd1
                break
        for tmpInd1 in range(len(chPeaksSm)):
            peaks = chPeaksSm[tmpInd1]
            if len(peaks) < lenPeaks:
                #Append to the correct side - when we have even number of cameras, median can split
                #the difference between the lists and yeild incorrect waveforms
                if np.sum(np.abs(peaks-chPeaksSm[mxLenPeakInd][:len(peaks)])) < \
                    np.sum(np.abs(peaks-chPeaksSm[mxLenPeakInd][-len(peaks):])):
                    for _ in range(lenPeaks-len(peaks)):
                        peaks = np.append(peaks.astype(float),np.nan)
                else:
                    for _ in range(lenPeaks-len(peaks)):
                        peaks = np.insert(peaks.astype(float),0,np.nan)
                chPeaksSm[tmpInd1] = peaks
        
        #Heartbeat should be consistent across cameras
        hbPeaks = np.nanmedian(np.asarray(chPeaksSm),axis=0)
        hbPeaks = hbPeaks.astype(int)
        
        #We search for the systolic peak around the start of the 
        yPrime    = y[1:,:] - y[:-1,:]
        #Primary method to find pulse
        probSystolic1 = [ find_peaks(yPrime[:,yP]**2,distance=minPeriod)[0] for yP in range(yPrime.shape[1]) ]
        #Backup method to find pulse
        probSystolic2 = [ find_peaks(yPrime[:,yP]**2)[0] for yP in range(yPrime.shape[1]) ]
        segStarts = np.zeros((ySmoothed.shape[1],len(hbPeaks)),dtype='i')
        for ind in range(len(hbPeaks)):
            inds1 = probSystolic1 - hbPeaks[ind]
            inds2 = probSystolic2 - hbPeaks[ind]
            for ch,chInds1 in enumerate(inds1):
                chInds1[chInds1>0] = -ySmoothed.shape[0]
                if np.abs(np.amax(chInds1))<5:
                    segStarts[ch,ind] = probSystolic1[ch][np.argmax(chInds1)]
                else:
                    chInds2 = copy.deepcopy(inds2[ch])
                    chInds2[chInds2>0] = -ySmoothed.shape[0]
                    segStarts[ch,ind] = probSystolic2[ch][np.argmax(chInds2)]

        segStarts -=1 #We want the index where the pulse starts
        chPeaks   = [ segSt for segSt in segStarts ]

        chLengths = [peaks[1:]-peaks[:-1] for peaks in chPeaks]

        ### get segments only pick ones that are close < +/-2 of approximation
        chLengthsLst = [l.tolist() for l in chLengths]
        chLengthsLst = sum(chLengthsLst, [])
        maxLength = np.ceil(np.median(chLengthsLst)+1.5)
        minLength = np.floor(np.median(chLengthsLst)-1.5)

        L = int(maxLength+1)

        chUseSegment = [ ( (minLength<=np.array(lengths)).astype(int) + (maxLength>=np.array(lengths)).astype(int) ) == 2
                        for lengths in chLengths ]

        self.goldenPulse[loc] = [{}] * len(chUseSegment)

        fits = [[]] * len(chUseSegment)
        for ch in range(len(chUseSegment)):
            nSegments =int(np.sum(chUseSegment[ch]))
            if nSegments == 0:
                self.goldenPulse[loc][ch] = np.nan*np.zeros((10,))
                self.pulseDeteced[loc,ch] = False
                noPulseCount += 1
                continue
            segments = np.zeros((L, nSegments))
            if nSegments<3:
                fewPulseCount += 1
            ind=0
            for j in range(len(chLengths[ch])):
                if not chUseSegment[ch][j]:
                    continue
                curSeg = y[chPeaks[ch][j]:chPeaks[ch][j+1]+1, ch]
                if curSeg[0]<curSeg[1]: #We might be off by 1 index from systolic start
                    curSeg = np.roll(curSeg, -1)
                segments[:chLengths[ch][j]+1, ind] = curSeg
                ind += 1
            segments[segments == 0] = np.nan
            
            numNan  = np.sum(np.isnan(segments),axis=1)
            delInd0 = np.where(numNan>int((segments.shape[1]+0.5)/2))
            #Remove segments that are too long provided there are enough segments to sample
            if delInd0[0].shape[0]>1:
                numNan  = np.sum(~np.isnan(segments[delInd0]),axis=0)
                delInd1 = np.where(numNan>=2)
                segments = np.delete(segments,delInd1,1)
            segments = np.delete(segments,delInd0,0)
            xArr = np.indices((segments.shape))[0].flatten()
            segsFlat = segments.flatten()
            xArr = xArr[~np.isnan(segsFlat)]
            segsFlat = segsFlat[~np.isnan(segsFlat)]
            timePoints = np.arange(segments.shape[0])

            if segments.shape[1]<=2:
                fitOrAvgSegment = np.nanmean(segments,axis=1)
            else:
                transformed_x2 = dmatrix("cr(xArr, df=16)",
                                        {"xArr": xArr}, return_type='matrix')
                fit2 = sm.RLM(segsFlat, transformed_x2).fit()
                fitOrAvgSegment = fit2.predict(dmatrix("cr(timePoints, df=16)",
                                        {"timePoints": timePoints}, return_type='matrix'))
                fits[ch] = fit2
            self.goldenPulse[loc][ch] = fitOrAvgSegment
            '''if ch==0:
                plt.plot(y,label=('C1','C2','C3','C4'));plt.title('Location '+str(loc+1)); plt.legend(); plt.show()
            if ch==0 or ch==3:
                plt.plot(segments);plt.title('Location '+str(loc+1)+' Channel '+str(ch));plt.plot(fitOrAvgSegment,color='k'); plt.show()'''
        return fewPulseCount, noPulseCount, fits

    
    def FindNoisyPulseAndAverage(self, y, period, loc):
        '''
            Runs pulse detection and crops, then averages the pulses in each camera to get the exemplar pulse 
            for each camera

            Parameters
            ----------
            y : 2D array
                First dimension is the time the second is the number of channels/cameras

            period : float
                Period(in number of samples) of pulse detected in the signal

            loc : int
                Index of location where y(the signal being processed) was acquired
            
            Returns
            -------
            fewPulseCount: int
                Count where fewer than three pulses are detected in a channel

            noPulseCount: int
                Count where no pulse is detected in a channel
        '''
        noPulseCount  = 0
        fewPulseCount = 0
        m = np.mean(y, axis=0)
        ySmoothed = gaussian_filter1d(y-m, 1, axis=0, mode='reflect')
        yPrime    = ySmoothed[1:,:] - ySmoothed[:-1,:]
        minPeriod = int(np.floor(period-1.5))-int(period*0.1)
        chPeaks   = [find_peaks(yPrime[:,ch]**2, distance=minPeriod)[0] for ch in range(ySmoothed.shape[1])]
        
        chLengths = [peaks[1:]-peaks[:-1] for peaks in chPeaks]

        ### get segments only pick ones that are close < +/-2 of approximation
        chLengthsLst = [l.tolist() for l in chLengths]
        chLengthsLst = sum(chLengthsLst, [])
        maxLength = np.ceil(np.median(chLengthsLst)+1.5)
        minLength = np.floor(np.median(chLengthsLst)-1.5)

        L = int(maxLength+1)

        chUseSegment = [ ( (minLength<=np.array(lengths)).astype(int) + (maxLength>=np.array(lengths)).astype(int) ) == 2
                        for lengths in chLengths ]

        self.goldenPulse[loc] = [{}] * len(chUseSegment)

        for ch in range(len(chUseSegment)):
            #Exclude peaks that are too close to the start/end of the scan
            boundary = int(np.ceil(maxLength/2))
            if not len(chUseSegment[ch]):
                self.goldenPulse[loc][ch] = np.nan*np.zeros((10,))
                self.pulseDeteced[loc,ch] = False
                noPulseCount += 1
                continue
            if chPeaks[ch][0]<boundary:
                chUseSegment[ch][0] = False
            if chPeaks[ch][-2]>y.shape[0]-boundary: #Last peak has no period associated with it -- may need to rethink logic
                chUseSegment[ch][-1] = False
            
            nSegments =int(np.sum(chUseSegment[ch]))
            if nSegments == 0:
                self.goldenPulse[loc][ch] = np.nan*np.zeros((10,))
                self.pulseDeteced[loc,ch] = False
                noPulseCount += 1
                continue

            segments = np.zeros((L, nSegments))
            if nSegments<3:
                fewPulseCount += 1

            ##Crop valid segments and aggregate them for alignment
            ind=0
            for j in range(len(chLengths[ch])):
                if not chUseSegment[ch][j]:
                    continue
                t2 = chPeaks[ch][j] + boundary
                if t2>=y.shape[0]:
                    t2 = y.shape[0]-1
                t1 = t2-segments.shape[0]
                if t1<0:
                    t1 = 0
                    t2 = t1+segments.shape[0]
                segments[:, ind] = y[t1:t2, ch]
                ind += 1
            
            ### subpixel alignment
            moving = segments[2:-2, -1] #align everything to the last one
            ncc = np.zeros((5,))
            alignedsegments = np.zeros((len(moving), nSegments))
            alignedsegments[:, -1] = moving
            xx=np.arange(L)

            for j in range(nSegments-1):
                fixed = segments[:, j]
                for k in range(5):
                    out = pearsonr(fixed[k:k+len(moving)], moving)
                    ncc[k]=out[0]
                indmax = np.argmax(ncc)
                if indmax>0 and indmax<4:
                    a = 0.5*(ncc[indmax-1]-2*ncc[indmax]+ncc[indmax+1])
                    b = 0.5*(ncc[indmax+1]-ncc[indmax-1])
                    dk = -b/a/2
                else:
                    dk = 0
                kmax = indmax + dk
                f = interp1d(xx, fixed, kind='linear')
                xnew = np.arange(kmax, kmax+len(moving))
                if len(xnew)>alignedsegments.shape[0]:
                    xnew = xnew[:-1]
                alignedsegments[:, j] = f(xnew)
        
            ##Average aligned pulses to get golden pulse
            alignedsegments[alignedsegments == 0] = np.nan
            alignedsegments = np.nanmean(alignedsegments, axis=1)
            alignedsegments = alignedsegments[~np.isnan(alignedsegments)]
            self.goldenPulse[loc][ch] = alignedsegments

        self.AlignCamerasGoldenPulse(loc)
        return fewPulseCount, noPulseCount
    
    def AlignCamerasGoldenPulse(self,loc):
        '''
            Align the exemplar pulse between cameras

            Parameters
            ----------
            attribute : loc
                Location index to process
        '''

        m = 5 #how many steps to sample
        nc = len(self.goldenPulse[loc])
        ncc = np.zeros((nc, nc))

        continueToAlignment = False
        
        # find best camera
        for i in range(nc):
            for j in range(nc):
                if i==j:
                    continue #Otherwise we will be selecting the one with the best AUC
                bobs = len(self.goldenPulse[loc][i])==len(self.goldenPulse[loc][j])
                your = all(np.isfinite(self.goldenPulse[loc][i]))
                uncle = all(np.isfinite(self.goldenPulse[loc][j]))
                if bobs and your and uncle:
                    ncc[i, j] = pearsonr(self.goldenPulse[loc][i], self.goldenPulse[loc][j])[0]
                    continueToAlignment = True
        fixedind = np.argmax(np.nansum(ncc, axis=0))
        fixed = self.goldenPulse[loc][fixedind]

        if not continueToAlignment:
            return
        
        # align the other cameras to the best camera
        for i in range(nc):
            if i != fixedind and all(np.isfinite(self.goldenPulse[loc][i])):
                moving = self.goldenPulse[loc][i]
                ncc = np.zeros((2*m,))
                for j in range(-m, m):
                    ncc[j+m] = pearsonr(fixed, np.concatenate((moving[j:], moving[:j])))[0]
                ncc[np.isnan(ncc)]=0
                ind = np.argmax(ncc)-m
                self.goldenPulse[loc][i] = np.concatenate((moving[ind:], moving[:ind]))
    
    def GetAreaUnderCurve(self, waveform):
        '''
            Returns the area under the curve for the input waveform

            Parameters
            ----------
            waveform : 1D numpy array
                Input waveform
        '''
        mini = np.amin(waveform)
        waveform -= mini
        maxi = np.amax(waveform)
        waveform /= maxi
        return np.mean(waveform)

    def ComputeWaveformAttributes(self, goldenPulse):
        '''
            Computes the following features:
            Canopy - the proportion of time spent above 25% of the systolic-diastolic range
            Onset  - the time taken to reach systolic maximum from onset 0-90% ranges are used to filter out fit/averaging noise
            Relative onset - time spent in systolic part relative to whole pulse

            Parameters
            ----------
            goldenPulse : 1D numpy array
                Fitted/average pulses in the scan
        '''
        pulseRange = np.nanmax(goldenPulse)-np.nanmin(goldenPulse)
        totalTime  = np.count_nonzero(~np.isnan(goldenPulse))
        indsInCanopy = sum(goldenPulse < np.nanmax(goldenPulse)-pulseRange*0.25)
        canopy    = indsInCanopy/totalTime
        startInd  = np.nanargmax(goldenPulse[:5])
        sysRange  = np.nanmax(goldenPulse)-pulseRange*0.9
        sysInds   = goldenPulse>sysRange
        endInd    = np.nonzero(~sysInds)[0][0]
        onset     = float(abs(endInd-startInd))*self.scanParams.dt
        onsetProp = float(abs(endInd-startInd))/totalTime
        return canopy, onset, onsetProp, moment(goldenPulse,2)

    def ComputeWaveformAttributesFromFit(self, fit, segLen):
        '''
            Computes the following features by sampling the spline fit:
            Canopy - the proportion of time spent above 25% of the systolic-diastolic range
            Onset  - the time taken to reach systolic maximum from onset 0-90% ranges are used to filter out fit/averaging noise
            Relative onset - time spent in systolic part relative to whole pulse

            Parameters
            ----------
            fit : 1D natural spline fit parameters
            segLen : length of the segment
        '''
        sampling    = 0.1
        totalTime  = segLen
        timePoints  = np.arange(0,segLen-1+0.001,sampling)
        goldenPulse = fit.predict(dmatrix("cr(timePoints, df=16)",
                                                {"timePoints": timePoints}, return_type='matrix'))
        pulseRange = np.nanmax(goldenPulse)-np.nanmin(goldenPulse)
        indsInCanopy = sum(goldenPulse < np.nanmax(goldenPulse)-pulseRange*0.25)
        canopy    = indsInCanopy/totalTime*sampling
        startInd  = np.nanargmax(goldenPulse[:5])
        sysRange  = np.nanmax(goldenPulse)-pulseRange*0.9
        sysInds   = goldenPulse>sysRange
        endInd    = np.nonzero(~sysInds)[0][0]
        onset     = float(abs(endInd-startInd))*self.scanParams.dt*sampling
        onsetProp = float(abs(endInd-startInd))/totalTime*sampling
        if onsetProp>0.5: #Convenience for breakpoint
            asd = 0
        return canopy, onset, onsetProp, moment(goldenPulse,2)