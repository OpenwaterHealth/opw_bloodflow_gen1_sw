# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 11:03:25 2021

@author: soren
"""

#how many scans (subject, number per s, total) 
#how many complete (subject, position, all)
#how many good

import pickle
import batchfile
import speckle_fcns_v7 as fcns
import numpy as np
import matplotlib.pyplot as plt
import os

plt.close('all')

############################################
############### functions ##################
############################################

def getaverageddata(title, nsubjects, nrepeats, npositions, fpath, scannames):

    #initialize variables
    something = ['']*nsubjects
    something_diff = ['']*nsubjects

    for i in range(nsubjects):
    
        nrepeats = scanspersubject[i]
        something[i] = np.zeros((nrepeats, npositions))
        something_diff[i] = np.zeros((nrepeats, int(npositions/2)))

        #load data
        for j in range(nrepeats):
            loadname = fpath + '/' + scannames[i][j] + '.pkl'
            f = open(loadname, "rb")
            data = pickle.load(f)
            f.close()
            
            something[i][j, :] = data[0][title]
            something_diff[i][j, :] = data[0][title + '_diff']
            
    return something, something_diff

def getdata(title, nsubjects, nrepeats, ncameras, removebad, fpath, scannames):

    #initialize variables
    something = ['']*nsubjects
    something_diff = ['']*nsubjects
    a_all = np.zeros((nsubjects, nrepeats, 4))
    b_all = np.zeros((nsubjects, nrepeats, 4))

    for i in range(nsubjects):
    
        nrepeats = scanspersubject[i]
        something[i] = np.zeros((nrepeats, npositions, ncameras))
        something_diff[i] = np.zeros((nrepeats, ndiffs, ncameras))

        #load data
        for j in range(nrepeats):
            loadname = fpath + '/' + scannames[i][j] + '.pkl'
            f = open(loadname, "rb")
            data = pickle.load(f)
            f.close()

            # individual locations
            if ncameras>1:
                something[i][j, :, :] = fcns.combinepositions(data, title)
            else:
                something[i][j, :, :] = fcns.combinepositions(data, title)
            if removebad:
                mask = fcns.combinepositions(data, 'pulse')
                if ncameras == 1:
                    mask = np.mean(mask, axis=1)[:, np.newaxis]    
                mask[np.isfinite(mask)]=1
                something[i][j, :, :] *= mask

            # left right difference
            normdiff, diff, d = fcns.getdifference(data, title, col, row)
            something_diff[i][j, :, :] = diff
            if removebad:
                normdiff, mask, d = fcns.getdifference(data, 'pulse', col, row) 
                if ncameras == 1:
                    mask = np.mean(mask, axis=1)[:, np.newaxis]
                mask[np.isfinite(mask)]=1
                something_diff[i][j, :, :] *= mask
            
            a_all[i, j, :] = data[0]['a_all']
            b_all[i, j, :] = data[0]['b_all']
        
    return something, something_diff, a_all, b_all
    

def simplestats_hhc30(x, title, ylim0=[0, 0], ylim1=[0, 0]):
    
    #(subject)(repeats, positions, cameras)
    #subjects2use = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    gs = 3 #group size
    for i in range(nsubjects):
        if i==0:
            xx = x[i]
        else:
            xx = np.concatenate((xx, x[i]), axis=0)
    nscans = xx.shape[0]        
    npositions = xx.shape[1]
    if len(xx.shape)<3:
        xx = xx[:, :, np.newaxis]
    ncameras = xx.shape[2]
    d = [15, 22, 28, 35]

    #average over subjects and repeats
    m = np.nanmean(xx, axis=0)    
    s = np.nanstd(xx, axis=0)
    
    # plot with position on x axis
    nrows = max([ncameras, 2])
    fig, ax = plt.subplots(nrows = nrows, ncols=2)
    if npositions<4:
        fig.suptitle(title + ' Left Right Difference')
        xticks = [0, 1, 2]
        xticklabels = ['Forehead', 'Fissure', 'Temple']
    elif npositions<7:
        fig.suptitle(title)
        xticks = [0, 1, 2, 3, 4, 5]
        xticklabels = ['Forehead L', 'Fissure L', 'Temple L', 'Forehead R', 'Fissure R', 'Temple R']
    elif npositions<16:
        fig.suptitle(title + ' Left Right Difference')
        xticks = [1, 4, 7, 10, 13]
        xticklabels = ['Forehead V', 'Forehead H', 'Fissure D', 'Temple V', 'Temple H']
    else:
        fig.suptitle(title)
        xticks = [1, 4, 7, 10, 13, 16, 19, 22, 25, 28]
        xticklabels = ['Forehead V', 'Forehead H', 'Fissure D', 'Temple V', 'Temple H', 'Forehead V', 'Forehead H', 'Fissure D', 'Temple V', 'Temple H']
    t = np.arange(npositions)
    for i in range(ncameras):
        if npositions==30:
            ax[i, 0].plot([14.5, 14.5], [np.nanmin(xx[:, :, i]), np.nanmax(xx[:, :, i])], '--r')
            ax[i, 1].plot([14.5, 14.5], [np.amin(s), np.amax(s)], '--r')
        for k in range(int(npositions/gs)):
            for j in range(nscans):    
                ax[i, 0].plot(t[gs*k:gs*(k+1)], xx[j, gs*k:gs*(k+1), i], '.-', color=[0.7, 0.7, 0.7])
            ax[i, 0].plot(t[gs*k:gs*(k+1)], m[gs*k:gs*(k+1), i], 's-k')
            ax[i, 1].plot(t[gs*k:gs*(k+1)], s[gs*k:gs*(k+1), i], 's-k')
            if ylim0[1]>ylim0[0]:
                ax[i, 0].set_ylim(ylim0)
            if ylim1[1]>ylim1[0]:
                ax[i, 1].set_ylim(ylim1)    
        if ncameras>1:
            ax[i, 0].set_ylabel(str(d[i]) + ' mm')
        ax[i, 0].set_xticks([])
        ax[i, 1].set_xticks([])
    ax[ncameras-1, 0].set_xticks(xticks)
    ax[ncameras-1, 1].set_xticks(xticks)
    ax[ncameras-1, 0].set_xticklabels(xticklabels, rotation='vertical')
    ax[ncameras-1, 1].set_xticklabels(xticklabels, rotation='vertical')
    ax[0, 0].set_title('Average All Scans')
    ax[0, 1].set_title('Std All Scans')
    if nrows == ncameras:
        plt.subplots_adjust(bottom=0.2)
    else:
        ax[nrows-1, 0].axis('off')
        ax[nrows-1, 1].axis('off')
        
    plt.savefig(os.getcwd() + '/temp/' + title + str(npositions) + '.png')    
        
    return m, s 
    
############################################    
######## user input ########################
############################################

npositions = 30
ncameras = 4
ndiffs = 15
nrepeats = 3
scanorder = 'hhc30'
firstsubject = 1
legend = ['15 mm', '22 mm', '28 mm', '35 mm']
col = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
row = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
batchname = 'HHC_1B_Schrodinger_Healthy'
removebad = False

#savepath = 'C:/Users/soren/Documents/Python Scripts/bloodflow_proto2/temp'
######################
### how many scans ###
######################

calnames, scannames, datapath = eval('batchfile.' + batchname + '()')
datapath = datapath = os.path.expanduser(datapath)
fpath = datapath + '/Results'
nsubjects = len(scannames)
scanspersubject = [0]*nsubjects
for i in range(nsubjects):
    scanspersubject[i]=len(scannames[i])
totalscans = np.sum(scanspersubject)
print('Number of Subjects: ' + str(nsubjects))    
print('Scans per Subject: ' + str(scanspersubject))
print('Total Scans: ' + str(totalscans))

##############################
######### load data ##########
##############################

mean_mean, mean_diff, a_all, b_all = getdata('mean_mean', nsubjects, nrepeats, ncameras, removebad, fpath, scannames)
contrast_mean, contrast_diff, a_all, b_all = getdata('contrast_mean', nsubjects, nrepeats, ncameras, removebad, fpath, scannames)
amplitude, amplitude_diff, a_all, b_all = getdata('amplitude', nsubjects, nrepeats, ncameras, removebad, fpath, scannames)
modulationdepth, modulationdepth_diff, a_all, b_all = getdata('modulationdepth', nsubjects, nrepeats, ncameras, removebad, fpath, scannames)
auc, auc_diff, a_all, b_all = getdata('contrast_auc', nsubjects, nrepeats, ncameras, removebad, fpath, scannames)
pulse, pulse_diff, a_all, b_all = getdata('pulse', nsubjects, nrepeats, ncameras, removebad, fpath, scannames)

contrastratio, contrastratio_diff, a_all, b_all = getdata('contrast_mean_ratio', nsubjects, nrepeats, 1, removebad, fpath, scannames)
amplituderatio, amplituderatio_diff, a_all, b_all = getdata('amplitude_ratio', nsubjects, nrepeats, 1, removebad, fpath, scannames)
modulationdepthratio, modulationdepthratio_diff, a_all, b_all = getdata('modulationdepth_ratio', nsubjects, nrepeats, 1, removebad, fpath, scannames)
aucratio, aucratio_diff, a_all, b_all = getdata('contrast_auc_ratio', nsubjects, nrepeats, 1, removebad, fpath, scannames)

'''contrast_mean_combined, contrast_mean_combined_diff = getaverageddata('contrast_mean_combined', nsubjects, nrepeats, 6, fpath, scannames)
amplitude_combined, amplitude_combined_diff = getaverageddata('amplitude_combined', nsubjects, nrepeats, 6, fpath, scannames)
modulationdepth_combined, modulationdepth_combined_diff = getaverageddata('modulationdepth_combined', nsubjects, nrepeats, 6, fpath, scannames)
contrast_auc_combined, contrast_auc_combined_diff = getaverageddata('contrast_auc_combined', nsubjects, nrepeats, 6, fpath, scannames)
'''
mua, mua_diff, a_all, b_all = getdata('mua_mean', nsubjects, nrepeats, 1, removebad, fpath, scannames)

################################################
############## display results #################
################################################

simplestats_hhc30(mean_mean, 'Mean')
simplestats_hhc30(mean_diff, 'Mean')
simplestats_hhc30(contrast_mean, 'Contrast', [0.2, 0.4], [0, 0.05])
simplestats_hhc30(contrast_diff, 'Contrast', [-0.1, 0.1], [0, 0.03])
simplestats_hhc30(amplitude, 'Amplitude', [0.006, 0.06], [0.004, 0.01])
simplestats_hhc30(amplitude_diff, 'Amplitude')
simplestats_hhc30(modulationdepth, 'Modulation Depth')
simplestats_hhc30(modulationdepth_diff, 'Modulation Depth')
simplestats_hhc30(auc, 'AUC')
simplestats_hhc30(auc_diff, 'AUC')
plt.close('all')
simplestats_hhc30(contrastratio, 'Contrast Ratio', [0.5, 1], [0, 0.1])
simplestats_hhc30(contrastratio_diff, 'Contrast Ratio', [-0.1, 0.1], [0, 0.1])
simplestats_hhc30(amplituderatio, 'Amplitude Ratio')
simplestats_hhc30(amplituderatio_diff, 'Amplitude Ratio')
simplestats_hhc30(modulationdepthratio, 'Modulation Depth Ratio')
simplestats_hhc30(modulationdepthratio_diff, 'Modulation Depth Ratio')
simplestats_hhc30(aucratio, 'AUC Ratio')
simplestats_hhc30(aucratio_diff, 'AUC Ratio')
plt.close('all')
'''simplestats_hhc30(contrast_mean_combined, 'Combined Contrast')
simplestats_hhc30(contrast_mean_combined_diff, 'Combined Contrast')
simplestats_hhc30(amplitude_combined, 'Combined Amplitude')
simplestats_hhc30(amplitude_combined_diff, 'Combined Amplitude')
simplestats_hhc30(modulationdepth_combined, 'Combined Modulation Depth')
simplestats_hhc30(modulationdepth_combined_diff, 'Combined Modulation Depth')
simplestats_hhc30(contrast_auc_combined, 'Combined AUC')
simplestats_hhc30(contrast_auc_combined_diff, 'Combined AUC')
'''
plt.close('all')
simplestats_hhc30(mua, 'Absorption (~CBV)')
simplestats_hhc30(mua_diff, 'Absorption (~CBV)')
plt.close('all')

####### good scans #######
print('#############################')     
print(' Fourier test ') 
allgood = np.zeros((nsubjects, nrepeats, npositions, ncameras))
for i in range(nsubjects):
    allgood[i, :, :, :] = pulse[i]
allgood[np.isfinite(allgood)]=1
good_subjectXcamera = np.nansum(np.nansum(allgood, axis=1), axis=1)
pgood_subjectXcamera = 100*good_subjectXcamera/nrepeats/npositions 
good_positionXcamera = np.nansum(np.nansum(allgood, axis=0), axis=0)
pgood_positionXcamera = 100*good_positionXcamera/nsubjects/nrepeats
c='rgbk'
print('########## good_subjectXcamera #########')
print(good_subjectXcamera.T)
print('########## good_positionXcamera #########')
print(good_positionXcamera.T)
fig, ax = plt.subplots(nrows=1, ncols=2)
fig.suptitle('Percent of Acquisitions Where Pulse is Dominant Feature')
for j in range(ncameras):
    ax[0].plot(np.arange(firstsubject, firstsubject+nsubjects), pgood_subjectXcamera[:, j], 'o-', color=c[j])
ax[0].legend(['15 mm', '22 mm', '28 mm', '35 mm'])
ax[0].set_xlabel('Subject Number')
ax[0].set_ylabel('Percent')
for i in range(10):
    for j in range(ncameras):
        ax[1].plot(np.arange(3*i, 3*i+3), pgood_positionXcamera[3*i:3*i+3, j], 'o-', color=c[j])   
ax[1].plot([14.5, 14.5], [75, 100], '--', color=[0.5, 0.5, 0.5])
ax[1].legend(['15 mm', '22 mm', '28 mm', '35 mm'])
xticks = [1, 4, 7, 10, 13, 16, 19, 22, 25, 28]
xticklabels = ['Forehead V', 'Forehead H', 'Fissure D', 'Temple V', 'Temple H', 'Forehead V', 'Forehead H', 'Fissure D', 'Temple V', 'Temple H']
ax[1].set_xticks(xticks)
ax[1].set_xticklabels(xticklabels, rotation='vertical')
ax[1].set_ylabel('Percent')

####### fitting coefficients #######
fig, ax = plt.subplots(nrows=1, ncols=2)
ax[0].plot(np.reshape(a_all, (nsubjects*nrepeats, ncameras)),'.-')
ax[0].legend(legend)
ax[0].set_title('a')
ax[1].plot(np.reshape(b_all, (nsubjects*nrepeats, ncameras)),'.-')
ax[1].legend(legend)
ax[1].set_title('b')




