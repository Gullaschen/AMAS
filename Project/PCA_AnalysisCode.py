import sys
import numpy as np
import re
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import os
import matplotlib
from os import listdir
from os.path import isfile, join
from matplotlib import cm
from matplotlib.cm import ScalarMappable
from scipy import interpolate
print(sys.executable)
print(sys.path)
sys.path.insert(0,'C:\\Users\\Villads\\anaconda3\\pkgs')
import sklearn
from sklearn.decomposition import PCA

matplotlib.rcParams.update({'font.size': 18})  
#%%
def read_data(file_name, num_columns=2, sep=','):
    i = 0
    data = [[] for _ in range(num_columns)]

    with open(os.path.join(file_name), mode='r') as csvfile:
        for row in csvfile:  
            if i > 0:
                row = row.split(sep)
                for j in range(num_columns):
                    data[j].append(float(row[j]))
            i += 1
    return np.array(data)
def line_function(x1,y1,x2,y2, x):
    
    # Calculate the slope and intercept of the line
    slope = (y2 - y1) / (x2 - x1)
    intercept = y1 - slope * x1
    
    # Calculate the y-values of the line for the given x-values
    y = slope * x + intercept
    return y
def read_powerfile(file_name, num_columns=2, sep=','):
    i = 0
    data = [[] for _ in range(num_columns)]

    with open(os.path.join(file_name), mode='r') as csvfile:
        for row in csvfile:  
            if i > 0 and i%2 ==0:
                row = row.split(sep)
                for j in range(num_columns):
                    data[j].append(float(row[j]))
            i += 1
    return np.array(data)

def extract_number(filename):
    # Extract the numeric part after "__ExcWL" using regular expression
    match = re.search(r'_(\d+)(?=__ExcWL)', filename)
    return int(match.group(1)) if match else float('inf')


def read_folder_ofSpectra_andPowerCorrect(path, num_columns=2, sep=','):
    files = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and "_PowerCalibrationFile_" not in f]
    
    #Sorting
    files = sorted(files, key=extract_number)

    all_data = []
    power_corr = read_powerfile(os.path.join(path, '_PowerCalibrationFile_.csv'), num_columns=2, sep=',')

    for file_name,ycorr in zip(files,power_corr[1]):
        print(file_name)
        if not file_name == '_PowerCalibrationFile_':
            data = read_data(file_name, num_columns, sep)
  #          plt.plot(data[0],data[num_columns-1]/ycorr)
            
            y_corr = y_corr_func(data[0])
            if correct:     
                all_data.append(data[num_columns-1]/ycorr / y_corr)
            if not correct:
                all_data.append(data[num_columns-1]/ycorr)
            
  #  plt.title('data_loaded')
  #  plt.show()
    return data[0],np.array(all_data)

def smoothing(x,data):
    y = data  + max(data)
   # plt.plot(x,y)
    
    
    percentage_diff = np.abs(np.diff(y))
    rolling_avg = np.convolve(y, np.ones(5)/5, mode='valid',)
    diff_from_avg = (y[4:] - rolling_avg) /rolling_avg
    
    
    outliers_mask = diff_from_avg > 0.1
    # Replace outliers with the average of neighboring points
    
    
    outliers = np.where(outliers_mask)[0] + 4
    n = 0
    while len(outliers) > 0:
        n +=1
        for i in outliers:
            
            if i < 15 or i > len(y)-15 :

                continue
            y[i] = (y[i -10 ] + y[i+10]) / 2
            
     #   plt.plot(x,y)
      #  plt.scatter(x[outliers],y[outliers],c='r')
        np.where(-1 == y)
        
        percentage_diff = np.abs(np.diff(y))
        rolling_avg = np.convolve(y, np.ones(5)/5, mode='valid',)
        diff_from_avg = (y[4:] - rolling_avg) /rolling_avg
    
    
        outliers_mask = diff_from_avg > 0.3
        # Replace outliers with the average of neighboring points
    
        outliers = np.where(outliers_mask)[0] + 4
        if n > 10:
            return data
    #plt.show()
    return y - max(data)

def savefile(filepath,x,y):
    with open(filepath, 'w') as file:
        np.savetxt(file, np.column_stack((x,y)), delimiter=',')
        file.write('\n')
#%%
# =============================================================================
# =============================================================================
# # Correction of 880 nm band, the 1060 dont need;
# =============================================================================
# =============================================================================
x_corr,y_corr = read_data(r'C:/Users/zct154/Documents/SuCo/750B-Quantum Effieciency Suco.csv')
def smooth_data(numbers, window_size=10):
    smoothed_numbers = []

    for i in range(len(numbers)):
        start_index = max(0, i - window_size // 2)
        end_index = min(len(numbers), i + window_size // 2 + 1)

        # Extract the neighboring elements
        neighbors = numbers[start_index:end_index]

        # Calculate the average and append to the smoothed list
        smoothed_numbers.append(sum(neighbors) / len(neighbors))

    return smoothed_numbers
y_corr = smooth_data(y_corr, window_size=10)
y_corr_func =  interpolate.interp1d(x_corr,y_corr,fill_value=0,bounds_error=False)


#%%
# =============================================================================
# 300 K
# =============================================================================
#path to folder with data:
correct = True
path = r'C:\Users\zct154\Documents\SuCo\77K\EuDPA_45_77K_HQ_EmissionFilesForExcitationSpectra'
x, data = read_folder_ofSpectra_andPowerCorrect(path, num_columns=2, sep=',')
plt.plot(x, data[40])
#If background data provided:
Background_data = False
#bkgpath = r'C:\Users\Villads\Documents\Project_NdpDO3A_Charlie\Optical data\DMSO\77K\BKGNdPDO3A_DMSO_10s_77k_c1070_EmissionFilesForExcitationSpectra'
#x, bkgdata = read_folder_ofSpectra_andPowerCorrect(bkgpath, num_columns=2, sep=',')

#Specify wavelength region 
wl = np.linspace(458,475,len(data))

#What is the range for Background subtraction:
condition1 = ((560 < x) & (x < 572)) 
condition2 = ((733 < x) & (x < 737))

#Area to integrate of the emission peak
condition = [i for i, value in enumerate(x) if 575 < value and value < 725]

#Wavelengths ranges on the excitation side:
Wavelengths_ranges = [
    (wl > 462) & (473 > wl) ,
    ]

Focus_band =0

name = '300K_c880'
#%%
# =============================================================================
# Data cleaning
# =============================================================================
#Crappy backround:

data_corrected = []

for i,dat in enumerate(data):     
    #determining bkg funtion
    x1,x2 = np.mean(x[condition1]),np.mean(x[condition2])
    y1,y2 = np.mean(dat[condition1]),np.mean(dat[condition2])
    
    bkgy = line_function(x1,y1,x2,y2, x)
    y_bkg_subtracted = dat - bkgy
    if Background_data:
        y1,y2 = np.mean(bkgdata[i][condition1]),np.mean(bkgdata[i][condition2])
        bkgy = line_function(x1,y1,x2,y2, x)
        bkgy_bkg_subtracted = bkgdata[i] - bkgy
        y_bkg_subtracted = y_bkg_subtracted - bkgy_bkg_subtracted

    #subtracting bkg:    
    y_bkg_subtracted = smoothing(x,y_bkg_subtracted)
    
    x1,x2 = np.mean(x[condition1]),np.mean(x[condition2])
    y1,y2 = np.mean(y_bkg_subtracted[condition1]),np.mean(y_bkg_subtracted[condition2])
    bkgy = line_function(x1,y1,x2,y2, x)
    y_corr = y_bkg_subtracted - bkgy
    plt.plot(x,y_corr)

    data_corrected.append(smoothing(x,y_corr))
data_corrected = np.array(data_corrected)
#%%
# =============================================================================
# Integrate for area and excitation plots from individual wavelengths in emission.
# ====================================================================

%matplotlib inline
area = []
area_tot = []
#condition_wl = [i for i, value in enumerate(wl) if 570 < value and value < 605]

for y in data_corrected:
    area.append(np.array(y[condition]).T)
    #Integrate peak with respect to condition:
    area_tot.append(-np.trapz(np.array(y[condition]),10**7/(np.array(x[condition]))))
    
    
area_tot = np.array(area_tot)
plt.plot(x,sum(data_corrected)/max(sum(data_corrected)),c = 'k',lw = 2)


#%%
# =============================================================================
# Excitation Spectrum:
# =============================================================================
fig, ax = plt.subplots(figsize=[8,6])
ax.plot(wl[Wavelengths_ranges[Focus_band]],area_tot[Wavelengths_ranges[Focus_band]],c = 'k',lw = 2)

len(Wavelengths_ranges[Focus_band])


ax.plot(wl,area_tot,c = 'k',lw = 2)
ax.set_xlabel('wavelength (nm)')
ax.set_ylabel('Intensity')
#ax.set_xlim(572,590)


#filepath = r'C:\Users\Villads\Documents\Project_NdpDO3A_Charlie\Optical data\DMSO\Generated data\300 K\Excitation{}.csv'.format(name)
#savefile(filepath,wl,area_tot)
#%%
# =============================================================================
# summed emission from each excitation band
# =============================================================================
Averaged_emission = []
for wl_r in Wavelengths_ranges:
    Averaged_emission.append(np.mean(data_corrected[wl_r],axis=0))
    

fig, ax = plt.subplots(figsize=[8,6])

plt.plot(x,Averaged_emission[0] / max(Averaged_emission[0]),label = 'band average {} nm'.format(round(np.mean(wl[Wavelengths_ranges[0]]))))

plt.legend()

ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Normalized Intensity (a.u.)')
plt.title('Averaged Emission')
#filepath = r'C:\Users\Villads\Documents\Project_NdpDO3A_Charlie\Optical data\DMSO\Generated data\300 K\emission_average{}.csv'.format(name)
#savefile(filepath,x,Averaged_emission[1] / max(Averaged_emission[1]))
#%%
# =============================================================================
# Individual emission spectra Plotted
# =============================================================================
fig, ax = plt.subplots(figsize=[8,6])
cmap = plt.get_cmap('viridis')
emission_data = data_corrected[Wavelengths_ranges[Focus_band]]
colors = [cmap(i / len(emission_data)) for i in range(len(emission_data))]
lines = []
for i,y in enumerate(emission_data):
    ax.plot(x, y, color=colors[i])

plt.plot(x,Averaged_emission[Focus_band],label = 'Averaged', c = 'k', lw = 2.5)
plt.legend()

high = int(wl[Wavelengths_ranges[Focus_band]][-1])
low = int(wl[Wavelengths_ranges[Focus_band]][0])
sm = ScalarMappable(cmap=cmap, norm=plt.Normalize(low,high))
sm.set_array([])

cbar = plt.colorbar(sm, ticks=range(low,high, 5))
# Add the color bar legend
cbar.set_label('Excitation Wavelength (nm)')
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Intensity (a.u.)')
#%%
fig, ax = plt.subplots(figsize=[8,6])
cmap = plt.get_cmap('viridis')
emission_data = data_corrected[Wavelengths_ranges[Focus_band]]
colors = [cmap(i / len(emission_data)) for i in range(len(emission_data))]
lines = []
for i,y in enumerate(emission_data):
    ax.plot(x, y/max(y), color=colors[i])

plt.plot(x,Averaged_emission[Focus_band] / max(Averaged_emission[Focus_band]),label = 'Averaged', c = 'k', lw = 2.5)
plt.legend()

high = int(wl[Wavelengths_ranges[Focus_band]][-1])
low = int(wl[Wavelengths_ranges[Focus_band]][0])
sm = ScalarMappable(cmap=cmap, norm=plt.Normalize(low,high))
sm.set_array([])

cbar = plt.colorbar(sm, ticks=range(low,high, 5))
# Add the color bar legend
cbar.set_label('Excitation Wavelength (nm)')
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Normalized Intensity (a.u.)')

#%%
# =============================================================================
# Prepare for PCA
# =============================================================================
fig, ax = plt.subplots(figsize=[8,6])

wl_r = Wavelengths_ranges[Focus_band]
ax.plot(wl[wl_r],area_tot[wl_r],c = 'k')
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Intensity (a.u.)')
plt.title('Excitation range for PCA')
#%%
Y = []
x_pca =x[condition]
fig, ax = plt.subplots(figsize=[8,6])
for y in data_corrected[wl_r]:
    y = y[condition]
    new_y = [(a - np.mean(y))/np.std(y) for a in y ]    
  #  new_y = y/max(y)
    plt.plot(x_pca,new_y)
    Y.append(new_y)

ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Intensity (a.u.)')
plt.title('input data for PCA')

#%%
# =============================================================================
# PCA Analysis
# =============================================================================


pca = PCA(n_components=3)
    
X_proj = pca.fit_transform(np.array(Y)) 
np.shape(X_proj)

cmap = plt.get_cmap('plasma')
colors = [cmap(i) for i in np.linspace(0, 1, len(X_proj))]
colors = [cmap(i) for i in (X_proj[:,1]-min(X_proj[:,1]))/max(X_proj[:,1]-min(X_proj[:,1]))]

#colors = [cmap(i) for i in np.linspace(0,1, len(yheat))]


pc_a,pc_b =0,1
# Create the scatter plot
fig, ax = plt.subplots(figsize=[12, 6])
ax.scatter(X_proj[:, pc_a], X_proj[:, pc_b], color=colors, alpha=0.8)

# Set labels for the axes
ax.set_xlabel('PC{} ({}%)'.format(pc_a + 1, round(pca.explained_variance_ratio_[pc_a], 3)))
ax.set_ylabel('PC{} ({}%)'.format(pc_b + 1, round(pca.explained_variance_ratio_[pc_b], 3)))

# Define the font size for the labels
label_fontsize = 12  # Adjust the font size as needed

# Add integer labels for each data point with a specified font size
for i, (a, b) in enumerate(zip(X_proj[:, pc_a], X_proj[:, pc_b])):
    ax.annotate(str(i), (a, b), textcoords="offset points", xytext=(0, 10), ha='center', fontsize=label_fontsize)
plt.title('PCA model')
plt.show()

individual_components = pca.components_ 
   
fig, ax = plt.subplots(figsize=[8,6])
plt.plot(x_pca,individual_components[0], label = 'PC1')
plt.plot(x_pca,individual_components[1], label = 'PC2')
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Intensity (a.u.)')
plt.legend()
plt.title('PCA components')
#%%
extreme_point_positive = np.zeros_like(X_proj[0])
extreme_point_negative = np.zeros_like(X_proj[0])

# Choose a multiplier for the extremes (you can adjust this)
extreme_multiplier =1
# Set the extreme values along the principal components
extreme_point_positive[0] = pca.explained_variance_[0] * extreme_multiplier
extreme_point_negative[0] = pca.explained_variance_[0] * -extreme_multiplier
extreme_point_positive[1] = pca.explained_variance_[1] * extreme_multiplier
extreme_point_negative[1] = pca.explained_variance_[1] * -extreme_multiplier
#extreme_point_positive[2] = pca.explained_variance_[2] * extreme_multiplier
#extreme_point_negative[2] = pca.explained_variance_[2] * -extreme_multiplier



# Use inverse_transform to map the extreme points back to the original feature space
extreme_point_positive_original = pca.inverse_transform([extreme_point_positive])
extreme_point_negative_original = pca.inverse_transform([extreme_point_negative])

# Print or use the extreme points in the original feature space
print("Extreme Point Positive (Original):", extreme_point_positive_original)
print("Extreme Point Negative (Original):", extreme_point_negative_original)

plt.figure(figsize=(9, 6))
plt.plot(x_pca,extreme_point_negative_original[0]- min(extreme_point_negative_original[0]) ,label = 'Spectra A')
plt.plot(x_pca,extreme_point_positive_original[0]- min(extreme_point_positive_original[0]) ,label = 'Spectra B')

plt.legend()
plt.xlabel('wavelengths (nm)')
plt.ylabel('Intensity (a.u.)')
plt.title('Reconstructed components by PCA')
#%%
# =============================================================================
# Reconstruction of spectra, Manual fingerprinting:
# =============================================================================
X_reconstructed = pca.inverse_transform(X_proj)
plt.figure(figsize=(9, 6))
for i in [0, 19]:  # Plot the first 5 spectra for illustration
    plt.plot(x_pca,Y[i], label='Original {}'.format(i), linewidth=2)
    plt.plot(x_pca,X_reconstructed[i], label='reconstructed = {}'.format(i))
    
    plt.legend()
   # plt.show()
plt.title('Spectral Reconstruction of 9 and 48')
plt.xlabel('Wavelengths (nm)')
plt.ylabel('Intensity (a.u.)')
ax.set_xlim(840,930)

#%%
# =============================================================================
# Generation of fingerprint structures1
# =============================================================================
#Specify two spectra that well represents the most different ones. Then we subtract them from each other.
Spectra1 = 15
Spectra2 = 18
#Specify the scale used to subtract spectra 2 from 1 with
scale2 = 0.3

bkgrange1, bkgrange2 = ((840 < x_pca) & (x_pca < 855)) ,((925 < x_pca) & (x_pca < 930)) 

#bkg subtracts the reconstructed spectra, as they are generated with the mean set to 0, not the mininum.
bkg1 = line_function(np.mean(x_pca[bkgrange1]),np.mean(X_reconstructed[Spectra1][bkgrange1]),np.mean(x_pca[bkgrange2]),np.mean(X_reconstructed[Spectra1][bkgrange2]), x_pca)
recon_1 = X_reconstructed[Spectra1] - bkg1
bkg2 = line_function(np.mean(x_pca[bkgrange1]),np.mean(X_reconstructed[Spectra2][bkgrange1]),np.mean(x_pca[bkgrange2]),np.mean(X_reconstructed[Spectra2][bkgrange2]), x_pca)
recon_2 = (X_reconstructed[Spectra2] - bkg2) *scale2

#Plots the reconstructucted spectra along with the first species we are manually generating
%matplotlib inline

plt.figure(figsize=(9, 6))
plt.plot(x_pca,recon_1, label='recon1', linewidth=2)
plt.plot(x_pca,recon_2, label='scaled recon2', linewidth=2)
#%
species1 = recon_1 - recon_2
plt.plot(x_pca,species1 , label='species1 = recon1 - scaled recon2', linewidth=2)

plt.xlabel('Wavelengths (nm)')
plt.ylabel('Normalized Intensity (a.u.)')
plt.title('Manual construction of species1')
plt.legend()
#%%
# =============================================================================
# Generation of fingerprint structures
# =============================================================================
#Specify the scale used to subtract species 1 from spectra 2 with
scale2 = 0.2

#Plots the reconstructucted spectra along with the second species we are manually generating
plt.figure(figsize=(9, 6))
plt.plot(x_pca,species1 * scale2 , label='scaled species1', linewidth=2)
plt.plot(x_pca,recon_2 , label='recon2', linewidth=2)

species2 = recon_2*1 - species1 * scale2

plt.plot(x_pca,species2 , label='species2 = recon2 - scaled species1', linewidth=2)

plt.xlabel('Wavelengths (nm)')
plt.ylabel('Normalized Intensity (a.u.)')
plt.title('Manual construction of species2')
plt.legend()
#%%
plt.figure(figsize=(9, 6))
plt.plot(x_pca,species1/max(species1),label = 'Spectra A')
plt.plot(x_pca,species2/max(species2),label = 'Spectra B')



plt.legend()
plt.xlabel('Wavelengths (nm)')
plt.ylabel('Intensity (a.u.)')
plt.title('Manual Reconstructed Components')
#%%
filepath = r'C:\Users\Villads\Documents\Project_NdpDO3A_Charlie\Optical data\DMSO\Generated data\300 K\Spectra_A_{}.csv'.format(name)
savefile(filepath,x[condition],extreme_point_negative_original[0] - min(extreme_point_negative_original[0] ))
filepath = r'C:\Users\Villads\Documents\Project_NdpDO3A_Charlie\Optical data\DMSO\Generated data\300 K\Spectra_B_{}.csv'.format(name)
savefile(filepath,x[condition],extreme_point_positive_original[0] - min(extreme_point_positive_original[0] ))
#%%
def Fit_LinearComp(x,y,func1,func2):     
    Guesses = np.array([0.6,0.2,0])

    def Loss(arguments):
        a= arguments[0]
        b= arguments[1]
        c =arguments[2]
        return np.sum((func1 * a + func2 * b + c - y )**2)
    

    # Define constraints to ensure a and b are non-negative
    bounds = [(0, 10000000000)] + [(0, 10000000000)]  + [(0, 10000000000)]

    fit_values = minimize(Loss, Guesses, method='Powel''l', bounds=bounds).x
    a,b,c = fit_values[0],fit_values[1],fit_values[2]
    return a,b,c

A,B,CRAP = [],[],[]

func1 = extreme_point_positive_original[0] - min(extreme_point_positive_original[0] )
func2 = extreme_point_negative_original[0]- min(extreme_point_negative_original[0] )




x_fit = x[condition]
total_area = []
error =[]
for i in range(len(data_corrected)):

    y = data_corrected[i][condition]
    
    len(func1)
    
    a,b,c = Fit_LinearComp(x_fit,y,func1,func2)

# =============================================================================
#     A.append(a)
#     B.append(b)
# =============================================================================
    
    componentA = func1 * a
    A.append(-np.trapz(np.array(componentA),10**7/(np.array(x_fit))))
    componentB = func2 * b
    B.append(-np.trapz(np.array(componentB),10**7/(np.array(x_fit))))
    
    total_area.append(-np.trapz(np.array(y),10**7/(np.array(x_fit))))
    error.append(np.sum((func1 * a + func2 * b + c - y )**2  ))

    
    if i > 0:
        plt.plot(x_fit,func1 * a + func2 * b +c ,label = 'fit')
        plt.plot(x_fit,y,label = 'data')
        plt.title('Excitaiton Wavelength = {}'.format(wl[i]))
        plt.xlabel('wavelengths (nm)')
        plt.ylabel('Intensity (a.u.)')
        plt.legend()
        plt.show()

A= np.array(A)
B= np.array(B)
sum(A)/sum(B)

#%%
fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 8))
ax.plot(wl,A, alpha = 0.5,label = 'component A')
ax.plot(wl,B, alpha = 0.5,label = 'component B')
ax.plot(wl,A+B,label = 'fit from PCA',c = 'r')

ax.plot(wl,total_area,c = 'k',label = 'measured excitation')
ax.legend()
ax2.plot(wl,np.array(np.sqrt(error)), label='data - fit_y')

ax2.legend()

ax2.set_xlabel('wavelengths (nm)')
ax.set_ylabel('Intensity (a.u.)')
ax2.set_ylabel('Error (a.u.)')

#%%
filepath = r'C:\Users\Villads\Documents\Project_NdpDO3A_Charlie\Optical data\DMSO\Generated data\300 K\Excitation_speciesA_{}.csv'.format(name)
savefile(filepath,wl,A)


filepath = r'C:\Users\Villads\Documents\Project_NdpDO3A_Charlie\Optical data\DMSO\Generated data\300 K\Excitation_speciesB_{}.csv'.format(name)
savefile(filepath,wl,B)

filepath = r'C:\Users\Villads\Documents\Project_NdpDO3A_Charlie\Optical data\DMSO\Generated data\300 K\Excitation_Error_{}.csv'.format(name)
savefile(filepath,wl,np.array(np.sqrt(error)))

#%%
# =============================================================================
# =============================================================================
# # Second range fit:
# =============================================================================
# =============================================================================
path = r'C:\Users\Villads\Documents\Project_NdpDO3A_Charlie\Optical data\DMSO\300K\1070em\NdPDO3A_580_DMSO_300K_10s_1070_EmissionFilesForExcitationSpectra'
x, data = read_folder_ofSpectra_andPowerCorrect(path, num_columns=2, sep=',')
wl2 = np.linspace(566,607,len(data))

data_corrected2 = []


for i,dat in enumerate(data):     
    #determining bkg funtion
    x1,x2 = np.mean(x[condition1]),np.mean(x[condition2])
    y1,y2 = np.mean(dat[condition1]),np.mean(dat[condition2])
    
    bkgy = line_function(x1,y1,x2,y2, x)
    y_bkg_subtracted = dat - bkgy
    if Background_data:
        y1,y2 = np.mean(bkgdata[i][condition1]),np.mean(bkgdata[i][condition2])
        bkgy = line_function(x1,y1,x2,y2, x)
        bkgy_bkg_subtracted = bkgdata[i] - bkgy
        y_bkg_subtracted = y_bkg_subtracted - bkgy_bkg_subtracted


    #subtracting bkg:    
    y_bkg_subtracted = smoothing(x,y_bkg_subtracted)
    
    x1,x2 = np.mean(x[condition1]),np.mean(x[condition2])
    y1,y2 = np.mean(y_bkg_subtracted[condition1]),np.mean(y_bkg_subtracted[condition2])
    bkgy = line_function(x1,y1,x2,y2, x)
    y_corr = y_bkg_subtracted - bkgy
   # plt.plot(x,y_corr)

    data_corrected2.append(smoothing(x,y_corr))
data_corrected2 = np.array(data_corrected2)
#%%

def Fit_LinearComp(x,y,func1,func2):     
    Guesses = np.array([0.6,0.2,0])

    def Loss(arguments):
        a= arguments[0]
        b= arguments[1]
        c =arguments[2]
        return np.sum((func1 * a + func2 * b + c - y )**2)
    

    # Define constraints to ensure a and b are non-negative
    bounds = [(0, 10000000000)] + [(0, 10000000000)]  + [(0, 10000000000)]

    fit_values = minimize(Loss, Guesses, method='Powel''l', bounds=bounds).x
    a,b,c = fit_values[0],fit_values[1],fit_values[2]
    return a,b,c

A,B,CRAP = [],[],[]

func1 = extreme_point_positive_original[0] - min(extreme_point_positive_original[0] )
func2 = extreme_point_negative_original[0]- min(extreme_point_negative_original[0] )




x_fit = x[condition]
total_area = []
error =[]
for i in range(len(data_corrected2)):

    y = data_corrected2[i][condition]
    
    len(func1)
    
    a,b,c = Fit_LinearComp(x_fit,y,func1,func2)

# =============================================================================
#     A.append(a)
#     B.append(b)
# =============================================================================
    
    componentA = func1 * a
    A.append(-np.trapz(np.array(componentA),10**7/(np.array(x_fit))))
    componentB = func2 * b
    B.append(-np.trapz(np.array(componentB),10**7/(np.array(x_fit))))
    
    total_area.append(-np.trapz(np.array(y),10**7/(np.array(x_fit))))
    error.append(np.sum((func1 * a + func2 * b + c - y )**2  ))

    
# =============================================================================
#     if (i  < 220 and i > 210):
#         plt.plot(x_fit,func1 * a + func2 * b +c ,label = 'fit')
#         plt.plot(x_fit,y,label = 'data')
#         plt.title('Excitaiton Wavelength = {}'.format(wl[i]))
#         plt.xlabel('wavelengths (nm)')
#         plt.ylabel('Intensity (a.u.)')
#         plt.legend()
#         plt.show()
# =============================================================================

A= np.array(A)
B= np.array(B)
#%%
fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 8))
ax.plot(wl2,A, alpha = 0.5,label = 'component A')
ax.plot(wl2,B, alpha = 0.5,label = 'component B')
ax.plot(wl2,A+B,label = 'fit from PCA',c = 'r')

ax.plot(wl2,total_area,c = 'k',label = 'measured excitation')
ax.legend()
ax2.plot(wl2,np.array(np.sqrt(error)), label='data - fit_y')

ax2.legend()

ax2.set_xlabel('wavelengths (nm)')
ax.set_ylabel('Intensity (a.u.)')
ax2.set_ylabel('Error (a.u.)')


#%%
filepath = r'C:\Users\Villads\Documents\Project_NdpDO3A_Charlie\Optical data\DMSO\Generated data\300 K\Excitation_speciesA_{}2.csv'.format(name)
savefile(filepath,wl,A)


filepath = r'C:\Users\Villads\Documents\Project_NdpDO3A_Charlie\Optical data\DMSO\Generated data\300 K\Excitation_speciesB_{}2.csv'.format(name)
savefile(filepath,wl,B)

filepath = r'C:\Users\Villads\Documents\Project_NdpDO3A_Charlie\Optical data\DMSO\Generated data\300 K\Excitation_Error_{}2.csv'.format(name)
savefile(filepath,wl,np.array(np.sqrt(error)))

