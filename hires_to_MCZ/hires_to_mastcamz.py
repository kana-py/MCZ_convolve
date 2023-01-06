#PROGRAM NAME: hires_to_mastcam_wbayer.pro
#WRITTEN BY: Ryan Anderson, Melissa Rice, Kathleen Hoza
#PURPOSE: to read in high resolution spectra and 
#convolve them to Mastcam bandpasses
# 
# Adapted to Mastcam from Melissa Rice's original hires_to_pancam_S.pro program: 10/26/2010 - Ryan Anderson
# 11/1/2010 - Changed effective wavelengths to centroid of the filter. - Ryan Anderson
# 11/22/2010 - Changed to use rd_tfile.pro to read in data. Added delimeter keyword. - Ryan Anderson
# 12/22/2010 - added microns keyword - Ryan Anderson
# 11/15/2011 - adapted to the new Mastcam passbands provided by Jim Bell
# 3/15/2013 - will now output convolved RGB Bayer Filter reflectances
# 6/18/2019 - translated to Python from IDL - Kathleen Hoza
# 10/3/2019 - modified for mastcam-z calibration measurements (no solar input)
# INPUTS: input_file (the ASCII file containing the spectra to convolve.
#                     The first column must be wavelength in nm, and
#                     the following columns must contain reflectance)
#          output_file (the program will save the convolved 
#                      wavelengths and reflectances to the given
#                      filename (must end with '.txt'))
#          ncol (the number of columns in the input file - can be
#                approximated to the nearest power of 10)
# 		  microns - set this keyword if the input spectra have wavelengths in microns. 
# 					Wavelengths will be multiplied by 1000 to convert to nm
# 
# OUTPUTS: a file is written with mastcam wavelengths in the first column and the convolved
#          mastcam spectra in the following columns.
# 
# 
from scipy import interpolate
from scipy import integrate
import numpy as np
import os


#READ IN THE HIGH-RESOLUTION SPECTRA AND WAVELENGTHS
microns=False
input_file='C:\\Users\\krist\\Documents\\Python_Scripts\\MCZ_convolve\\data2convolve\\2022_02_09_Louise_Dunites_Andesites_ReconGeom.csv'
os.chdir('C:\\Users\\krist\\Documents\\Python_Scripts\\MCZ_convolve\\hires_to_MCZ')
input_data=np.genfromtxt(input_file, dtype=float, delimiter=',',skip_header=5, unpack=True) #makes .csv readable and 'unpack=True flips column and row
spectrum_wav=np.array(input_data[0]) #makes first row

if microns:
     spectrum_wav=spectrum_wav*1000
#spectrum_wav=spectrum_wav(where(spectrum_wav gt 0))

n=len(input_data) #number of spectra we are convolving+1 for wavelengths

#DEFINE AN ARRAY FOR THE OUTPUTS
output_data=[[]]

for m in range(1,n):
    #spectrum_orig=reform(double(input_data[m,*]))
    #spectrum_orig=spectrum_orig(where(spectrum_wav gt 0))
    spectrum_orig=np.array(input_data[m])
    # spectrum_orig=[1,-1,2,-2]
    # indices_less_zero=[(i,r) for (i,r) in zip(spectrum_orig,x) if i >= 0]
    # print(indices_less_than_zero)
    
    #READ IN THE SOLAR SPECTRUM AND PUT INTO INPUT SPECTRUM WAVELENGTHS --> no solar spectrum for MCZ (calib measurements indoors)

    #solar_input=np.genfromtxt('C:\\Users\\krist\\OneDrive\\Documents\\Python Scripts\\mastcam_convolve\\hires_to_mastcam\\sun_input.txt', unpack=True)

    #interpol = interpolate.interp1d(solar_input[0],solar_input[1], kind='linear')
    #solar_new=interpol(spectrum_wav) #interpolate values of solar function at spectrum wavelengths
    
    #multiply high-res spectrum by the interpolated solar spectrum --> not needed for MCZ indoor calibration, variable reassigned
    spectrum=spectrum_orig#*solar_new
    
    #fix any crazy negative values in the spectrum - written in IDL, not currently implemented in python
    #zero=where(spectrum LT 0.)
    #IF zero[0] GT 0 THEN spectrum(where(spectrum LT 0.))=0.
    
    #define wavelength array (300-1100 nm with 5nm steps) to use with Mastcam filter data
    wvl=np.arange(350,1050,5)#make_array(161,/index)*5.+300.
    
    #restore merged / normalized Mastcam filter data
    [bayerL0B_wvl, bayerL0B]=np.genfromtxt('mastcamz/normalized_MCZ_L0B.txt',unpack=True)
    [bayerL0R_wvl, bayerL0R]=np.genfromtxt('mastcamz/normalized_MCZ_L0R.txt',unpack=True)
    [bayerL0G_wvl, bayerL0G]=np.genfromtxt('mastcamz/normalized_MCZ_L0G.txt',unpack=True)
    [filterL1_wvl,filterL1]=np.genfromtxt('mastcamz/normalized_MCZ_L1.txt',unpack=True)
    [filterL2_wvl,filterL2]=np.genfromtxt('mastcamz/normalized_MCZ_L2.txt',unpack=True)
    [filterL3_wvl,filterL3]=np.genfromtxt('mastcamz/normalized_MCZ_L3.txt',unpack=True)
    [filterL4_wvl,filterL4]=np.genfromtxt('mastcamz/normalized_MCZ_L4.txt',unpack=True)
    [filterL5_wvl,filterL5]=np.genfromtxt('mastcamz/normalized_MCZ_L5.txt',unpack=True)
    [filterL6_wvl,filterL6]=np.genfromtxt('mastcamz/normalized_MCZ_L6.txt',unpack=True)
    
    [bayerR0B_wvl, bayerR0B]=np.genfromtxt('mastcamz/normalized_MCZ_R0B.txt',unpack=True)
    [bayerR0R_wvl, bayerR0R]=np.genfromtxt('mastcamz/normalized_MCZ_R0R.txt',unpack=True)
    [bayerR0G_wvl, bayerR0G]=np.genfromtxt('mastcamz/normalized_MCZ_R0G.txt',unpack=True)
    [filterR1_wvl,filterR1]=np.genfromtxt('mastcamz/normalized_MCZ_R1.txt',unpack=True)
    [filterR2_wvl,filterR2]=np.genfromtxt('mastcamz/normalized_MCZ_R2.txt',unpack=True)
    [filterR3_wvl,filterR3]=np.genfromtxt('mastcamz/normalized_MCZ_R3.txt',unpack=True)
    [filterR4_wvl,filterR4]=np.genfromtxt('mastcamz/normalized_MCZ_R4.txt',unpack=True)
    [filterR5_wvl,filterR5]=np.genfromtxt('mastcamz/normalized_MCZ_R5.txt',unpack=True)
    [filterR6_wvl,filterR6]=np.genfromtxt('mastcamz/normalized_MCZ_R6.txt',unpack=True)
    
    
    #interpolate the the narrow-band filters and solar spectrum to the same wavelength scale
    interpol = interpolate.interp1d(bayerL0B_wvl,bayerL0B, kind='linear')
    filterL0_blue=interpol(wvl)
    interpol = interpolate.interp1d(bayerL0R_wvl,bayerL0R, kind='linear')
    filterL0_red=interpol(wvl)
    interpol = interpolate.interp1d(bayerL0G_wvl, bayerL0G, kind='linear')
    filterL0_green=interpol(wvl)
    
    interpol = interpolate.interp1d(filterL1_wvl,filterL1, kind='linear')
    filterL1=interpol(wvl)
    
    interpol = interpolate.interp1d(filterL2_wvl,filterL2, kind='linear')
    filterL2=interpol(wvl)
    
    interpol = interpolate.interp1d(filterL3_wvl,filterL3, kind='linear')
    filterL3=interpol(wvl)
    
    interpol = interpolate.interp1d(filterL4_wvl,filterL4, kind='linear')
    filterL4=interpol(wvl)
    
    interpol = interpolate.interp1d(filterL5_wvl,filterL5, kind='linear')
    filterL5=interpol(wvl)
    
    interpol = interpolate.interp1d(filterL6_wvl,filterL6, kind='linear')
    filterL6=interpol(wvl)



    interpol = interpolate.interp1d(bayerR0B_wvl,bayerR0B, kind='linear')
    filterR0_blue=interpol(wvl)
    interpol = interpolate.interp1d(bayerR0R_wvl, bayerR0R, kind='linear')
    filterR0_red=interpol(wvl)
    interpol = interpolate.interp1d(bayerR0G_wvl, bayerR0G, kind='linear')
    filterR0_green=interpol(wvl)
    
    interpol = interpolate.interp1d(filterR1_wvl,filterR1, kind='linear')
    filterR1=interpol(wvl)
    
    interpol = interpolate.interp1d(filterR2_wvl,filterR2, kind='linear')
    filterR2=interpol(wvl)
    
    interpol = interpolate.interp1d(filterR3_wvl,filterR3, kind='linear')
    filterR3=interpol(wvl)
    
    interpol = interpolate.interp1d(filterR4_wvl,filterR4, kind='linear')
    filterR4=interpol(wvl)
    
    interpol = interpolate.interp1d(filterR5_wvl,filterR5, kind='linear')
    filterR5=interpol(wvl)
    
    interpol = interpolate.interp1d(filterR6_wvl,filterR6, kind='linear')
    filterR6=interpol(wvl)


    ## Solar input removed for indoor MCZ calibration
    #interpol=interpolate.interp1d(spectrum_wav,solar_new,kind='linear')
    #solar_interp=interpol(wvl)
    
    
    #Integrate over each filter
    filterL0_red_int=integrate.trapz(filterL0_red,wvl)
    filterL0_green_int=integrate.trapz(filterL0_green,wvl)
    filterL0_blue_int=integrate.trapz(filterL0_blue,wvl)
    filterL1_int=integrate.trapz(filterL1,wvl)
    filterL2_int=integrate.trapz(filterL2,wvl)
    filterL3_int=integrate.trapz(filterL3,wvl)
    filterL4_int=integrate.trapz(filterL4,wvl)
    filterL5_int=integrate.trapz(filterL5,wvl)
    filterL6_int=integrate.trapz(filterL6,wvl)
    
    filterR0_red_int=integrate.trapz(filterR0_red,wvl)
    filterR0_green_int=integrate.trapz(filterR0_green,wvl)
    filterR0_blue_int=integrate.trapz(filterR0_blue,wvl)
    filterR1_int=integrate.trapz(filterR1,wvl)
    filterR2_int=integrate.trapz(filterR2,wvl)
    filterR3_int=integrate.trapz(filterR3,wvl)
    filterR4_int=integrate.trapz(filterR4,wvl)
    filterR5_int=integrate.trapz(filterR5,wvl)
    filterR6_int=integrate.trapz(filterR6,wvl)
    
    
    #normalize each filter so that the integral is 1
    filterL0_red=filterL0_red/filterL0_red_int
    filterL0_green=filterL0_green/filterL0_green_int
    filterL0_blue=filterL0_blue/filterL0_blue_int
    filterL1=filterL1/filterL1_int
    filterL2=filterL2/filterL2_int
    filterL3=filterL3/filterL3_int
    filterL4=filterL4/filterL4_int
    filterL5=filterL5/filterL5_int
    filterL6=filterL6/filterL6_int
    
    filterR0_red=filterR0_red/filterR0_red_int
    filterR0_green=filterR0_green/filterR0_green_int
    filterR0_blue=filterR0_blue/filterR0_blue_int
    filterR1=filterR1/filterR1_int
    filterR2=filterR2/filterR2_int
    filterR3=filterR3/filterR3_int
    filterR4=filterR4/filterR4_int
    filterR5=filterR5/filterR5_int
    filterR6=filterR6/filterR6_int
    

    #interpolate the high-res spectrum to the same wavelength scale as the narrow filters
    interpol=interpolate.interp1d(spectrum_wav,spectrum,kind='linear')
    new_spec_rad = interpol(wvl)
    
    
    #multiply high-res spectrum by each filter, integrate,(and divide out solar spectrum)
    filterL0_red=integrate.trapz(filterL0_red*new_spec_rad, wvl)#/integrate.trapz(filterL0_red*solar_interp,wvl)

    filterL0_green=integrate.trapz(filterL0_green*new_spec_rad, wvl)#/integrate.trapz(filter0_green*solar_interp,wvl)
           
    filterL0_blue=integrate.trapz(filterL0_blue*new_spec_rad, wvl)#/integrate.trapz(filter0_blue*solar_interp,wvl)

    filterL1=integrate.trapz(filterL1*new_spec_rad, wvl)#/integrate.trapz(filter1_525nm*solar_interp,wvl)

    filterL2=integrate.trapz(filterL2*new_spec_rad, wvl)#/integrate.trapz(filter2_440nm*solar_interp,wvl)

    filterL3=integrate.trapz(filterL3*new_spec_rad, wvl)#/integrate.trapz(filter3_750nm*solar_interp,wvl)

    filterL4=integrate.trapz(filterL4*new_spec_rad, wvl)#/integrate.trapz(filter4_905nm*solar_interp,wvl)

    filterL5=integrate.trapz(filterL5*new_spec_rad, wvl)#/integrate.trapz(filter5_865nm*solar_interp,wvl)

    filterL6=integrate.trapz(filterL6*new_spec_rad, wvl)#/integrate.trapz(filter6_1035nm*solar_interp,wvl)


    filterR0_red=integrate.trapz(filterR0_red*new_spec_rad, wvl)#/integrate.trapz(filterL0_red*solar_interp,wvl)

    filterR0_green=integrate.trapz(filterR0_green*new_spec_rad, wvl)#/integrate.trapz(filter0_green*solar_interp,wvl)
           
    filterR0_blue=integrate.trapz(filterR0_blue*new_spec_rad, wvl)#/integrate.trapz(filter0_blue*solar_interp,wvl)

    filterR1=integrate.trapz(filterR1*new_spec_rad, wvl)#/integrate.trapz(filter1_525nm*solar_interp,wvl)

    filterR2=integrate.trapz(filterR2*new_spec_rad, wvl)#/integrate.trapz(filter2_440nm*solar_interp,wvl)

    filterR3=integrate.trapz(filterR3*new_spec_rad, wvl)#/integrate.trapz(filter3_750nm*solar_interp,wvl)

    filterR4=integrate.trapz(filterR4*new_spec_rad, wvl)#/integrate.trapz(filter4_905nm*solar_interp,wvl)

    filterR5=integrate.trapz(filterR5*new_spec_rad, wvl)#/integrate.trapz(filter5_865nm*solar_interp,wvl)

    filterR6=integrate.trapz(filterR6*new_spec_rad, wvl)#/integrate.trapz(filter6_1035nm*solar_interp,wvl)

    #average left and right RGB filters 
    filter_red=(filterL0_red + filterR0_red)/2
    filter_green=(filterL0_green + filterR0_green)/2
    filter_blue=(filterL0_blue + filterR0_blue)/2
    filter1=(filterL1 + filterR1)/2

    #define arrays
    mastcamz_spect=[filterL6,filter_blue,filterL5,filter_green,filterL4,filter_red,filterL3,filterL2,filter1,filterR2,filterR3,filterR4,filterR5,filterR6]
    
    #define mastcam effective wavelengths 
    mastcamz_wav = [440.,480.,528.,544.,605.,630.5,677.,754.,800.,866.,910.,939.,978.,1022.]
    
    #sort the mastcamz arrays in ascending order
    mastcamz=[x for _, x in sorted(zip(mastcamz_wav,mastcamz_spect), key=lambda pair: pair[0])]
    
    output_data.append(mastcamz)

mastcamz_wav=sorted(mastcamz_wav)
output_data[0]=mastcamz_wav
headers=[]
with open(input_file,'r') as f:
     for i in range(5):
          headers.append(f.readline())
print(os.getcwd())
with open(input_file.strip('.csv')+'_convolved.csv','w+') as f:
     for line in headers:
          f.write(line)
     for i in range(len(output_data[0])):
          for j in range(len(output_data)):
               f.write(str(output_data[j][i])+',')
          f.write('\n')

