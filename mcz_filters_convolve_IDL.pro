; PROGRAM NAME: mcz_filters_convolve.pro
; WRITTEN BY: Melissa Rice
; DATE: 8/19/2019
; PURPOSE: To read in high resolution spectra and convolve them to Mastcam-Z bandpasses
; NOTES: Bandpasses are taken from the spectral throughput data from MCZ calibration (May 2019) - NOT YET FINAL 
;            Bayer RGB and 805 bands all taken from LEFT camera (R0R,R0G,R0B and R1 excluded from convolutions below) 
;            Bandpass data taken as InBand + OutBand merged data for the following bayer channels:
;                L6: blue channel
;                L5: green channel
;                L4-1: red channel
;                R2-6: green channel
;        This version does NOT include the solar spectrum - that will be added before surface ops 
;             
; 
;This script is adapted from a 2010 Mastcam script, which was adapted from a 2003(?) Pancam script

; INPUTS: input_file (the ASCII .txt file containing the spectra to convolve.
;                    The first column must be wavelength in nm, and
;                    the following columns must contain reflectance.
;					 Note that there should NOT be any header row(s).)
;         output_file (the program will save the convolved 
;                     wavelengths and reflectances to the given
;                     filename (must end with '.txt'))
;         ncol (the number of columns in the input file - can be
;               approximated to the nearest power of 10)
;		  microns - set this keyword if the input spectra have wavelengths in microns. 
;					Wavelengths will be multiplied by 1000 to convert to nm
;
; OUTPUTS: a file is written with Mastcam-Z wavelengths in the first column and the convolved
;          Mastcam-Z spectra in the following columns.


; DEFINE THE PROGRAM
PRO mcz_filters_convolve,input_file,output_file,ncol,microns=microns

; READ IN THE HIGH-RESOLUTION SPECTRA AND WAVELENGTHS
input_data=READ_ASCII(input_file)
IF ncol LE 10 THEN BEGIN 
input_data=input_data.FIELD1
ENDIF ELSE BEGIN 
input_data=input_data.FIELD01
ENDELSE  
IF ncol GT 100 THEN BEGIN
input_data=input_data.FIELD001
ENDIF

spectrum_wav=reform(double(input_data[0,*]))
if keyword_set(microns) then spectrum_wav=spectrum_wav*1000.
spectrum_wav=spectrum_wav(where(spectrum_wav gt 0))

s=size(input_data)
n=s[1]
PRINT,n-1

; DEFINE AN ARRAY FOR THE OUTPUTS
output_data=fltarr(n,14)

;define wavelength array (300-1100 nm with 5nm steps) to use with Mastcam-Z filter data
wave=make_array(161,/index)*5.+300.

;restore normalized, merged in-band and out-of-band Mastcam-Z filter data
restore,'mcamz_cal/normalized_MCZ_R2.sav'
R2wvl=wvl
R2=tputnorm
restore,'mcamz_cal/normalized_MCZ_R3.sav'
R3wvl=wvl
R3=tputnorm
restore,'mcamz_cal/normalized_MCZ_R4.sav'
R4wvl=wvl
R4=tputnorm
restore,'mcamz_cal/normalized_MCZ_R5.sav'
R5wvl=wvl
R5=tputnorm
restore,'mcamz_cal/normalized_MCZ_R6.sav'
R6wvl=wvl
R6=tputnorm
restore,'mcamz_cal/normalized_MCZ_L0R.sav'
L0Rwvl=wvl
L0R=tputnorm
restore,'mcamz_cal/normalized_MCZ_L0G.sav' 
L0Gwvl=wvl
L0G=tputnorm
restore,'mcamz_cal/normalized_MCZ_L0B.sav'
L0Bwvl=wvl
L0B=tputnorm
restore,'mcamz_cal/normalized_MCZ_L1.sav'
L1wvl=wvl
L1=tputnorm
restore,'mcamz_cal/normalized_MCZ_L2.sav'
L2wvl=wvl
L2=tputnorm
restore,'mcamz_cal/normalized_MCZ_L3.sav'
L3wvl=wvl
L3=tputnorm
restore,'mcamz_cal/normalized_MCZ_L4.sav'
L4wvl=wvl
L4=tputnorm
restore,'mcamz_cal/normalized_MCZ_L5.sav'
L5wvl=wvl
L5=tputnorm
restore,'mcamz_cal/normalized_MCZ_L6.sav'
L6wvl=wvl
L6=tputnorm

;interpolate the the narrow-band filters and solar spectrum to the same wavelength scale
L0R=interpol(L0R,L0Rwvl,wave)
L0G=interpol(L0G,L0Gwvl,wave)
L0B=interpol(L0B,L0Bwvl,wave)
L1=interpol(L1,L1wvl,wave)
L2=interpol(L2,L2wvl,wave)
L3=interpol(L3,L3wvl,wave)
L4=interpol(L4,L4wvl,wave)
L5=interpol(L5,L5wvl,wave)
L6=interpol(L6,L6wvl,wave)
R2=interpol(R2,R2wvl,wave)
R3=interpol(R3,R3wvl,wave)
R4=interpol(R4,R4wvl,wave)
R5=interpol(R5,R5wvl,wave)
R6=interpol(R6,R6wvl,wave)


;Integrate over each filter
L0R_int=int_tabulated(wave,L0R)
L0G_int=int_tabulated(wave,L0G)
L0B_int=int_tabulated(wave,L0B)
L1_int=int_tabulated(wave,L1)
L2_int=int_tabulated(wave,L2)
L3_int=int_tabulated(wave,L3)
L4_int=int_tabulated(wave,L4)
L5_int=int_tabulated(wave,L5)
L6_int=int_tabulated(wave,L6)
R2_int=int_tabulated(wave,R2)
R3_int=int_tabulated(wave,R3)
R4_int=int_tabulated(wave,R4)
R5_int=int_tabulated(wave,R5)
R6_int=int_tabulated(wave,R6)

;normalize each filter so that the integral is 1
L0R=L0R/L0R_int
L0G=L0G/L0G_int
L0B=L0B/L0B_int
L1=L1/L1_int
L2=L2/L2_int
L3=L3/L3_int
L4=L4/L4_int
L5=L5/L5_int
L6=L6/L6_int
R2=R2/R2_int
R3=R3/R3_int
R4=R4/R4_int
R5=R5/R5_int
R6=R6/R6_int

FOR m=1,n-1 DO BEGIN

spectrum_orig=reform(double(input_data[m,*]))
spectrum_orig=spectrum_orig(where(spectrum_wav gt 0))

spectrum=spectrum_orig

; fix any crazy negative values in the spectrum (set them to the adjacent reflectance value)
zero=where(spectrum LT 0.)
IF zero[0] GE 0 THEN spectrum(where(spectrum LT 0.))=spectrum(where(spectrum LT 0.)+1)


; interpolate the high-res spectrum to the same wavelength scale as the narrow filters
new_spec_rad = interpol(spectrum,spectrum_wav,wave)


;multiply high-res spectrum by each filter, integrate
L0R_out=int_tabulated(wave,L0R*new_spec_rad)
L0G_out=int_tabulated(wave,L0G*new_spec_rad)
L0B_out=int_tabulated(wave,L0B*new_spec_rad)
L1_out=int_tabulated(wave,L1*new_spec_rad)
L2_out=int_tabulated(wave,L2*new_spec_rad)
L3_out=int_tabulated(wave,L3*new_spec_rad)
L4_out=int_tabulated(wave,L4*new_spec_rad)
L5_out=int_tabulated(wave,L5*new_spec_rad)
L6_out=int_tabulated(wave,L6*new_spec_rad)
R2_out=int_tabulated(wave,R2*new_spec_rad)
R3_out=int_tabulated(wave,R3*new_spec_rad)
R4_out=int_tabulated(wave,R4*new_spec_rad)
R5_out=int_tabulated(wave,R5*new_spec_rad)
R6_out=int_tabulated(wave,R6*new_spec_rad)


; define arrays
mastcam_spect=[L6_out,L0B_out,L5_out,L0G_out,L4_out,L0R_out,L3_out,L2_out,L1_out,R2_out,R3_out,R4_out,R5_out,R6_out]

; define mastcam-z effective wavelengths 
mastcam_wav = [441.,471.,529.,546.,605.,638.,678.,754.,801.,866.,910.,940.,979.,1012.]
plot,mastcam_wav,mastcam_spect
output_data[m,*]=mastcam_spect

ENDFOR

output_data[0,*]=mastcam_wav

openw,lun,output_file,/get_lun,width=n*13
printf,lun,output_data
free_lun,lun


  END

