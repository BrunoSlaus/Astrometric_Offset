# -*- coding: utf-8 -*-
"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Created on Fri May 12 16:03:09 2017


            
Author: Bruno Slaus
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
from Topcat_Match import skymatch
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import math
import astropy.units as u

############################################################################
#                    Sliding Parameters:                                   # 
############################################################################
Input_Folder     = 'Input/'      #Name of the folder with the input fields
Output_Folder    = 'Output/'     #Name of the folder where the matched field is stored

Field_Name_1     = 'XXL-N_gmrt_out_CorrIRACch1.fits'
Field_Name_2     = 'XXL-N_irac_ch1.fits'

Match_Radius     = 1             #Match radius in arcsec

#WARNING: Write all coulmn names in lower case
#Warning: The unit should be deg
Ra_1_Column  = 'ra'
Dec_1_Column = 'dec'
Ra_2_Column  = 'ra'
Dec_2_Column = 'dec'

N_Bins    = 28
Range_Min = 'auto'    #Histogram min range
Range_Max = 'auto'    #Histogram max range

############################################################################
print('\n******************************')
print('Starting the Astrometric_Corr.py code.')
print('******************************\n')

Field_1 = Input_Folder+Field_Name_1
Field_2 = Input_Folder+Field_Name_2
print('Field_1      = ',Field_1)
print('Field_2      = ',Field_2)
print('Match_Radius = ',Match_Radius)
print('Ra_1_Column  = ',Ra_1_Column)
print('Dec_1_Column = ',Dec_1_Column)
print('Ra_2_Column  = ',Ra_2_Column)
print('Dec_2_Column = ',Dec_2_Column)


#Matching the two fields
skymatch(
    Field_1,
    [Ra_1_Column, Dec_1_Column],
    Field_2,
    [Ra_2_Column, Dec_2_Column],
    Match_Radius, Output_Folder+'Matched.fits')

Matched_Field = fits.open(Output_Folder+'Matched.fits')[1].data

if Ra_1_Column == Ra_2_Column:
    Ra_1_Column = Ra_1_Column + '_1'
    Ra_2_Column = Ra_2_Column + '_2'
if Dec_1_Column == Dec_2_Column:
    Dec_1_Column = Dec_1_Column + '_1'
    Dec_2_Column = Dec_2_Column + '_2'

if Range_Min == 'auto':
    Range_Min = -Match_Radius-(Match_Radius*2 / N_Bins *4)
if Range_Max == 'auto':
    Range_Max =  Match_Radius+(Match_Radius*2 / N_Bins *4)    

Ra_Offset  = (Matched_Field[Ra_1_Column]  - Matched_Field[Ra_2_Column])  * np.cos( (Matched_Field[Dec_1_Column]*u.deg).to(u.rad).value )
Dec_Offset =  Matched_Field[Dec_1_Column] - Matched_Field[Dec_2_Column]
Ra_Offset  = Ra_Offset  * u.deg
Dec_Offset = Dec_Offset * u.deg

Ra_Offset  = Ra_Offset.to(u.arcsec).value   #To je puta 3600
Dec_Offset = Dec_Offset.to(u.arcsec).value
 
Ra_Offset_F,  Ra_Offset_Edges  = np.histogram(Ra_Offset, bins=N_Bins, range=(Range_Min,Range_Max))
Dec_Offset_F, Dec_Offset_Edges = np.histogram(Dec_Offset, bins=N_Bins, range=(Range_Min,Range_Max))

Mean_Ra_Offset =  np.mean(Ra_Offset)
Mean_Dec_Offset = np.mean(Dec_Offset)

if Mean_Ra_Offset>=0:
    Mean_Ra_Offset_String =  str(Mean_Ra_Offset)[:6]
else:
    Mean_Ra_Offset_String =  str(Mean_Ra_Offset)[:7]
if Mean_Dec_Offset>=0:
    Mean_Dec_Offset_String = str(Mean_Dec_Offset)[:6]
else:
    Mean_Dec_Offset_String = str(Mean_Dec_Offset)[:7]

#Offset without projection (i.e. cosine-element)
Ra_Offset_No_Projection  = Matched_Field[Ra_1_Column]  - Matched_Field[Ra_2_Column]
Dec_Offset_No_Projection = Matched_Field[Dec_1_Column] - Matched_Field[Dec_2_Column]
Mean_Ra_No_Projection    = np.mean(Ra_Offset_No_Projection)
Mean_Dec_No_Projection   = np.mean(Dec_Offset_No_Projection)



fig = plt.figure(figsize=(4,6))

ax = fig.add_axes([0.1,0.1,0.40,0.3])
plt.step(Ra_Offset_Edges[:-1], Ra_Offset_F, where='post')
plt.xlabel('Ra_Offset')
plt.axvline(0, color='black', ls='--')
plt.tick_params(axis='y', which='both', labelsize=7)

ax = fig.add_axes([0.5,0.1,0.40,0.3])
plt.step(Dec_Offset_Edges[:-1], Dec_Offset_F, where='post')
ax.yaxis.set_ticks_position('right')
plt.xlabel('Dec_Offset')
plt.axvline(0, color='black', ls='--')
plt.tick_params(axis='y', which='both', labelsize=7)

ax = fig.add_axes([0.1,0.4,0.8,0.5])
plt.scatter(Ra_Offset, Dec_Offset, s=0.3)
plt.scatter(Mean_Ra_Offset, Mean_Dec_Offset, s=26, color='red', marker="+")
plt.axhline(0, color='black')
plt.axvline(0, color='black')
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position("top")
ax.yaxis.set_label_position("right")
#plt.xlim(-Match_Radius-(Match_Radius/10), Match_Radius+(Match_Radius/10))
#plt.ylim(-Match_Radius-(Match_Radius/10), Match_Radius+(Match_Radius/10))
ax.text(0.99, 0.99, r'$\left<\Delta RA\right> = $'+Mean_Ra_Offset_String+'\n'+r'$\left<\Delta DEC\right> = $'+Mean_Dec_Offset_String, color='red', transform=ax.transAxes, fontsize=8, horizontalalignment='right', verticalalignment='top')
plt.xlabel('Ra_Offset')
plt.ylabel('Dec_Offset')

plt.savefig(Output_Folder + 'Astrometric_Errors.png', dpi=300)
plt.close(fig)


print('\n***********************************************')
print('FINISH: The plot and the matched \nfield are saved in Output.\n')
print('Mean_Ra_Offset  == ', Mean_Ra_Offset,  ' arcsec')
print('Mean_Dec_Offset == ', Mean_Dec_Offset, ' arcsec')
print('\nOffsets without projection corrections:')
print('Mean_Ra_Offset_No_Projection  == ', Mean_Ra_No_Projection,  'deg')
print('Mean_Dec_Offset_No_Projection == ', Mean_Dec_No_Projection, 'deg')
print('***********************************************\n')

