# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 16:00:46 2021
Case OST_SHM_cL
@author: romai
"""

import numpy as np
import matplotlib.pyplot as plt
import os
#import GeneralProperties_SE as globalProperties
#import mySimpleModel as QSModel

#Nt = 1000     # [-] # Number of time steps

#globalProperties.init() # initialize global structure
#globalProperties.properties['wing']['surface'] = np.trapz(QSModel.chordDistribution(globalProperties.properties['wing']['shape']),dx = globalProperties.properties['wing']['span']/(discretisationRate-1))
#S =np.mean(QSModel.chordDistribution(globalProperties.properties['wing']['shape']))*globalProperties.properties['wing']['span']
#cmean = np.mean(QSModel.chordDistribution(globalProperties.properties['wing']['shape']))

#%% Main user settings

mmM = 1000 #  mm to m
factorChord = 2
#Xcut = 0.0001
thickness = 0.5
discretisationRate = 5000

t2 = thickness/2

#fname = '../../../All/system/domainSetting'
fname = '../All/system/domainSetting'
def read_domainSetting(fname):
    file = np.loadtxt(fname,dtype=str)
    #print(file)
    for i in range(len(file)):
        if file[i,0] == 'deltaR':
            deltaR = float(file[i,1][:-1])
        elif file[i,0] == 'span':
            span = float(file[i,1][:-1])
        elif file[i,0] == 'chord':
            chord = float(file[i,1][:-1])
        elif file[i,0] == 'flap_freq':
            flap_freq = float(file[i,1][:-1])
        elif file[i,0] == 'A_phi':
            A_phi = float(file[i,1][:-1])
        elif file[i,0] == 'A_alpha':
            A_alpha = float(file[i,1][:-1])
        elif file[i,0] == 'alpha_dev':
            alpha_dev = float(file[i,1][:-1])
        elif file[i,0] == 'K_phi':
            K_phi = float(file[i,1][:-1])
        elif file[i,0] == 'K_alpha':
            K_alpha = float(file[i,1][:-1])

            
    return deltaR, span, chord, flap_freq, A_phi, A_alpha, alpha_dev,K_phi, K_alpha

#read input file
deltaR, span, chord, flap_freq, A_phi, A_alpha, alpha_dev,K_phi, K_alpha = read_domainSetting(fname)

span = span/1000
#span = 0.050 #globalProperties.properties['wing']['span']; 
chord = chord/1000
chord_root = chord*4/np.pi;
thickOverset = factorChord*chord;

#%%

factorSpanCut= 0.3;
Xtip = span - span*factorSpanCut;
Ztip = np.sqrt((1 - np.square(Xtip/span)))*chord_root/2;
XpeakLos = Xtip - Ztip;
chord_root2 = chord_root/2 + thickOverset
s2 = span + thickOverset
X2tip = (2*XpeakLos + np.sqrt( 4*np.square(XpeakLos) - 4*(np.square(XpeakLos) - np.square(chord_root2))*(1+np.square(chord_root2)/np.square(s2)) )) / (2*(1+np.square(chord_root2)/np.square(s2)))
Z2tip = X2tip - XpeakLos;

Xvecmm = np.linspace(0,Xtip,discretisationRate)
Xvec = np.linspace(0,Xtip,discretisationRate)*mmM
z = np.sqrt(np.abs(1 - (Xvecmm/span)**2))*chord_root/2*mmM

Xvecmm_P2 = np.linspace(Xtip,span,discretisationRate)
Xvec_P2 = np.linspace(Xtip,span,discretisationRate)*mmM
z_P2 = np.sqrt(np.abs(1 - (Xvecmm_P2/span)**2))*chord_root/2*mmM

x2mm = np.linspace(0,X2tip,discretisationRate)
x2   = np.linspace(0,X2tip,discretisationRate)*mmM
z2 = (np.sqrt(np.abs(1 - (x2mm/s2)**2))*chord_root2)*mmM

x2mm_P2 = np.linspace(X2tip,s2,discretisationRate)
x2_P2   = np.linspace(X2tip,s2,discretisationRate)*mmM
z2_P2 = (np.sqrt(np.abs(1 - (x2mm_P2/s2)**2))*chord_root2)*mmM


yp = np.ones(discretisationRate)*t2
ym = -np.ones(discretisationRate)*t2


y2p = np.ones(discretisationRate)*t2
y2m = -np.ones(discretisationRate)*t2

ypp = np.ones(discretisationRate)*(t2+factorChord*chord*mmM)
ymm = np.ones(discretisationRate)*(-t2-factorChord*chord*mmM)


#change directory to system/edges
os.chdir('system/edges')

#%%
###########################
## YLevel0  - z+ 
f= open("2_12.txt","w+")
f.write("polyLine 2 12\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec[i],yp[i],z[i]))
f.write(") \r")   
f.close();
f= open("6_14.txt","w+")
f.write("polyLine 6 14\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec[i],ym[i],z[i]))
f.write(") \r")   
f.close();
## YLevel0  - z++
f= open("8_13.txt","w+")
f.write("polyLine 8 13\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2[i],y2p[i],z2[i]))
f.write(") \r")   
f.close();
f= open("10_15.txt","w+")
f.write("polyLine 10 15\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2[i],y2m[i],z2[i]))
f.write(") \r")   
f.close();

############  P2
## YLevel0  - z+ 
f= open("12_16.txt","w+")
f.write("polyLine 12 16\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec_P2[i],yp[i],z_P2[i]))
f.write(") \r")   
f.close();
f= open("14_18.txt","w+")
f.write("polyLine 14 18\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec_P2[i],ym[i],z_P2[i]))
f.write(") \r")   
f.close();
## YLevel0  - z++
f= open("13_17.txt","w+")
f.write("polyLine 13 17\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2_P2[i],y2p[i],z2_P2[i]))
f.write(") \r")   
f.close();
f= open("15_19.txt","w+")
f.write("polyLine 15 19\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2_P2[i],y2m[i],z2_P2[i]))
f.write(") \r")   
f.close();



#%% #########################
## YLevel0  - z- 
f= open("21_28.txt","w+")
f.write("polyLine 21 28\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec[i],yp[i],-z[i]))
f.write(") \r")   
f.close();
f= open("23_30.txt","w+")
f.write("polyLine 23 30\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec[i],ym[i],-z[i]))
f.write(") \r")   
f.close();
## YLevel0  - z--
f= open("25_29.txt","w+")
f.write("polyLine 25 29\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2[i],y2p[i],-z2[i]))
f.write(") \r")   
f.close();
f= open("27_31.txt","w+")
f.write("polyLine 27 31\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2[i],y2m[i],-z2[i]))
f.write(") \r")   
f.close();

############  P2
## YLevel0  - z- 
f= open("28_16.txt","w+")
f.write("polyLine 28 16\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec_P2[i],yp[i],-z_P2[i]))
f.write(") \r")   
f.close();
f= open("30_18.txt","w+")
f.write("polyLine 30 18\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec_P2[i],ym[i],-z_P2[i]))
f.write(") \r")   
f.close();
## YLevel0  - z--
f= open("29_17.txt","w+")
f.write("polyLine 29 17\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2_P2[i],y2p[i],-z2_P2[i]))
f.write(") \r")   
f.close();
f= open("31_19.txt","w+")
f.write("polyLine 31 19\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2_P2[i],y2m[i],-z2_P2[i]))
f.write(") \r")   
f.close();

#%% 2
###########################
## YLevel1  - z+ 
f= open("34_38.txt","w+")
f.write("polyLine 34 38\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec[i],ypp[i],z[i]))
f.write(") \r")   
f.close();
## YLevel0  - z++
f= open("36_39.txt","w+")
f.write("polyLine 36 39\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2[i],ypp[i],z2[i]))
f.write(") \r")   
f.close();

############  P2
## YLevel0  - z+ 
f= open("38_40.txt","w+")
f.write("polyLine 38 40\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec_P2[i],ypp[i],z_P2[i]))
f.write(") \r")   
## YLevel0  - z++
f= open("39_41.txt","w+")
f.write("polyLine 39 41\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2_P2[i],ypp[i],z2_P2[i]))
f.write(") \r")   
f.close();



#%% #########################
## YLevel1  - z- 
f= open("43_46.txt","w+")
f.write("polyLine 43 46\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec[i],ypp[i],-z[i]))
f.write(") \r")   
f.close();
## YLevel0  - z--
f= open("45_47.txt","w+")
f.write("polyLine 45 47\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2[i],ypp[i],-z2[i]))
f.write(") \r")   
f.close();


############  P2
## YLevel0  - z- 
f= open("46_40.txt","w+")
f.write("polyLine 46 40\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec_P2[i],ypp[i],-z_P2[i]))
f.write(") \r")   
f.close();
## YLevel0  - z--
f= open("47_41.txt","w+")
f.write("polyLine 47 41\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2_P2[i],ypp[i],-z2_P2[i]))
f.write(") \r")   
f.close();

#%% 3
###########################
## YLevel-1  - z+ 
f= open("53_57.txt","w+")
f.write("polyLine 53 57\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec[i],ymm[i],z[i]))
f.write(") \r")   
f.close();
## YLevel0  - z++
f= open("55_58.txt","w+")
f.write("polyLine 55 58\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2[i],ymm[i],z2[i]))
f.write(") \r")   
f.close();

############  P2
## YLevel0  - z+ 
f= open("57_59.txt","w+")
f.write("polyLine 57 59\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec_P2[i],ymm[i],z_P2[i]))
f.write(") \r")   
## YLevel0  - z++
f= open("58_60.txt","w+")
f.write("polyLine 58 60\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2_P2[i],ymm[i],z2_P2[i]))
f.write(") \r")   
f.close();



#%% #########################
## YLevel1  - z- 
f= open("62_65.txt","w+")
f.write("polyLine 62 65\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec[i],ymm[i],-z[i]))
f.write(") \r")   
f.close();
## YLevel0  - z--
f= open("64_66.txt","w+")
f.write("polyLine 64 66\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2[i],ymm[i],-z2[i]))
f.write(") \r")   
f.close();


############  P2
## YLevel0  - z- 
f= open("65_59.txt","w+")
f.write("polyLine 65 59\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (Xvec_P2[i],ymm[i],-z_P2[i]))
f.write(") \r")   
f.close();
## YLevel0  - z--
f= open("66_60.txt","w+")
f.write("polyLine 66 60\n")
f.write("( \n")
for i in range(discretisationRate):
     f.write("( %.8f %.8f %.8f)\n" % (x2_P2[i],ymm[i],-z2_P2[i]))
f.write(") \r")   
f.close();





#plt.figure();
#plt.plot(Xvec,z,'b')
#plt.plot(Xvec_P2,z_P2,'r')
#plt.plot(x2,z2,'b')
#plt.plot(x2_P2,z2_P2,'r')


