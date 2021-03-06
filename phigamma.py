# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 13:42:18 2015

@author: hareyakana
"""
import matplotlib.pylab as pl
#import time 
import math as m
import numpy as np
#import itertools as i

ID=list(range(1,38))
"""Event type"""
S,T,C=1,0,2 
topology=[S,S,T,S,T,S,S,T,S,S,S,S,T,S,S,S,S,T,S,S,S,S,T,S,S,S,S,T,S,S,S,C,S,S,S,S,T,T,S,S,S,S,T,T,T,S,T,S,S,S,S,S,T,S]

"""Deposited Enery in TeV obtained from the 4 years data """
E=np.array([47.6, 117, 78.7, 165.4, 71.4, 28.4, 34.3, 32.6, 63.2, 97.2, 88.4, 104.1, 252.7,
   1040.7, 57.5, 30.6, 199.7, 31.5, 71.5, 1140.8, 30.2, 219.5, 82.2, 30.5, 33.5, 210,
   60.2, 46.1, 32.7, 128.7, 42.5, 0,384.7, 42.1, 2003.7, 28.9, 30.8, 200.5, 101.3, 157.3,
   87.6, 76.3, 46.5, 84.6, 429.9, 158, 74.3, 104.7, 59.9, 22.2, 66.2, 158.1, 27.6,
   54.5])
e=E*1000.

"""Background"""
bkmuoni=np.array([0,0.05,1.1,0.6,0.33,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10])
bkmuon=10**bkmuoni
bkatmi=np.array([0,0,0.35,0.34,0.25,-0.1,-0.6,-1.1,-1.9,-10,-10,-10,-10,-10,-10])
bkatm=10**bkatmi
atm=np.array([0,0.2,0.62,0.75,0.72,0.48,0.2,-0.3,-0.65,-1.1,-1.52,-10,-10,-10,-10])
atmnu=10**atm
background=bkmuon+bkatm+atmnu
#####
ET1=np.linspace(4.0,6.8,num=15)
ET2=np.linspace(4.2,7.0,num=15)
et1=10**ET1
et2=10**ET2
ET=np.linspace(4.1,6.9,num=15)
et=(et1+et2)/2


bins=np.zeros(len(ET1))  ###DATA

for i in range(len(E)):
    for k in range(len(ET1)):
        if abs(e[i])>abs(et1[k]) and abs(e[i])<abs(et2[k]):
            bins[k]=bins[k]+1
    
"""Theory model/etc"""
time=1347*24*60*60  #total IceCube sample time

"""POWER LAW"""

def P(E,n,g,e,a,bk):  # of unit s-1 cm-2
    p=n*1e-18*((E/1e5)**-g)
    g=p*time*4*m.pi*a*1e4*(E-e)+bk
    return g
    
#PROM/CONVENTIONAL SAME FORMULA BUT DIFFERENT PARAMETER
#th=P(et)*time*4*m.pi


"""Effective area"""
#ae=np.array([])
#am=np.array([])
#at=np.array([])

"""from public data"""
import csv
e1=[]
e2=[]
e3=[]

#with open('effective_area.nu_e.txt','r')as f:
#    for row in f:
#        x,y,h,k,z,p=row.split()
#        e1.append(float(x))
#        e2.append(float(y))
#        e3.append(float(z))

with open('nue_4pi.txt','r')as f:
    reader=csv.DictReader(f)
    for row in reader:
        e1.append(float(row['x']))
        e2.append(float(row['y']))
        e3.append(float(row['z']))

m1=[]
m2=[]
m3=[]

#with open('effective_area.nu_mu.txt','r')as f:
#    for row in f:
#        x,y,h,k,z,p=row.split()
#        m1.append(float(x))
#        m2.append(float(y))
#        m3.append(float(z))

with open('numu_4pi.txt','r')as f:
    reader=csv.DictReader(f)
    for row in reader:
        m1.append(float(row['x']))
        m2.append(float(row['y']))
        m3.append(float(row['z']))        

t1=[]
t2=[]
t3=[]

#with open('effective_area.nu_tau.txt','r')as f:
#    for row in f:
#        x,y,h,k,z,p=row.split()
#        t1.append(float(x))
#        t2.append(float(y))
#        t3.append(float(z))

with open('nutau_4pi.txt','r')as f:
    reader=csv.DictReader(f)
    for row in reader:
        t1.append(float(row['x']))
        t2.append(float(row['y']))
        t3.append(float(row['z']))

""""""        
area=np.zeros(len(et1),dtype=float)
repeat=np.zeros(len(et1),dtype=float)

for i in range(len(et1)):
    for k in range(len(e1)):
        if abs(e2[k])>abs(et1[i]) and abs(e2[k])<abs(et2[i]):
            area[i]=area[i]+e3[k]
            repeat[i]=repeat[i]+1
        if abs(m2[k])>abs(et1[i]) and abs(m2[k])<abs(et2[i]):
            area[i]=area[i]+m3[k]
            repeat[i]=repeat[i]+1
        if abs(t2[k])>abs(et1[i]) and abs(t2[k])<abs(et2[i]):
            area[i]=area[i]+t3[k]
            repeat[i]=repeat[i]+1

for x in range(0,7):
    area[2*x]=area[2*x]*2.5/2.

for x in range(0,8):
    area[x*2-1]=area[x*2-1]*2.5/3.
    
#########################
"""best fit parameter"""
#bins[-4]=bins[-5]=0

#phi=2.2 #normalization index=8.4, astro flux=2.2
#gam=np.linspace(2,5,num=1001)
#
#test=[]
#import scipy as sp
#
#for i in range(0,len(gam)):
#    theory=P(et2,phi,gam[i])*time*4*np.pi*area*(et2-et1)+background
#    chis=((bins-theory)*(bins-theory))/theory
#    test.append(sum(chis))
#    
#pl.plot(gam,test)
#pl.show()
#print("minimum chi squared",min(test))
#print("minimum spectral index gamma=",gam[test.index(min(test))])
#
#gamma=gam[test.index(min(test))]
#flux=np.linspace(0,20,num=1001)
#
#test2=[]
#for j in range(0,len(flux)):
#    theory=P(et2,flux[j],gamma)*time*4*np.pi*area*(et2-et1)+background
#    chis=((bins-theory)*(bins-theory))/theory
#    test2.append(sum(chis))
#    
#pl.subplot
#pl.plot(flux,test2)
#pl.show()
#print("minimum chi squared",min(test2))
#print("minimum flux, phi=",flux[test2.index(min(test2))])
    
#########################
#########################
phi=np.linspace(0.1,10,num=501)
gamma=np.linspace(2,4,num=501)

p,g=np.meshgrid(phi,gamma)

def chi(observed,expected):
    ch=(observed-expected)*(observed-expected)/expected   
    return ch

blank=[]
for i in range(0,len(et)):
    theo=P(et2[i],p,g,et1[i],area[i],background[i])
    theory=chi(bins[i],theo)
    blank.append(theory)

total=blank[0]+blank[1]+blank[2]+blank[3]+blank[4]+blank[5]+blank[6]+blank[7]+blank[8]+blank[9]+blank[10]+blank[11]+blank[12]+blank[13]+blank[14]

import matplotlib.gridspec as gridspec

h, w = total.shape
gs = gridspec.GridSpec(2,2,width_ratios=[w,w*.2], height_ratios=[h,h*.2])
ax = [pl.subplot(gs[0]),]
ax.append(pl.subplot(gs[1], sharey=ax[0]))
ax.append(pl.subplot(gs[2], sharex=ax[0]))
bounds = [phi.min(),phi.max(),gamma.min(),gamma.max()]
ax[0].imshow(total, cmap='Blues_r', extent=bounds,origin='lower')
ax[1].plot(total[:,int(w/2)],g[:,int(w/2)],'.',total[:,int(w/2)],g[:,int(w/2)])
#ax[1].axis([total[:,int(w/2)].max(), total[:,int(w/2)].min(), g.min(), g.max()])
ax[2].plot(p[int(h/2),:],total[int(h/2),:],'.',p[int(h/2),:],total[int(h/2),:])
pl.show()

p1=[]
invalue=[]
x=[]
y=[]



for i in range(0,len(phi)):
    for k in range(0,len(gamma)):
        theory=P(et2,phi[i],gamma[k],et1,area,background)
        chisquared=chi(bins,theory)
        p1.append(sum(chisquared))
        invalue.append([phi[i],gamma[k]])
        x.append(phi[i])
        y.append(gamma[k])

print(invalue[p1.index(min(p1))],min(p1))
