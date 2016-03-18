# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 13:42:18 2015

@author: hareyakana
"""
import pylab as pl
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
bkmuoni=np.array([0,0.1,1.1,0.6,0.33,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10])
bkmuon=10**bkmuoni
bkatmi=np.array([0,0,0.35,0.345,0.25,-0.1,-0.6,-1.1,-1.9,-10,-10,-10,-10,-10,-10])
bkatm=10**bkatmi
atm=np.array([0,0.15,0.62,0.75,0.72,0.48,0.2,-0.3,-0.65,-1.1,-1.52,-10,-10,-10,-10])
atmnu=10**atm
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

"""2 POWER LAW
best phi = 2.575
gamma<1PeV= 2.88
gamma>PeV=2.6

"""  

def P1(e):  # of unit s-1 cm-2
    g=2.88  #1PeV
    n=2.275e-18  
    p=n*((e/1e5)**(-g))
    return p
    
def P2(e): 
    g=2.6  #>PeV
    n=2.275e-18  
    p=n*((e/1e5)**(-g))
    return p
    
def Q(e):
    g=2
    n=0.64e-18 ##phi best fit from IceCube
    p=n*((e/1e5)**(-g))
    return p
    
#PROM/CONVENTIONAL SAME FORMULA BUT DIFFERENT PARAMETER
#th=P(et)*time*4*m.pi

"""Effective area"""
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

for i in range(len(et1)):
    for k in range(len(e1)):
        if abs(e2[k])>abs(et1[i]) and abs(e2[k])<abs(et2[i]):
            area[i]=area[i]+e3[k]
        if abs(m2[k])>abs(et1[i]) and abs(m2[k])<abs(et2[i]):
            area[i]=area[i]+m3[k]
        if abs(t2[k])>abs(et1[i]) and abs(t2[k])<abs(et2[i]):
            area[i]=area[i]+t3[k]

for x in range(0,7):
    area[2*x]=area[2*x]*1/2

for x in range(0,8):
    area[x*2-1]=area[x*2-1]*1/3

theory=[]
for i in range(0,len(et2)):
    if et2[i] <= 1e6:
        theory.append(P1(et2[i])*time*4*np.pi*area[i]*1e4*(et2[i]-et1[i])+bkmuon[i]+bkatm[i]+atmnu[i])
    if et2[i] > 1e6:
        theory.append(P2(et2[i])*time*4*np.pi*area[i]*1e4*(et2[i]-et1[i])+bkmuon[i]+bkatm[i]+atmnu[i])
theoryh=Q(et2)*time*4*m.pi*area*1e4*(et2-et1)+bkmuon+bkatm+atmnu

###plotting
pl.loglog()
pl.xlabel('Energy Deposited, GeV')
pl.ylabel('Events over 1347 days')
pl.grid()
pl.scatter(et ,bins)
pl.step(et2,theory,label=r'background+$\gamma=2.718$, $\phi=2.655$')
pl.step(et2,theoryh,label=r'background+$\gamma=2$,$\phi=0.77$')
pl.step(et2,bkmuon,label=r'background atm $\mu$')
pl.step(et2,bkatm,label=r'background atm $\nu$')
pl.step(et2,atmnu,label=r'atm $\nu$ 90% CL Charm limit')
pl.xlim(10**4.2,1e7)
pl.ylim(10**-1.5,1e2)
pl.legend()
pl.show()
pl.subplot() #effective area
pl.loglog()
pl.step(e2,e3,label=r'$\nu_e$')
pl.step(m2,m3,label=r'$\nu_\mu$')
pl.step(t2,t3,label=r'$\nu_\tau$')
pl.step(et2,area,label='total')
pl.xlabel('Energy Deposited, GeV')
pl.ylabel(r'Effective area, $m^2$')
pl.legend(loc=4)
pl.show()

bkmuon2=bkmuon*(et2-et1)/(time*4*m.pi*area*1e4)
bkatm2=bkatm*(et2-et1)/(time*4*m.pi*area*1e4)
atmnu2=atmnu*(et2-et1)/(time*4*m.pi*area*1e4)
theory2=[]
theoryh2=Q(et2)*(et2-et1)**2+bkmuon2+bkatm2+atmnu2
obs=bins*(et2-et1)/(time*4*m.pi*area*1e4)

for i in range(0,len(et2)):
    if et2[i] <= 1e6:
        theory2.append(P1(et2[i])*(et2[i]-et1[i])**2+bkmuon2[i]+bkatm2[i]+atmnu2[i])
    if et2[i] > 1e6:
        theory2.append(P2(et2[i])*(et2[i]-et1[i])**2+bkmuon2[i]+bkatm2[i]+atmnu2[i])

pl.subplot() 
pl.loglog()
pl.xlabel('Energy Deposited, GeV')
pl.ylabel(r'$E^2\phi, GeV cm^{-1} s^{-1} sr^{-1}$')
pl.grid()
pl.scatter(et,obs)
pl.step(et2,theory2,label=r'background+$\gamma=2.718$, $\phi=2.655$')
pl.step(et2,theoryh2,label=r'background+$\gamma=2$,$\phi=0.77$')
pl.step(et2,bkmuon2,label=r'background atm $\mu$')
pl.step(et2,bkatm2,label=r'background atm $\nu$')
pl.step(et2,atmnu2,label=r'atm $\nu$ 90% CL Charm limit')
pl.xlim(10**4.2,1e7)
pl.ylim(1e-12,1e-6)
pl.legend()
pl.show()

#"""Computation of chi-2"""
#
##traditional chi squared
#kai=((bins-theory)*(bins-theory))/theory
#print("chi squared value")
#print(kai)
#print("sum",sum(kai),'14 energy bins')
#
#pvalue=[]
#for i in range(0,len(bins)):
#    if bins[i]==0:
#        chi=m.exp(-theory[i])
#    else:
#        chi=m.exp(bins[i]-theory[i])*((theory[i]/bins[i])**bins[i]) #as in papares
#    pvalue.append(chi)
#value=np.array(pvalue)
#twice=-2*np.log(value)
#print("p-value")
#print(twice)
#print("mean p", sum(twice),'14 bins')
#
#import scipy as sp
#Pchi=theory**bins*np.exp(-theory)/(sp.misc.factorial(bins))
#print("poissonian probabilty in each bin")
#print(Pchi)
#print("mean probability =",np.mean(Pchi))
