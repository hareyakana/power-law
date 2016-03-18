# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 13:42:18 2015

@author: hareyakana
BROKEN POWER LAW
2POWER LAW CUT OFF AT ... GeV
"""
import matplotlib.pylab as pl
from matplotlib import mlab,cm
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
bkmuoni=np.array([0,0.1,1.1,0.6,0.33,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10])
bkmuon=10**bkmuoni
bkatmi=np.array([0,0,0.35,0.345,0.25,-0.1,-0.6,-1.1,-1.9,-10,-10,-10,-10,-10,-10])
bkatm=10**bkatmi
atm=np.array([0,0.15,0.62,0.75,0.72,0.48,0.2,-0.3,-0.65,-1.1,-1.52,-10,-10,-10,-10])
atmnu=10**atm
background=bkmuon+bkatm+atmnu
#####
ET1=np.linspace(4.0,6.8,num=15)
ET2=np.linspace(4.2,7.0,num=15)
et1=10**ET1
et2=10**ET2
et=(et1+et2)/2


bins=np.zeros(len(ET1))  ###DATA

for i in range(len(E)):
    for k in range(len(ET1)):
        if abs(e[i])>abs(et1[k]) and abs(e[i])<abs(et2[k]):
            bins[k]=bins[k]+1
    
"""Theory model/etc"""
time=1347*24*60*60  #total IceCube sample time

"""2 POWER LAW/broken power law"""
def P(E,n,g,e,a,bk):  # of unit s-1 cm-2
    p=n*1e-18*(((E)/1e5)**-g)
    f1=p*time*4*m.pi*a*1e4*(E-e)+bk
    return f1
    
def P2(E,n,g,e,a,bk):
    p=n*1e-18*(((E)/1e5)**-g)
    f2=p*time*4*m.pi*a*1e4*(E-e)+bk
    return f2

"""Effective area"""
"""from public data"""
import csv
e1=[]
e2=[]
e3=[]

with open('nue_4pi.txt','r')as f:
    reader=csv.DictReader(f)
    for row in reader:
        e1.append(float(row['x']))
        e2.append(float(row['y']))
        e3.append(float(row['z']))

m1=[]
m2=[]
m3=[]

with open('numu_4pi.txt','r')as f:
    reader=csv.DictReader(f)
    for row in reader:
        m1.append(float(row['x']))
        m2.append(float(row['y']))
        m3.append(float(row['z']))        

t1=[]
t2=[]
t3=[]

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

for x in range(0,round(len(area)/2-.9)):
    area[2*x]=area[2*x]*3/5.

for x in range(0,round(len(area)/2)):
    area[x*2-1]=area[x*2-1]*2/5.
    
#########################
"""best fit parameter"""
#bins[-4]=bins[-5]=0

#########################
#########################
phi=np.linspace(0,5,num=21)
phi2=np.linspace(0,30,num=21)
gamma=np.linspace(2,4,num=21)
gamma2=np.linspace(2,4,num=21)

for i in range(0,len(bins)):
    if bins[i]==0:
        bins[i]=1e-20

p,g,g2=np.meshgrid(phi,gamma,gamma2)
#import scipy as sp
def chi(observed,expected):
#    ch=(observed-expected)*(observed-expected)/expected
#    ch=(expected**observed)*np.exp(-expected)/(sp.misc.factorial(observed)) #poisson probabilty
    ch=(np.exp(observed-expected))*((expected/observed)**observed)
#    ch=-2*np.log(ch1)    
    return ch

blank=[]
#for i in range(0,len(et)):
#    theo=P(et2[i],p,g,et1[i],area[i],background[i],g2)
#    theory=chi(bins[i],theo)
#    blank.append(theory)
#
#total=np.zeros(shape=(len(phi),len(gamma),len(gamma2)))
#
#for i in range(0,len(blank)):
#    total+=blank[i]
    
#import matplotlib.gridspec as gridspec

#h, w = total.shape
#gs = gridspec.GridSpec(2,2,width_ratios=[w,w*.1], height_ratios=[h,h*.1])
#ax = [pl.subplot(gs[0]),]
#ax.append(pl.subplot(gs[1], sharey=ax[0]))
#ax.append(pl.subplot(gs[2], sharex=ax[0]))
#ax.append(pl.subplot(gs[3]))
#bounds = [phi.min(),phi.max(),gamma.min(),gamma.max()]
#ax[0].imshow(total, cmap='Blues_r', extent=bounds,origin='lower')
#ax[1].plot(total[:,int(w/2)],g[:,int(w/2)],'.',total[:,int(w/2)],g[:,int(w/2)])
#ax[1].axis([total[:,int(w/2)].max(), total[:,int(w/2)].min(), g.min(), g.max()])
#ax[2].plot(p[int(h/2),:],total[int(h/2),:],'.',p[int(h/2),:],total[int(h/2),:])

p1=[]
p2=[]
x=[]
y=[]
z=[]
w=[]
test=[]

#cutoff at 1PeV
for p in range(0,len(phi)):
    for p2 in range(0,len(phi2)):
        for q in range(0,len(gamma)):
            for q2 in range(0,len(gamma2)):
                x.append(phi[p])
                y.append(phi2[p2])
                z.append(gamma[q])
                w.append(gamma2[q2])
                for i in range(0,len(et2)):
                    if et2[i] <= 1e6:
                        th11=P(et2[i],phi[p],gamma[q],et1[i],area[i],background[i])
                        th1=chi(bins[i],th11)                        
                        test.append(th1)
                    if et2[i] > 1e6:
                        th22=P(et2[i],phi2[p2],gamma2[q2],et1[i],area[i],background[i])
                        th2=chi(bins[i],th22)                        
                        test.append(th2)

H=len(phi)*len(phi2)*len(gamma)*len(gamma2)
stat2=np.reshape(test,(H,-1))
stat=np.sum(stat2,axis=1)
print(stat,min(stat),max(stat))
print(x[np.argmax(stat)],y[np.argmax(stat)],
        z[np.argmax(stat)],w[np.argmax(stat)])

#for i in range(0,len(phi)):
#    for k in range(0,len(gamma)):
#        for a in range(0,len(gamma2)):
#            x.append(phi[i])
#            y.append(gamma[k])
#            z.append(gamma2[a])
#            for b in range(0,len(et2)):
#                theory=P(et2[b],phi[i],gamma[k],et1[b],area[b],background[b])
#                chisquared=chi(bins[b],theory)
#                p1.append(chisquared)
#                theoryw=P(et2[b],phi[i],gamma2[a],et1[b],area[b],background[b])
#                chisquaredw=chi(bins[b],theoryw)
#                p2.append(chisquaredw)
#                test.append(et2[b])
#    
#stat2=[]
#for i in range(0,len(p1)):
#    if test[i]<10**6:
#        stat2.append(p1[i])
#    else:
#        stat2.append(p2[i])
#
#shapeL=int(len(stat2)/len(et2))
#shapeH=len(et2)
#stat3=np.reshape(stat2,(shapeL, -1))
#stat=np.sum(stat3,axis=1)
#stat=np.array(stat)

#print(stat3)
#print(len(stat),len(stat3))
    
#ma=max(p1)/max(p1)
#mi=min(p1)/max(p1)
#levels=np.arange(mi,ma,(ma-mi)/3)
#z=total/max(p1)
#cax=ax[0].contour(z,levels,cmap='viridis',extent=bounds,origin='lower')
#cbar=pl.colorbar(cax,orientation='vertical')
#pl.tight_layout()
#pl.show()


#dummy=[]
#for i in range(0,len(phi)):
#    g=2.
#    th=P(et2,phi[i],g,et1,area,background)
#    ch=chi(bins,th)
#    k=sum(ch)
#    dummy.append(k)

#print("maximum phi, gamma for chi squared ")
#print(min(stat),max(stat),len(stat))
#print(np.argmax(stat),
#      'phi',x[np.argmax(stat)],
#      'gamma<1PeV',y[np.argmax(stat)],
#      'gamma>1PeV',z[np.argmax(stat)])

#print("gamma =2, max phi=",phi[dummy.index(max(dummy))])

