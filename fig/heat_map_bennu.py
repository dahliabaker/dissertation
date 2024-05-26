#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 08:44:28 2023

@author: dahliabaker
"""
#read plain yorp file

import numpy as np
import matplotlib.pyplot as plt

import asteroid_fortuna as asteroid
import MathOp as mo
import pandas as pd

v,f,af,Uf = asteroid.FileImport('Bennu_v20_degraded_6k.objp')
del af, Uf
#collect the vertices and connectivity lists for the original shape
verts = np.array(v)
facets = np.array(f)

#we care about the n hat vector of every facet
nhat = []
rvec = []
rvec_mag = []
#area_shapefacets = []
#for all facets, calculate nhat in body frame
for i in facets:
    v1 = np.array(verts[i[0]-1])
    v2 = np.array(verts[i[1]-1])
    v3 = np.array(verts[i[2]-1])
    u_vec = mo.hat(np.array(v2-v1))
    v_vec = mo.hat(np.array(v3-v1))
        
    n_dir = np.cross(u_vec,v_vec)
    n_dir = mo.hat(n_dir)
    nhat.append(n_dir)
    
    #also find center
    center_i = (v1+v2+v3)/3
    rvec.append(center_i)
    rvec_mag.append(np.linalg.norm(center_i[0:2]))
    #center coordinate corresponds to r vector in body frame

 
longlat_facet = np.zeros((len(facets),2))
 
 
for cnt in range(0,len(facets)):
    bigR = np.linalg.norm(rvec[cnt])
    longlat_facet[cnt,1] = np.rad2deg(np.arcsin(rvec[cnt][2]/bigR))#latitude
    longlat_facet[cnt,0] = np.rad2deg(np.arctan2(rvec[cnt][1],rvec[cnt][0]))#longitude

body_name = 'bennu'
face_num = len(facets)



A0i_db = np.zeros((face_num,3))
C0i_db = np.zeros((face_num,3))
A0i_body_file = open('body_file_'+body_name+'.txt','w+')
    

f1 = open('A0i_yorp_'+body_name,'r')
lines = f1.readlines()
A0i_x = []
A0i_y = []
A0i_z = []
facet_num = []
for item in lines[0:len(facets)]:
    sep = item.split()
    facet_num.append(float(sep[0]))
    A0i_x.append(float(sep[1]))
    A0i_y.append(float(sep[2]))
    A0i_z.append(float(sep[3]))

cnt=0
zipobj = zip(A0i_x,A0i_y,A0i_z)
for x,y,z in zipobj:
    A0i_db[cnt,:] = np.array([x,y,z])
    C0i_db[cnt,:] = np.cross(rvec[cnt],np.array([x,y,z]))
    cnt+=1   
    
#A0i_boulders = A0i[5898:,:]



for num in range(0,face_num):
    A0i_body_file.write(str(num)+" "+str(A0i_db[num,0])+ " "+str(A0i_db[num,1])+ " "+str(A0i_db[num,2])+ "\n")

               
  
A0i_body_file.close()



#%%
plt.figure(figsize=(9,6))
plt.hist(C0i_db[:,2],100)
plt.xlabel('YORP Coeff Value',fontsize=24)
plt.ylabel('Frequency',fontsize=24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.title('Base Shape Facet YORP Distribution',fontsize=30)



#%%
plt.figure(figsize=(9,6))
plt.scatter(range(0,len(facets)),C0i_db[:,2])
plt.xlabel('Facet Number',fontsize=24)
plt.ylabel('YORP Value',fontsize=24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.title('Base Shape Facet YORP Distribution',fontsize=30)
plt.text(-100,0.5e-5,'North Pole',fontsize=16)
plt.text(5100,-1e-5,'^ South Pole',fontsize=16)



#%%

plt.rcParams['figure.dpi']=500
plt.rcParams['xtick.labelsize']=24
plt.figure(figsize=(12,10))
#plt.scatter(longlat_facet[:,0],longlat_facet[:,1],c=np.sqrt(abs(C0i_db[:,2]))*np.sign(C0i_db[:,2]),cmap='coolwarm',vmin=-3e-3,vmax=3e-3,alpha=0.9)
plt.scatter(longlat_facet[:,0],longlat_facet[:,1],c=C0i_db[:,2],cmap='seismic',vmin=min(C0i_db[:,2]),vmax=-min(C0i_db[:,2]))
cb = plt.colorbar(orientation='horizontal')
cb.set_label(label="$\mathdefault{C_{0,z},(km^3)} $",fontsize=28)
plt.xlabel('Longitude',fontsize=24)
plt.ylabel('Latitude',fontsize=24)
plt.xticks((-180,-90,0,90,180),fontsize=24)
plt.yticks((-90,-45,0,45,90),fontsize=24)
plt.xlim([-180,180])
# plt.ylim([-90,90])


#%%
plt.figure(figsize=(12,10))
plt.scatter(longlat_facet[:,0],longlat_facet[:,1],c=rvec_mag,cmap='rainbow')
cb = plt.colorbar()

plt.xlabel('Longitude',fontsize=24)
plt.ylabel('Latitude',fontsize=24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.title('Facet Distance from $\hat{Z}$',fontsize=30)





