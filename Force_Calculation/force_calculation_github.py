import sys, logging, random, math, os, re
import numpy as np 
import matplotlib.pyplot as plt

def pdbAtom(a):
    if a.startswith("END"):
        return 0
    atom = [a[12:16].strip(), a[17:20].strip(), int(a[22:26]), a[21],
            float(a[30:38]), float(a[38:46]), float(a[46:54])]
    if atom[3].strip() == '':
        atom[3] = None
    return list(atom) 

#charge of atoms
ARG=[-0.47,0.31,0.07,0.09,-0.18,0.09,0.09,-0.18,0.09,0.09,0.20,0.09,0.09,-0.70,0.44,0.64,-0.80,0.46,0.46,-0.80,0.46,0.46,0.51,-0.51]

cGAMP=[0.160000,0.090000,-0.500000,0.160000,0.090000,0.280000,-0.710000,0.340000,0.120000,-0.050000,-0.740000,0.500000,0.130000,-0.750000,0.430000,0.460000,-0.770000,
0.380000,0.380000,0.140000,0.090000,-0.660000,0.430000,0.010000,0.090000,-0.570000,-0.080000,0.090000,0.090000,-0.570000,1.500000,-0.780000,-0.780000,0.160000,
0.090000,-0.500000,0.160000,0.090000,-0.020000,0.260000,-0.680000,0.320000,0.350000,-0.740000,0.750000,-0.340000,0.260000,0.540000,-0.510000,0.000000,
-0.600000,0.250000,0.160000,0.010000,0.090000,-0.570000,0.140000,0.090000,-0.660000,0.430000,-0.080000,0.090000,0.090000,-0.570000,1.500000,-0.780000,-0.780000]


def CENTER(x):
    C=[]
    for i in range (0, len(x)):
        C.append(x[i][4:7]);
    C=np.array(C)
    center=[np.average(C[:,0]),np.average(C[:,1]),np.average(C[:,2])];
    return center

def distance(p1,p2):
    t1=np.array(p1);t2=np.array(p2)
    vec=t1-t2;
    distx=vec[0];disty=vec[1];distz=vec[2];
    r=pow(distx*distx+disty*disty+distz*distz,0.5);
    return r,distx,disty,distz

def force(g1,g2):
    Total_force=[]
    for i in range (0, len(g1)):
        p1=g1[i][4:7];
        for j in range (0,len(g2)):
            p2=g2[j][4:7]
            dist=distance(p1,p2); 
            forcex=ARG[i]*cGAMP[j]*1.0/pow(dist[0],3)*dist[1]; forcey=ARG[i]*cGAMP[j]*1.0/pow(dist[0],3)*dist[2]; forcez=ARG[i]*cGAMP[j]*1.0/pow(dist[0],3)*dist[3]
            forcer=[forcex,forcey,forcez];
            Total_force.append(forcer)
    return Total_force        


#calculate the force between cGAMP and ARG1 / ARG2

FF1,FF2=[],[]; PP1,PP2=[],[] #force and position
for frame in range (0, 100):    
    s = open ("frame"+str(frame)+".pdb","r")
    title, atoms, box = [], [], []
    for i in s:
        if i.startswith("ENDMDL"):
            print ("starting")
        elif i.startswith("TITLE"):
            title.append(i)
        elif i.startswith("CRYST1"):
             print ("pass")
        elif i.startswith("ATOM") or i.startswith("HETATM"):
            atoms.append(pdbAtom(i));

    atom1=atoms[0:24]
    atom2=atoms[24:48]
    atomref=atoms[48:115]; 
          

    Force1=np.array(force(atom1,atomref)); 
    AF1=[np.average(Force1[:,0]),np.average(Force1[:,1]),np.average(Force1[:,2])]; 
    FF1.append(AF1);
    Force2=np.array(force(atom2,atomref)); 
    AF2=[np.average(Force2[:,0]),np.average(Force2[:,1]),np.average(Force2[:,2])];
    FF2.append(AF2);
    
    C1=np.array(CENTER(atom1));C2=np.array(CENTER(atom2));Cref=np.array(CENTER(atomref));    
    PP1.append(C1-Cref);
    PP2.append(C2-Cref); print (PP1)

#units coversion    
Constant1=138.935458
Constant2=6.023

FF1=np.array(FF1);FF2=np.array(FF2);PP1=np.array(PP1);PP2=np.array(PP2);
Ave_FF1=[np.average(FF1[:,0]),np.average(FF1[:,1]),np.average(FF1[:,2])]; 
Ave_FF2=[np.average(FF2[:,0]),np.average(FF2[:,1]),np.average(FF2[:,2])]; 
Ave_PP1=[np.average(PP1[:,0]),np.average(PP1[:,1]),np.average(PP1[:,2])]; 
Ave_PP2=[np.average(PP2[:,0]),np.average(PP2[:,1]),np.average(PP2[:,2])]; 

#calculate tangentional force

XY_vec=np.array(Ave_PP1[0:2])-np.array(Ave_PP2[0:2])
XY_FF1=Constant1/Constant2*1000*np.array(Ave_FF1[0:2]); 
par_force=np.dot(XY_vec,XY_FF1)/np.linalg.norm(XY_vec)
tan_force=pow(XY_FF1[0]*XY_FF1[0]+XY_FF1[1]*XY_FF1[1]-par_force*par_force,0.5);

XY_vec2=np.array(Ave_PP2[0:2])-np.array(Ave_PP1[0:2])
XY_FF2=Constant1/Constant2*1000*np.array(Ave_FF2[0:2]); 
par_force2=np.dot(XY_vec2,XY_FF2)/np.linalg.norm(XY_vec2); 
tan_force2=pow(XY_FF2[0]*XY_FF2[0]+XY_FF2[1]*XY_FF2[1]-par_force2*par_force2,0.5); 
print (tan_force,tan_force2)


