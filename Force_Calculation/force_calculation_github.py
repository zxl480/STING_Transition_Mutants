import sys
import numpy as np
import matplotlib.pyplot as plt

# Function to parse PDB atom lines
def pdbAtom(a):
    if a.startswith("END"):
        return 0
    atom = [a[12:16].strip(), a[17:20].strip(), int(a[22:26]), a[21],
            float(a[30:38]), float(a[38:46]), float(a[46:54])]
    if atom[3].strip() == '':
        atom[3] = None
    return list(atom)

# Function to calculate center of a set of atoms
def CENTER(x):
    C = np.array([atom[4:7] for atom in x])
    return np.mean(C, axis=0)

# Function to calculate distance between two points
def distance(p1, p2):
    t1 = np.array(p1)
    t2 = np.array(p2)
    vec = t1 - t2
    distx = vec[0]
    disty = vec[1]
    distz = vec[2]
    r = np.sqrt(distx * distx + disty * disty + distz * distz)
    return r, distx, disty, distz

# Function to calculate force between atom groups
def force(g1, g2):
    total_force = []
    for i in range(len(g1)):
        p1 = g1[i][4:7]
        for j in range(len(g2)):
            p2 = g2[j][4:7]
            dist = distance(p1, p2)
            forcex = ARG[i] * cGAMP[j] * 1.0 / pow(dist[0], 3) * dist[1]
            forcey = ARG[i] * cGAMP[j] * 1.0 / pow(dist[0], 3) * dist[2]
            forcez = ARG[i] * cGAMP[j] * 1.0 / pow(dist[0], 3) * dist[3]
            forcer = [forcex, forcey, forcez]
            total_force.append(forcer)
    return total_force

#charge of atoms
ARG=[-0.47,0.31,0.07,0.09,-0.18,0.09,0.09,-0.18,0.09,0.09,0.20,0.09,0.09,-0.70,0.44,0.64,-0.80,0.46,0.46,-0.80,0.46,0.46,0.51,-0.51]

cGAMP=[0.160000,0.090000,-0.500000,0.160000,0.090000,0.280000,-0.710000,0.340000,0.120000,-0.050000,-0.740000,0.500000,0.130000,-0.750000,0.430000,0.460000,-0.770000,
0.380000,0.380000,0.140000,0.090000,-0.660000,0.430000,0.010000,0.090000,-0.570000,-0.080000,0.090000,0.090000,-0.570000,1.500000,-0.780000,-0.780000,0.160000,
0.090000,-0.500000,0.160000,0.090000,-0.020000,0.260000,-0.680000,0.320000,0.350000,-0.740000,0.750000,-0.340000,0.260000,0.540000,-0.510000,0.000000,
-0.600000,0.250000,0.160000,0.010000,0.090000,-0.570000,0.140000,0.090000,-0.660000,0.430000,-0.080000,0.090000,0.090000,-0.570000,1.500000,-0.780000,-0.780000]

# Lists to store forces and positions
FF1, FF2 = [], []
PP1, PP2 = [], []

# Simulation loop for different frames
for frame in range(100):
    s = open("frame" + str(frame) + ".pdb", "r")
    title, atoms, box = [], [], []

    for i in s:
        # Parsing PDB data
        if i.startswith("ENDMDL"):
            print("Starting frame", frame)
        elif i.startswith("TITLE"):
            title.append(i)
        elif i.startswith("CRYST1"):
            print("Passing CRYST1")
        elif i.startswith("ATOM") or i.startswith("HETATM"):
            atoms.append(pdbAtom(i))

    atom1 = atoms[0:24]
    atom2 = atoms[24:48]
    atomref = atoms[48:115]

    # Calculating forces for atom groups
    Force1 = np.array(force(atom1, atomref))
    AF1 = np.mean(Force1, axis=0)
    FF1.append(AF1)
    Force2 = np.array(force(atom2, atomref))
    AF2 = np.mean(Force2, axis=0)
    FF2.append(AF2)

    # Calculating center positions
    C1 = CENTER(atom1)
    C2 = CENTER(atom2)
    Cref = CENTER(atomref)
    PP1.append(C1 - Cref)
    PP2.append(C2 - Cref)
    print("PP1:", PP1)

# Constants for unit conversion
Constant1 = 138.935458
Constant2 = 6.023

# Converting forces and positions to arrays
FF1 = np.array(FF1)
FF2 = np.array(FF2)
PP1 = np.array(PP1)
PP2 = np.array(PP2)

# Calculating average forces and positions
Ave_FF1 = np.mean(FF1, axis=0)
Ave_FF2 = np.mean(FF2, axis=0)
Ave_PP1 = np.mean(PP1, axis=0)
Ave_PP2 = np.mean(PP2, axis=0)

# Calculating tangential forces
XY_vec = Ave_PP1[0:2] - Ave_PP2[0:2]
XY_FF1 = Constant1 / Constant2 * 1000 * Ave_FF1[0:2]
par_force = np.dot(XY_vec, XY_FF1) / np.linalg.norm(XY_vec)
tan_force = np.sqrt(XY_FF1[0] * XY_FF1[0] + XY_FF1[1] * XY_FF1[1] - par_force * par_force)

XY_vec2 = Ave_PP2[0:2] - Ave_PP1[0:2]
XY_FF2 = Constant1 / Constant2 * 1000 * Ave_FF2[0:2]
par_force2 = np.dot(XY_vec2, XY_FF2) / np.linalg.norm(XY_vec2)
tan_force2 = np.sqrt(XY_FF2[0] * XY_FF2[0] + XY_FF2[1] * XY_FF2[1] - par_force2 * par_force2)

# Printing calculated tangential forces
print("Tangential forces:", tan_force, tan_force2)
