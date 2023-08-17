import sys
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

dcd_file_with = "xkr4_with"
dcd_file_without = "xkr4_without"
filename = "Supfig6"

def make_corarray(x):
    dcd_file = str(x)
    u = mda.Universe('%s.dcd' %dcd_file)
    
    tra = u.trajectory

    time = len(tra)
    num_atoms = len(u.atoms)

    #The first three are the primary binding residues
    residues = ["E123", "D127", "E310"]
    
    coor = {}
    for i in residues:
        coor[i] = np.zeros((time, 3))
    for i in residues:
        for j in range(time):
            num = i[1:]
            num = int(num)
            coor[i][j] = tra[j][num]
    #calculate the center of the primary binding amino acids
    center = np.zeros((time, 3))
    for i in range(time):
        center[i] = (coor[residues[0]][i]+coor[residues[1]][i]+coor[residues[2]][i])/3

    #calculate the distances
    dis = np.zeros((num_atoms, time))
    for i in range(num_atoms):
        for j in range(time):
            distance = np.sqrt((tra[j][i][0]-center[j][0])**2 + (tra[j][i][1]-center[j][1])**2 + (tra[j][i][2]-center[j][2])**2)
            dis[i][j] = distance
    output = np.corrcoef(dis)
    return output

coef_with = make_corarray(dcd_file_with)
coef_without = make_corarray(dcd_file_without)

fig = plt.figure(figsize = (30,12))

ax = fig.add_subplot(1, 2, 1)
sns.heatmap(coef_without, ax=ax, vmin = -1.0, vmax = 1.0)

ax = fig.add_subplot(1, 2, 2)
sns.heatmap(coef_with, ax=ax, vmin = -1.0, vmax = 1.0)

plt.savefig("%s_matrix.pdf" %filename)
plt.show()
