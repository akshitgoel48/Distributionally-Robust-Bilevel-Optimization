# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 21:59:40 2022

@author: Akshit
"""

#%% C = V = 0
import numpy as np
import matplotlib.pyplot as plt

path = "sensitivity_plots_CVeq0/"

LT = ['-', '-', '-', '--', ':', '-.']
Mrk = ['s', 's', 's', 'd', '+', 'x']

leader3 = np.loadtxt(open(path+"C_Veq0_insamp_obj.csv", "rb"), delimiter=",")
row, cloumn = np.shape(leader3)
gamma1_mat = np.loadtxt(open(path+"C_Veq0_gamma1_mat.csv", "rb"), delimiter=",")
gamma2_mat = np.loadtxt(open(path+"C_Veq0_gamma2_mat.csv", "rb"), delimiter=",")

row, column = np.shape(leader3)

f=plt.figure(figsize=(8,7)) 
 
for i in range(row):
    plt.plot(gamma2_mat, -leader3[i,:], linestyle = LT[i], marker=Mrk[i],
             label='$\gamma_1=$'+str(np.round(gamma1_mat[i],2)))

plt.xlabel(r'$\gamma_2$')
plt.ylabel(r"$\mathrm{Total\ profit\ (\$)}$")
plt.legend(loc='lower left')
f.savefig(path+"C_Veq0_gamma2.pdf", bbox_inches='tight')
plt.close()

#%% C = V = 0

path = "sensitivity_plots_CVeq0/"

import numpy as np
import matplotlib.pyplot as plt

leader3 = np.loadtxt(open(path+"C_Veq0_gamma1_1.3_insamp_obj.csv", "rb"), delimiter=",")
Lb = np.loadtxt(open(path+"C_Veq0_Lb_mat.csv", "rb"), delimiter=",")
row, column = np.shape(leader3)

gamma1 = 1.3
gamma2_mat = np.loadtxt(open(path+"C_Veq0_gamma2_mat_lb.csv", "rb"), delimiter=",")

row, column = np.shape(leader3)

f = plt.figure(figsize=(8,6)) 

LT = ['-', '-', '--', ':', '-.']
Mrk = ['s', 's', 'd', '+', 'x']

for i in range(column):
    if gamma2_mat[i] != 2:
        plt.plot(Lb, -leader3[:,i], linestyle = LT[i],
                 label='$\gamma_1=$'+str(gamma1)+'$\ \gamma_2=$'+str(np.round(gamma2_mat[i],2)), marker=Mrk[i])

plt.xlabel(r'$\mathrm{Lower\ bound}$')
plt.ylabel(r"$\mathrm{Total\ profit\ (\$)}$")
plt.legend(loc='upper left')
f.savefig(path+"C_Veq0_lb.pdf", bbox_inches='tight')
plt.close()
