# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 21:46:17 2022

@author: Akshit
"""

#%% C, V \neq 0 

import numpy as np
import matplotlib.pyplot as plt

path = "sensitivity_plots_CVneq0/"

LT = ['-', '-', '-', '--', ':', '-.', '-.']
Mrk = ['s', 's', 'o', 'o', '+', '+', 'x']

leader3 = np.loadtxt(open(path+"C_Vneq0_insamp_obj.csv", "rb"), delimiter=",")
row, cloumn = np.shape(leader3)
gamma1_mat = np.loadtxt(open(path+"C_Vneq0_gamma1_mat.csv", "rb"), delimiter=",")
gamma2_mat = np.loadtxt(open(path+"C_Vneq0_gamma2_mat.csv", "rb"), delimiter=",")

row, column = np.shape(leader3)

f=plt.figure(figsize=(8,7)) 

index_set = np.arange(0, 37, 4)
for i in range(row):
    if gamma1_mat[i] in [0.2, 0.3, 0.4, 0.5]:
        plt.plot(gamma2_mat[index_set], -leader3[i,index_set], label='$\gamma_1=$'+str(np.round(gamma1_mat[i],2)), 
             linestyle = LT[min(i,5)], marker=Mrk[i])
   
plt.xlabel(r'$\gamma_2$')
plt.ylabel(r"$\mathrm{Total\ profit\ (\$)}$")
plt.legend()
f.savefig(path+"C_Vneq0_gamma2.pdf", bbox_inches='tight')
plt.close()

#%% C, V \neq 0 

import numpy as np
import matplotlib.pyplot as plt

path = "sensitivity_plots_CVneq0/"

leader3 = np.loadtxt(open(path+"C_Vneq0_gamma1_0.2_insamp_obj.csv", "rb"), delimiter=",")
Lb = np.loadtxt(open(path+"C_Vneq0_Lb_mat.csv", "rb"), delimiter=",")
row, column = np.shape(leader3)

gamma1 = 0.2
gamma2_mat = np.loadtxt(open(path+"C_Vneq0_gamma2_mat_lb.csv", "rb"), delimiter=",")

f=plt.figure(figsize=(9,8)) 

LT = ['-', '-', '--', ':', '-.']
Mrk = ['s', 's', '+', '^', 'x', 'd']

for i in range(column):
    if gamma2_mat[i] in [3, 5, 7]:
       plt.plot(Lb, -leader3[:,i], linestyle = LT[min(i,4)],
                label='$\gamma_1=$'+str(gamma1)+'$\ \gamma_2=$'+str(np.round(gamma2_mat[i],2)), marker=Mrk[i])

plt.xlabel(r'$\mathrm{Lower\ bound}$')
plt.ylabel(r"$\mathrm{Total\ profit\ (\$)}$")
plt.legend(loc='lower right')
plt.ylim(300, 700)
f.savefig(path+"C_Vneq0_lb.pdf", bbox_inches='tight')
plt.close()
