# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 21:50:11 2022

@author: Akshit
"""

#%% Boxplots (NEW)

import numpy as np
import matplotlib.pyplot as plt

def set_box_color(bp, color, color1):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color1)
    plt.setp(bp['means'], color='yellow')

def boxplot_smp(out_TSBP, out_DRBP1, out_DRBP2, title_plot, choice, ylim_vec=[]):
    data_1 = np.concatenate((out_TSBP[:,[0]], out_TSBP[:,[1]], out_TSBP[:,[2]]), axis = 1)
    data_2 = np.concatenate((out_DRBP1[:,[0]], out_DRBP1[:,[1]], out_DRBP1[:,[2]]), axis = 1)
    data_3 = np.concatenate((out_DRBP2[:,[0]], out_DRBP2[:,[1]], out_DRBP2[:,[2]]), axis = 1)     
    fig = plt.figure(figsize =(8, 4)) 
    bp1 = plt.boxplot(-data_1, positions=[1, 5, 9], widths=0.6, showmeans=choice, showfliers=True, whis=(5,95), 
                      flierprops={'marker': '+', 'markeredgecolor': 'lime'}, 
                      meanprops={'markerfacecolor': 'brown', 'markeredgecolor': 'brown'})
    bp2 = plt.boxplot(-data_2, positions=[2, 6, 10], widths=0.6, showmeans=choice, showfliers=True, whis=(5,95), 
                      flierprops={'marker': '+', 'markeredgecolor': 'lime'},
                      meanprops={'markerfacecolor': 'brown', 'markeredgecolor': 'brown'})
    bp3 = plt.boxplot(-data_3, positions=[3, 7, 11], widths=0.6, showmeans=choice, showfliers=True, whis=(5,95), 
                      flierprops={'marker': '+', 'markeredgecolor': 'lime'},
                      meanprops={'markerfacecolor': 'brown', 'markeredgecolor': 'brown'})
    color1 = 'black'
    set_box_color(bp1, 'red', color1) # colors are from http://colorbrewer2.org/    
    set_box_color(bp2, 'blue', color1)
    set_box_color(bp3, 'gold', color1)
    plt.plot([], c='red', label='TS-BP')
    plt.plot([], c='blue', label='DRBP($\gamma_1=0,\gamma_2=1$)')
    plt.plot([], c='gold', label='DRBP($\gamma_1=0.5,\gamma_2=1$)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks([2, 6, 10], ticks)
    plt.tight_layout()
    plt.xlabel(r'$\mathrm{Sample\ size\ }(N)$')
    plt.ylabel(r"$\mathrm{Total\ profit\ (\$)}$")
    if np.shape(ylim_vec)[0] > 0:
        plt.ylim(ylim_vec[0],ylim_vec[1])
    fig.savefig(path+title_plot+".pdf", bbox_inches='tight')  # saves figure to pdf file. 
    plt.close()
 

ticks = ['$N=10$', '$N=100$', '$N=1000$']

#%  C=V=0 : set path to C_Veq0_out_samp_obj
path = "Boxplots/C_Veq0_out_samp_obj/"

N1mat = np.loadtxt(open(path+"N1mat.csv", "rb"), delimiter=",")
BOX_TSBP_uniform = np.loadtxt(open(path+"BOX_TSBP_uniform.csv", "rb"), delimiter=",")
BOX_TSBP_uniform2 = np.loadtxt(open(path+"BOX_TSBP_uniform2.csv", "rb"), delimiter=",")
BOX_TSBP_normal = np.loadtxt(open(path+"BOX_TSBP_normal.csv", "rb"), delimiter=",")

BOX_DRBP1_uniform = np.loadtxt(open(path+"BOX_DRBP1_uniform.csv", "rb"), delimiter=",")
BOX_DRBP1_uniform2 = np.loadtxt(open(path+"BOX_DRBP1_uniform2.csv", "rb"), delimiter=",")
BOX_DRBP1_normal = np.loadtxt(open(path+"BOX_DRBP1_normal.csv", "rb"), delimiter=",")

BOX_DRBP2_uniform = np.loadtxt(open(path+"BOX_DRBP2_uniform.csv", "rb"), delimiter=",")
BOX_DRBP2_uniform2 = np.loadtxt(open(path+"BOX_DRBP2_uniform2.csv", "rb"), delimiter=",")
BOX_DRBP2_normal = np.loadtxt(open(path+"BOX_DRBP2_normal.csv", "rb"), delimiter=",")

boxplot_smp(BOX_TSBP_uniform, BOX_DRBP1_uniform, BOX_DRBP2_uniform, 'Uniform True', True, [])
boxplot_smp(BOX_TSBP_uniform2, BOX_DRBP1_uniform2, BOX_DRBP2_uniform2, 'Uniform Misspecified', True, [])
boxplot_smp(BOX_TSBP_normal, BOX_DRBP1_normal, BOX_DRBP2_normal, 'Normal Misspecified', True, [])

#%  C,Vneq0 : set path to C_Vneq0_out_samp_obj
path = "Boxplots/C_Vneq0_out_samp_obj/"

N1mat = np.loadtxt(open(path+"N1mat.csv", "rb"), delimiter=",")
BOX_TSBP_uniform = np.loadtxt(open(path+"BOX_TSBP_uniform.csv", "rb"), delimiter=",")
BOX_TSBP_uniform2 = np.loadtxt(open(path+"BOX_TSBP_uniform2.csv", "rb"), delimiter=",")
BOX_TSBP_normal = np.loadtxt(open(path+"BOX_TSBP_normal.csv", "rb"), delimiter=",")

BOX_DRBP1_uniform = np.loadtxt(open(path+"BOX_DRBP1_uniform.csv", "rb"), delimiter=",")
BOX_DRBP1_uniform2 = np.loadtxt(open(path+"BOX_DRBP1_uniform2.csv", "rb"), delimiter=",")
BOX_DRBP1_normal = np.loadtxt(open(path+"BOX_DRBP1_normal.csv", "rb"), delimiter=",")

BOX_DRBP2_uniform = np.loadtxt(open(path+"BOX_DRBP2_uniform.csv", "rb"), delimiter=",")
BOX_DRBP2_uniform2 = np.loadtxt(open(path+"BOX_DRBP2_uniform2.csv", "rb"), delimiter=",")
BOX_DRBP2_normal = np.loadtxt(open(path+"BOX_DRBP2_normal.csv", "rb"), delimiter=",")

boxplot_smp(BOX_TSBP_uniform, BOX_DRBP1_uniform, BOX_DRBP2_uniform, 'Uniform True', True, [850, 2000])
boxplot_smp(BOX_TSBP_uniform2, BOX_DRBP1_uniform2, BOX_DRBP2_uniform2, 'Uniform Misspecified', True, [600, 1650])
boxplot_smp(BOX_TSBP_normal, BOX_DRBP1_normal, BOX_DRBP2_normal, 'Normal Misspecified', True, [850, 2000])