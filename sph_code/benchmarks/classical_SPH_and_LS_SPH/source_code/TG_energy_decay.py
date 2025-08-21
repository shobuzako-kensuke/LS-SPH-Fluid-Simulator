#=========================#
#  module                 # 
#=========================#
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import ana_read_fortran_files as read_f90
plt.rcParams['mathtext.fontset'] = 'cm' # mathfont for figure

#=========================#
#  main                   # 
#=========================#
def main(Re, N_sys, N_inn, N_file, last_file, write_step, Delta_t_eff, parent_path, save_path):
    count  = 0
    time   = np.zeros((N_file))  # for time record
    energy = np.zeros((N_file)) # for energy record

    #=========================#
    #  read & cal             # 
    #=========================#
    for i in range(0, last_file+1, write_step):
        time[count] = Delta_t_eff * float(i)
        
        file_name = parent_path + '/data/{}/SP_uvw.dat'.format(i)  # read SP_uvw
        SP_uvw    = read_f90.binary_files(file_name, N_sys, 2)
        energy_tmp = 0.0
        for me in range(0, N_inn):
            energy_tmp += SP_uvw[me,0]**2.0 + SP_uvw[me,1]**2.0

        energy[count] = 0.5* energy_tmp / float(N_inn)

        count += 1

    #=========================#
    #  figure                 # 
    #=========================#
    ana = np.zeros((N_file))
    ana = 0.25* np.exp(-4* math.pi**2* time[:]/ Re)

    fig = plt.figure(figsize=(6.7,6), facecolor='white')
    plt.subplots_adjust(left=0.25, right=0.95, bottom=0.18, top=0.95, wspace=0.35, hspace=0.4)
    ax = fig.add_subplot(111, facecolor='white')
    # plot
    ax.plot(time, energy, c='r', label='result', linewidth=1.5)
    ax.plot(time, ana, c='k', label='Analytical', linewidth=1.5)
    # ax.set_yscale('log')
    # axis label
    ax.set_xlabel('$t$', fontsize=30, labelpad=10)
    ax.set_ylabel('$E_{k}$', fontsize=30, labelpad=16)
    # ticks
    ax.tick_params(axis='both', which='major', direction='out', length=6, width=1, labelsize=22)
    ax.minorticks_off() # minor ticks
    # grid
    ax.grid(which='major', color='none', linewidth=0.5) # None
    ax.set_axisbelow(True) # grid back
    # outer frame
    ax.spines["bottom"].set_linewidth(1.2)
    ax.spines["top"].set_linewidth(1.2)
    ax.spines["right"].set_linewidth(1.2)
    ax.spines["left"].set_linewidth(1.2)
    # save
    fig.savefig('{}/energy_decay.png'.format(save_path), format='png', dpi=300, transparent=False)
    plt.close()

    # save
    np.savez(parent_path + '/ana_energy_decay.npz', t=time, E=energy)  # save

# END #