import pandas as pd
from general_plot import general_plot
import numpy as np
import matplotlib.pyplot as plt

if __name__=="__main__":

    timesteps1 = np.arange(1,42,1)
    
    fig,ax = plt.subplots(1,1)
    
    a = 4
    if a == 1:
        error_array1 = np.loadtxt('error_local3.txt')
        ax.plot(timesteps1, error_array1, '-o', label='Local radius 3')
        #ax.plot(timesteps1, error_array1, '-ro', label='Local radius 5')

        error_array1 = np.loadtxt('error_local5.txt')
        ax.plot(timesteps1, error_array1,  '-o', label='Local radius 5')
        #ax.plot(timesteps1, error_array1, '-ro', label='Local radius 5')

        error_array1 = np.loadtxt('error_local10.txt')
        #ax.plot(timesteps1, error_array1, '-go', label='Local radius 10')
        ax.plot(timesteps1, error_array1,  '-o', label='Local radius 10')

        error_array1 = np.loadtxt('error_local15.txt')
        #ax.plot(timesteps1, error_array1, '-bo', label='Local radius 15')
        ax.plot(timesteps1, error_array1,  '-o', label='Local radius 15')

        error_array1 = np.loadtxt('error_local20.txt')
        #ax.plot(timesteps1, error_array1, '-co', label='Local radius 20')
        ax.plot(timesteps1, error_array1,  '-o', label='Local radius 20')
    
        error_array1 = np.loadtxt('error_local25.txt')
        #ax.plot(timesteps1, error_array1, '-co', label='Local radius 25')
        ax.plot(timesteps1, error_array1,  '-o', label='Local radius 25')
    
        error_array1 = np.loadtxt('error_local35.txt')
        #ax.plot(timesteps1, error_array1, '-co', label='Local radius 25')
        ax.plot(timesteps1, error_array1,  '-o', label='Local radius 35')
        
        error_array1 = np.loadtxt('error_global.txt')
        #ax.plot(timesteps1, error_array1, '-co', label='Local radius 25')
        ax.plot(timesteps1, error_array1,  '-o', label='Global')
    elif a==2:
        error_array1 = np.loadtxt('reduced_error_local3.txt')
        ax.plot(timesteps1, error_array1, '-o', label='Local radius 3')
        #ax.plot(timesteps1, error_array1, '-ro', label='Local radius 5')

        error_array1 = np.loadtxt('reduced_error_local5.txt')
        ax.plot(timesteps1, error_array1,  '-o', label='Local radius 5')
        #ax.plot(timesteps1, error_array1, '-ro', label='Local radius 5')

        error_array1 = np.loadtxt('reduced_error_local10.txt')
        #ax.plot(timesteps1, error_array1, '-go', label='Local radius 10')
        ax.plot(timesteps1, error_array1,  '-o', label='Local radius 10')

        error_array1 = np.loadtxt('reduced_error_local15.txt')
        #ax.plot(timesteps1, error_array1, '-bo', label='Local radius 15')
        ax.plot(timesteps1, error_array1,  '-o', label='Local radius 15')

        error_array1 = np.loadtxt('reduced_error_local20.txt')
        #ax.plot(timesteps1, error_array1, '-co', label='Local radius 20')
        ax.plot(timesteps1, error_array1,  '-o', label='Local radius 20')
    
        error_array1 = np.loadtxt('reduced_error_local25.txt')
        #ax.plot(timesteps1, error_array1, '-co', label='Local radius 25')
        ax.plot(timesteps1, error_array1,  '-o', label='Local radius 25')
    
        error_array1 = np.loadtxt('reduced_error_local35.txt')
        #ax.plot(timesteps1, error_array1, '-co', label='Local radius 25')
        ax.plot(timesteps1, error_array1,  '-o', label='Local radius 35')
    
        error_array1 = np.loadtxt('reduced_error_global.txt')
        #ax.plot(timesteps1, error_array1, '-co', label='Local radius 25')
        ax.plot(timesteps1, error_array1,  '-o', label='Global')
    elif a==3:
        error_array1 = np.loadtxt('reduced_error_local20_ens64.txt')
        #ax.plot(timesteps1, error_array1, '-co', label='Local radius 20')
        ax.plot(timesteps1, error_array1,  '-o', label='Ens 64; Local radius 20')
        
        error_array1 = np.loadtxt('reduced_error_local20_ens16.txt')
        #ax.plot(timesteps1, error_array1, '-co', label='Local radius 20')
        ax.plot(timesteps1, error_array1,  '-o', label='Ens 16; Local radius 20')
    
    elif a==4:
        error_array1 = np.loadtxt('../Assimilated_results/localization_effect/4_obs/reduced_error_local25.txt')
        #ax.plot(timesteps1, error_array1, '-co', label='Local radius 20')
        ax.plot(timesteps1, error_array1,  '-o', label='4 obs Local radius 25')
        
        error_array1 = np.loadtxt('../Assimilated_results/localization_effect/12_obs/reduced_error_local15.txt')
        #ax.plot(timesteps1, error_array1, '-co', label='Local radius 20')
        ax.plot(timesteps1, error_array1,  '-o', label='12 obs Local radius 15')
        
        error_array1 = np.loadtxt('../Assimilated_results/localization_effect/35_obs/reduced_error_local15.txt')
        #ax.plot(timesteps1, error_array1, '-co', label='Local radius 20')
        ax.plot(timesteps1, error_array1,  '-o', label='35 obs Local radius 15')
        
        error_array1 = np.loadtxt('../Assimilated_results/localization_effect/42_obs/reduced_error_local15.txt')
        #ax.plot(timesteps1, error_array1, '-co', label='Local radius 20')
        ax.plot(timesteps1, error_array1,  '-o', label='42 obs Local radius 15')

    
    error_array1 = np.loadtxt('error_freerun.txt')
    #ax.plot(timesteps1, error_array1, '-mo', label='Free run')
    ax.plot(timesteps1, error_array1,  '-o', label='Free run')

    ax.set_xlabel('Assimilation step')
    ax.set_ylabel('Error norm')
    plt.title('Error comparison for twin experiment')
    plt.legend(loc='best')
    plt.savefig('error_check_ensALL.pdf')
    print "Plotted error norm plot"
    
    
    if 0:
        fig,ax = plt.subplots(1,1)

        error_array1 = np.loadtxt('error_32ens.txt')
        ax.plot(timesteps1, error_array1, '-ro', label='ESTKF')
   
        error_array1 = np.loadtxt('error_32ens_enkf.txt')
        ax.plot(timesteps1, error_array1, '-go', label='EnKF')

        error_array1 = np.loadtxt('error_freerun.txt')
        ax.plot(timesteps1, error_array1, '-mo', label='Free run')

        ax.set_xlabel('Assimilation step')
        ax.set_ylabel('Error norm')
        plt.title('Error comparison for twin experiment')
        plt.legend(loc='best')
        plt.savefig('error_filter_comp.pdf')
        print "Plotted error_filter_comp.pdf"
