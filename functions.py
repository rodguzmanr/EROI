#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions needed by the EROI_simple_model.py script to compute
   and plot the main features of the indealized EROI scenarios.

This file contains the mathematical formulas needed in order to
compute the EROI and ERDE variables corresponding to the model
hypothesis and the values set in the variables.py file.
It also contains the functions to plot and save graphically the
results from the simulation.

November 2021, v1 by Rodrigo Guzman

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def set_df_time(var_list, time_length):
    """Function to initialize the time variable considered for
       the simulation for each variable.

    Parameters
    ----------   
    var_list : list (Python object)
        List with all the dataframes needing the time variable.

    time_length : int
        Number of time steps for the time vector.

    Returns
    -------
    var_list : list (Python object)
        List with all the updated dataframes (time variable
        has been added).

    """

    # Defining the time vector
    time_vect = np.arange(time_length)

    # Setting the time variable for each dataframe
    for var in var_list:
        var['time'] = time_vect
#    energy_res = energy_res.set_index('time')

    return var_list


def compute_eroi_and_net(energy_res):
    """Function to compute eroi and net for any given energy
       dataframe already having erde and gross energy.

    Parameters
    ----------   
    energy_res : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        for a given type of energy resource.

    Returns
    -------
    energy_res : DataFrame (Python pandas object)
        Updated energy dataframe, eroi and net energy have been
        derived from erde and gross variables.

    """

    # > < : ces signes ne marchent pas sur mon clavier...

    # Computing eroi for this finite energy source
    energy_res['eroi'] = energy_res['gross']/energy_res['erde']

    # Computing net for this finite energy source
    energy_res['net'] = energy_res['gross'] - energy_res['erde']
#    energy_res['net'] = energy_res['erde']*( energy_res['eroi']-1 )

    return energy_res


def compute_fin_res(fin_res, fin_res_dens, fin_res_avail, fin_res_peak_time, fin_res_min_erde, fin_res_abun):
    """Function to compute the time-dependent variables corresponding
       to the finite energy resource following the model hypotheses.

    Parameters
    ----------   
    fin_res : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        (initialized) for a finite energy resource.

    fin_res_dens : float
        Intrinsic energy density.

    fin_res_avail : float
        Resource availability.

    fin_res_peak_time : float
        Arbitrary moment in time where the gross energy curve peaks.

    fin_res_min_erde : float
        Minimum energy required to start delivering energy.

    fin_res_abun : float
        Worldwide abundance for this finite energy resource.

    Returns
    -------
    fin_res : DataFrame (Python pandas object)
        Updated structured dataset of different energy variables time
        series (erde, gross, eroi and net) for a finite energy
        resource.

    """

    # Computing erde for this finite energy source
    fin_res['erde'] = fin_res_min_erde*np.exp(fin_res['time']/fin_res_abun)
    
    # Computing gross for this finite energy source
    fin_res['gross'] = fin_res_dens*( 1/(fin_res_avail*np.sqrt(2*np.pi)) )*np.exp( -((fin_res['time']-fin_res_peak_time)**2)/(2*fin_res_avail**2) )

    # Computing eroi and net for this finite energy source
    fin_res = compute_eroi_and_net(fin_res)

    # Setting to zero unrealistic values (eroi and net << 0)
    for i in np.arange(len(fin_res['time'])):
        if ( ((fin_res['eroi'][i+1]-fin_res['eroi'][i]) < 0) and (fin_res['net'][i] < 0) ):
            break
    # Except for the time variable
    for var in fin_res.columns:
        if var != 'time':
            fin_res.loc[i:, var] = 0

    return fin_res


def compute_inf_res(inf_res, inf_res_max_ener, inf_res_time_infra_deploy, inf_res_min_erde_sys, inf_res_time_best_tech, inf_res_infra_deploy, inf_res_inter_coop, inf_res_time_infra_tech):
    """Function to compute the time-dependent variables corresponding
       to the infinite energy resource following the model hypotheses.

    Parameters
    ----------   
    inf_res : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        (initialized) for an infinite energy resource.

    inf_res_max_ener : float
        Maximum energy per time unit willing to be reached.

    inf_res_time_infra_deploy : float
        Time required to fully deploy the worldwide infrastructure.

    inf_res_min_erde_sys : float
        Minimum system energy required.

    inf_res_time_best_tech : float
        Minimum time required to achieve best technological
        performance.

    inf_res_infra_deploy : float
        Worldwide infrastructure deployment required.

    inf_res_inter_coop : float
        International cooperation level.

    inf_res_time_infra_tech : float
        Minimum time required to achieve best infrastructure
        technological performance.

    Returns
    -------
    inf_res : DataFrame (Python pandas object)
        Updated structured dataset of different energy variables time
        series (eroi and erde have been added) for an infinite
        energy resource.

    """

    # Computing erde for this infinite energy source
    inf_res['erde'] = inf_res_min_erde_sys*(1-np.exp(-inf_res['time']/inf_res_time_best_tech)) + inf_res_infra_deploy*np.sin(2*np.pi*inf_res['time']/inf_res_inter_coop)*np.exp(-inf_res['time']/inf_res_time_infra_tech)

    # Computing gross production for this infinite energy source
    inf_res['gross'] = inf_res_max_ener*( 1-np.exp(-inf_res['time']/inf_res_time_infra_deploy) )

    # Computing eroi and net for this finite energy source
    inf_res = compute_eroi_and_net(inf_res)
    
    return inf_res


def compute_inf_res_lag(inf_res_lag, inf_res, inf_res_time_lag):
    """Function to compute the time-dependent variables corresponding
       to the infinite energy resource following the model hypotheses.

    Parameters
    ----------   
    inf_res_lag : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        (initialized) for a time-lagged infinite energy resource.

    inf_res : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        for an infinite energy resource.

    inf_res_time_lag : int
        Time lag with respect to the beginning of the simulation.

    Returns
    -------
    inf_res_lag : DataFrame (Python pandas object)
        Updated structured dataset of different energy variables time
        series (eroi and erde have been added) for a time_lagged
        infinite energy resource.

    """

    # Computing erde for this time-lagged infinite energy source
    a = np.zeros(len(inf_res['time']))
    a[inf_res_time_lag:] = inf_res['erde'][:-inf_res_time_lag]
    inf_res_lag['erde'] = a
    
    # Computing gross for this time-lagged infinite energy source
    b = np.zeros(len(inf_res['time']))
    b[inf_res_time_lag:] = inf_res['gross'][:-inf_res_time_lag]
    inf_res_lag['gross'] = b
    
    # Computing eroi and net for this finite energy source
    inf_res_lag = compute_eroi_and_net(inf_res_lag)
    
    return inf_res_lag


def compute_ghg_atm(ghg_atm, fin_res, ghg_max_unit_emis, ghg_emis_intens, ghg_disint_cst, ghg_preind_conc):
    """Function to compute the time-dependent variables corresponding
       to the GHG variables following the model hypotheses.

    Parameters
    ----------   
    ghg_atm : DataFrame (Python pandas object)
        Structured dataset of the GHG variables (initialized).

    fin_res : DataFrame (Python pandas object)
        Structured dataset of the Finite energy resource variables.

    ghg_max_unit_emis : float
        Maximum GHG emission per finite energy unit consumed.

    ghg_emis_intens : float
        Intrinsic emission intensity.

    ghg_desint_cst : float
        GHG disintegration constant in the atmosphere.

    ghg_preind_conc : float
        Preindustrial GHG concentration in the atmosphere.

    Returns
    -------
    ghg_atm : DataFrame (Python pandas object)
        Updated structured dataset of the GHG variables (emis and
        conc have been added).

    """

    # Computing emissions coming from the finite energy source
    ghg_atm['emis'] = ghg_max_unit_emis*np.exp(-ghg_atm['time']/ghg_emis_intens)*fin_res['gross']

    # Computing the GHG atmospheric concentration derived from
    # the emissions
    # Initializing first time step
    ghg_atm.loc[0, 'conc'] = ghg_preind_conc
    # Recursive computation for the full time series
    for i in ghg_atm['time'][1:]:
        ghg_atm.loc[i, 'conc'] = ghg_atm.loc[i-1, 'conc']*np.exp(-(ghg_atm.loc[i-1, 'conc']-ghg_preind_conc)/ghg_disint_cst) + ghg_atm.loc[i, 'emis']
    
    return ghg_atm


def compute_inf_res_var(inf_res_var, inf_res, var_erde_loss_rate):
    """Function to compute the time-dependent variables corresponding
       to the infinite energy resource following the model hypotheses.

    Parameters
    ----------   
    inf_res_var : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        (initialized) for an infinite energy resource influenced
        by the natural atmospheric variability.

    inf_res : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        (initialized) for an infinite energy resource.

    var_erde_loss_rate : float
        Percentage of the erde increase due to the atmospheric
        variability (unpredictability) in order to keep delivering
        the same amount of gross energy .

    Returns
    -------
    inf_res_var : DataFrame (Python pandas object)
        Updated structured dataset of different energy variables time
        series (erde, gross, eroi and net have been added) for an
        infinite energy resource influenced by the atmospheric
        variability.

    """

    # Computing erde for this infinite energy source
    inf_res_var['erde'] = inf_res['erde']*( 1+(var_erde_loss_rate/100) )

    # Computing gross production for this infinite energy source
    inf_res_var['gross'] = inf_res['gross']

    # Computing eroi and net for this natural variability scenario
    inf_res_var = compute_eroi_and_net(inf_res_var)
    
    return inf_res_var


def compute_inf_res_cc(inf_res_cc, inf_res_var, ghg_atm, ghg_impact_factor):
    """Function to compute the time-dependent variables corresponding
       to the infinite energy resource following the model hypotheses.

    Parameters
    ----------   
    inf_res_cc : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        (initialized) for an infinite energy resource influenced by
        both the natural atmospheric variability and climate change.

    inf_res_var : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        (initialized) for an infinite energy resource influenced
        by the natural atmospheric variability.

    ghg_atm : DataFrame (Python pandas object)
        Structured dataset of the GHG variables.

    ghg_impact_factor : float
        Constant estimating the impact on the environment caused by
        the GHG concentration increase in the atmosphere.

    Returns
    -------
    inf_res_cc : DataFrame (Python pandas object)
        Updated structured dataset of different energy variables time
        series (erde, gross, eroi and net have been added) for an
        infinite energy resource influenced by the atmospheric
        variability and climate change.

    """

    # Computing erde for this infinite energy source
    inf_res_cc['erde'] = inf_res_var['erde']*( 1+(ghg_impact_factor*(ghg_atm['conc']-ghg_atm['conc'][0])/ghg_atm['conc'][0]) )

    # Computing gross production for this infinite energy source
    inf_res_cc['gross'] = inf_res_var['gross']

    # Computing eroi and net for this climate change scenario
    inf_res_cc = compute_eroi_and_net(inf_res_cc)
    
    return inf_res_cc


def compute_coupled_res(coupled_res, fin_res, inf_res_lag):
    """Function to compute the time-dependent variables corresponding
       to the coupled (finite + infinite) energy resource following
       the model hypotheses.

    Parameters
    ----------   
    coupled_res : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        combining a finite and an infinite energy resource where
        the latter starts to be developed with a time lag.

    fin_res : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        for an finite energy resource.

    inf_res_lag : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        for a time-lagged infinite energy resource.

    Returns
    -------
    coupled_res : DataFrame (Python pandas object)
        Updated coupled_res variable (erde, gross, eroi and net
        have been added).

    """

    # Computing erde for this coupled energy source
    coupled_res['erde'] = fin_res['erde'] + inf_res_lag['erde']

    # Computing gross production for this infinite energy source
    coupled_res['gross'] = fin_res['gross'] + inf_res_lag['gross']

    # Computing eroi and net for this climate change scenario
    coupled_res = compute_eroi_and_net(coupled_res)
    
    return coupled_res


def plot_energy_evol(fig_num, energy_type, energy_res, current_time):
    """Function to plot the energy variables time evolution.

    Parameters
    ----------   

    fig_num : str
        Suptitle's first word, figure number.

    energy_type : str
        Energy type name, finite or infinite.

    energy_res : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        (eroi, erde, net) for a given type of energy resource.

    current_time : int
        Specific time (by default associated to our current time)
        where to plot a vertical line on the figures.

    """

    fig = plt.figure(figsize=(5, 8))
    plt.suptitle(fig_num+': '+energy_type+' energy resource time evolution', fontsize=12)
    
    # Top subplot 
    plt.subplot(4, 1, 1)
    plt.plot(energy_res['time'], energy_res['erde'], color='blue')
    plt.plot(np.ones(int(np.max(energy_res['erde'])+1))*current_time, np.arange(int(np.max(energy_res['erde'])+1)), color='black',lw=0.5)
    plt.title('a) Energy Required to Deliver Energy', fontsize=10)
#    plt.text(nb_cube_mod[-1]*0.5, alt[0]+1000, 'Total number of\ncasings within the\ntower : '+str(nb_cube_tower), fontsize=10)
    plt.ylabel('ERDE [energy unit]')
    plt.xlabel('Time')

    # Middle upper subplot 
    plt.subplot(4, 1, 2)
    plt.plot(energy_res['time'], energy_res['gross'], color='black')
    plt.plot(np.ones(int(np.max(energy_res['gross'])+1))*current_time, np.arange(int(np.max(energy_res['gross'])+1)), color='black',lw=0.5)
    plt.title('b) Gross Energy Production', fontsize=10)
    plt.ylabel('GEP [energy unit]')
    plt.xlabel('Time')

    # Middle lower subplot 
    plt.subplot(4, 1, 3)
    plt.plot(energy_res['time'], energy_res['eroi'], color='green')
    plt.plot(energy_res['time'], np.ones(len(energy_res['time'])), color='black', linewidth=0.5)
    plt.plot(np.ones(int(np.max(energy_res['eroi'])+1))*current_time, np.arange(int(np.max(energy_res['eroi'])+1)), color='black',lw=0.5)
    plt.title('c) Energy Return On energy Investment', fontsize=10)
    plt.ylabel('EROI [no unit]')
    plt.xlabel('Time')
    
    # Bottom subplot
    plt.subplot(4, 1, 4)
    plt.plot(energy_res['time'], energy_res['net'], color='red')
    plt.plot(energy_res['time'], np.zeros(len(energy_res['time'])), color='black', linewidth=0.5)
    plt.plot(np.ones(int(np.max(energy_res['net'])+1))*current_time, np.arange(int(np.max(energy_res['net'])+1)), color='black',lw=0.5)
    plt.title('d) Net Available Energy', fontsize=10)
    plt.ylabel('NAE [energy unit]')
    plt.xlabel('Time')

#    plt.legend(fontsize=9)
    plt.tight_layout()
    
    # Save figure
    if energy_type == 'Finite':
        plt.savefig('Fig1_'+energy_type+'_EROI_simple_model.png')
    elif energy_type == 'Infinite':
        plt.savefig('Fig3_'+energy_type+'_EROI_simple_model.png')
    elif energy_type == 'Coupled':
        plt.savefig('Fig5_'+energy_type+'_EROI_simple_model.png')
    

def plot_ghg_evol(fig_num, ghg_atm, current_time):
    """Function to plot the atmospheric and LTB profiles.

    Parameters
    ----------   

    fig_num : str
        Suptitle's first word, figure number.

    ghg_atm : DataFrame (Python pandas object)
        Structured dataset of the different GHG variables time series
        (emis, conc).

    current_time : int
        Specific time (by default associated to our current time)
        where to plot a vertical line on the figures.

    """

    fig = plt.figure(figsize=(5, 6))
    plt.suptitle(fig_num+': Time evolution of the GHG variables', fontsize=12)
    
    # Top subplot 
    plt.subplot(2, 1, 1)
    plt.plot(ghg_atm['time'], ghg_atm['emis'], color='black')
    plt.plot(np.ones(int(np.max(ghg_atm['emis'])+1))*current_time, np.arange(int(np.max(ghg_atm['emis'])+1)), color='black',lw=0.5)
    plt.title('a) GreenHouse Gas Emissions in the atmosphere', fontsize=10)
    plt.ylabel('Emissions [mass unit]')
    plt.xlabel('Time')
    
    # Bottom subplot
    plt.subplot(2, 1, 2)
    plt.plot(ghg_atm['time'], ghg_atm['conc'], color='red')
    plt.plot(np.ones(int(np.max(ghg_atm['conc'])-ghg_atm['conc'][0]))*current_time, np.arange(int(np.max(ghg_atm['conc']-ghg_atm['conc'][0])))+ghg_atm['conc'][0], color='black',lw=0.5)
    plt.title('b) GreenHouse Gas Concentration in the atmosphere', fontsize=10)
    plt.ylabel('Concentration [ppmv]')
    plt.xlabel('Time')

#    plt.legend(fontsize=9)
    plt.tight_layout()
    
    # Save figure
    plt.savefig('Fig2_GHG_EROI_simple_model.png')


def plot_inf_res_var(fig_num, energy_type, energy_res, energy_res_var, current_time):
    """Function to plot the energy variables time evolution.

    Parameters
    ----------   

    fig_num : str
        Suptitle's first word, figure number.

    energy_type : str
        Energy type name, finite or infinite.

    energy_res : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        (eroi, erde, net) for a given type of energy resource.

    energy_res_var : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        (eroi, erde, net) for the energy resource + variability.

    current_time : int
        Specific time (by default associated to our current time)
        where to plot a vertical line on the figures.

    """

    fig = plt.figure(figsize=(5, 8))
    plt.suptitle(fig_num+': '+energy_type+' energy resource time evolution', fontsize=12)
    
    # Top subplot 
    plt.subplot(4, 1, 1)
    plt.plot(energy_res['time'], energy_res['erde'], color='blue', label='Ideal')
    plt.plot(energy_res_var['time'], energy_res_var['erde'], '--', color='blue', label='Natural variability')    
    plt.plot(np.ones(int(np.max(energy_res['erde'])+1))*current_time, np.arange(int(np.max(energy_res['erde'])+1)), color='black',lw=0.5)
    plt.title('a) Energy Required to Deliver Energy', fontsize=10)
#    plt.text(nb_cube_mod[-1]*0.5, alt[0]+1000, 'Total number of\ncasings within the\ntower : '+str(nb_cube_tower), fontsize=10)
    plt.ylabel('ERDE [energy unit]')
    plt.xlabel('Time')
    plt.legend(loc='upper right', fontsize=7)

    # Middle upper subplot 
    plt.subplot(4, 1, 2)
    plt.plot(energy_res['time'], energy_res['gross'], color='black', label='Ideal')
    plt.plot(energy_res_var['time'], energy_res_var['gross'], '--', color='black', label='Natural variability')    
    plt.plot(np.ones(int(np.max(energy_res['gross'])+1))*current_time, np.arange(int(np.max(energy_res['gross'])+1)), color='black',lw=0.5)
    plt.title('b) Gross Energy Production', fontsize=10)
    plt.ylabel('GEP [energy unit]')
    plt.xlabel('Time')

    # Middle lower subplot 
    plt.subplot(4, 1, 3)
    plt.plot(energy_res['time'], energy_res['eroi'], color='green', label='Ideal')
    plt.plot(energy_res_var['time'], energy_res_var['eroi'], '--', color='green', label='Natural variability')    
    plt.plot(energy_res['time'], np.ones(len(energy_res['time'])), color='black', linewidth=0.5)
    plt.plot(np.ones(int(np.max(energy_res['eroi'])+1))*current_time, np.arange(int(np.max(energy_res['eroi'])+1)), color='black',lw=0.5)
    plt.title('c) Energy Return On energy Investment', fontsize=10)
    plt.ylabel('EROI [no unit]')
    plt.xlabel('Time')
    plt.legend(loc='lower right', fontsize=7)
    
    # Bottom subplot
    plt.subplot(4, 1, 4)
    plt.plot(energy_res['time'], energy_res['net'], color='red', label='Ideal')
    plt.plot(energy_res_var['time'], energy_res_var['net'], '--', color='red', label='Natural variability')    
    plt.plot(energy_res['time'], np.zeros(len(energy_res['time'])), color='black', linewidth=0.5)
    plt.plot(np.ones(int(np.max(energy_res['net'])+1))*current_time, np.arange(int(np.max(energy_res['net'])+1)), color='black',lw=0.5)
    plt.title('d) Net Available Energy', fontsize=10)
    plt.ylabel('NAE [energy unit]')
    plt.xlabel('Time')
    plt.legend(loc='lower right', fontsize=7)
    
    plt.tight_layout()
    
    # Save figure
    plt.savefig('Fig4_'+energy_type+'_var_EROI_simple_model.png')

    
    
def plot_case_study(fig_num, ghg_atm, inf_res, inf_res_var, inf_res_cc):
    """Function to plot the different Idealized scenarios.

    Parameters
    ----------   

    fig_num : str
        Suptitle's first word, figure number.

    ghg_atm : DataFrame (Python pandas object)
        Structured dataset of the different GHG variables time series
        (emis, conc).

    inf_res : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        for an infinite energy resource.

    inf_res_var : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        for an infinite energy resource influenced by the natural 
        atmospheric variability.

    inf_res_cc : DataFrame (Python pandas object)
        Structured dataset of different energy variables time series
        for an infinite energy resource influenced by both the 
        natural atmospheric variability and climate change.

    """

    fig = plt.figure(figsize=(6, 8))
    plt.suptitle(fig_num+': Time evolution of the GHG concentration\nand its influence on ERDE, EROI and Net energy production', fontsize=12)
        
    # Top subplot
    plt.subplot(4, 1, 1)
    plt.plot(ghg_atm['time'], ghg_atm['conc'], color='red')
    plt.title('a) GreenHouse Gas Concentration in the atmosphere', fontsize=10)
    plt.ylabel('Concentration [ppmv]')
    plt.xlabel('Time')

    # Middle upper subplot 
    plt.subplot(4, 1, 2)
    plt.plot(inf_res['time'], inf_res['erde'], color='blue', label='Ideal')
    plt.plot(inf_res_var['time'], inf_res_var['erde'], '--', color='blue', label='Natural variability')
    plt.plot(inf_res_cc['time'], inf_res_cc['erde'], '-.', color='blue', label='Climate Change')
    plt.title('b) Energy Required to Deliver Energy', fontsize=10)
    plt.ylabel('ERDE [energy unit]')
    plt.xlabel('Time')
    plt.legend(loc='best', fontsize=7)

    # Middle lower subplot 
    plt.subplot(4, 1, 3)
    plt.plot(inf_res['time'], inf_res['eroi'], color='green', label='Ideal')
    plt.plot(inf_res_var['time'], inf_res_var['eroi'], '--', color='green', label='Natural variability')
    plt.plot(inf_res_cc['time'], inf_res_cc['eroi'], '-.', color='green', label='Climate Change')
    plt.plot(inf_res['time'], np.ones(len(inf_res['time'])), color='black', linewidth=0.5)
    plt.title('c) Energy Return On energy Investment', fontsize=10)
    plt.ylabel('EROI [no unit]')
    plt.xlabel('Time')
    plt.legend(loc='best', fontsize=7)
    
    # Bottom subplot
    plt.subplot(4, 1, 4)
    plt.plot(inf_res['time'], inf_res['net'], color='red', label='Ideal')
    plt.plot(inf_res_var['time'], inf_res_var['net'], '--', color='red', label='Natural variability')
    plt.plot(inf_res_cc['time'], inf_res_cc['net'], '-.', color='red', label='Climate Change')
    plt.plot(inf_res['time'], np.zeros(len(inf_res['time'])), color='black', linewidth=0.5)
    plt.title('d) Net Available Energy', fontsize=10)
    plt.ylabel('NAE [energy unit]')
    plt.xlabel('Time')    
    plt.legend(loc='lower right', fontsize=7)
    
    plt.tight_layout()
    
    # Save figure
    plt.savefig('Fig6_var+CC_EROI_simple_model.png')
