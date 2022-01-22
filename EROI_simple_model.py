#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Script to compute the main features of the indealized EROIs

This script needs a few parameters to be set in the variable.py
file in order to compute the time evolution of the main variables
(EROI, ERDE, GHG) for the different types of energy source and
scenarios considered.
Plots to illustrate the output are also created at the end of it.

November 2021, v1 by Rodrigo Guzman

"""

import numpy as np
import pandas as pd

#import constants as cst
import variables as var
import functions as func

if __name__ == "__main__":
    """Script to compute main features of the EROI simple model.

    """

    ###############################################################
    ####### I - Defining and initializing all the variables #######
    ###############################################################

    # Creating all the dataframe variables
    # Finite energy resource
    fin_res = pd.DataFrame()

    # Infinite energy resource (ideal)
    inf_res = pd.DataFrame()

    # Infinite energy resource (lagged in time)
    inf_res_lag = pd.DataFrame()

    # GHG variables influenced by finite energy resource
    ghg_atm = pd.DataFrame()

    # Infinite energy resource (natural variability)
    inf_res_var = pd.DataFrame()
    
    # Infinite energy resource (natural variability +
    # GHG concentration increase-induced climate change influence)
    inf_res_cc = pd.DataFrame()

    # Coupled (Finite + lagged Infinite) energy resource
    coupled_res = pd.DataFrame()

    # Coupled (Finite + lagged Infinite) energy resource
    # (natural variability)
    coupled_res_var = pd.DataFrame()

    # Coupled (Finite + lagged Infinite) energy resource
    # (natural variability + GHG concentration influence)
    coupled_res_cc = pd.DataFrame()

    # Creating a list with all the variables
    var_list = [fin_res, inf_res, inf_res_lag, ghg_atm, inf_res_var, inf_res_cc, coupled_res, coupled_res_var, coupled_res_cc]

    # Initializing all the variables with a time vector
    var_list = func.set_df_time(var_list, var.time_length)
    
    ###############################################################
    #### II - Non-renewable variables derived from variables.py ###
    ###############################################################

    # Computing all the energy variables (erde, gross, eroi and net)
    fin_res = func.compute_fin_res(fin_res, var.fin_res_dens, var.fin_res_avail, var.fin_res_peak_time, var.fin_res_min_erde, var.fin_res_abun)
    
    ###############################################################
    ##### III - Renewable variables derived from variables.py #####
    ###############################################################

    # Computing all the energy variables (erde, gross, eroi and net)
    inf_res = func.compute_inf_res(inf_res, var.inf_res_max_ener, var.inf_res_time_infra_deploy, var.inf_res_min_erde_sys, var.inf_res_time_best_tech, var.inf_res_infra_deploy, var.inf_res_inter_coop, var.inf_res_time_infra_tech)

    # Computing time-lagged energy variables
    inf_res_lag = func.compute_inf_res_lag(inf_res_lag, inf_res, var.inf_res_time_lag)
    
    ###############################################################
    ######### IV - GHG variables derived from variables.py ########
    ###############################################################

    # Compute ghg variables (emis and conc)
    ghg_atm = func.compute_ghg_atm(ghg_atm, fin_res, var.ghg_max_unit_emis, var.ghg_emis_intens, var.ghg_disint_cst, var.ghg_preind_conc)
    
    ##############################################################
    #####  V - Renewable energy scenario with uncertainties  #####
    ##############################################################

    # Compute energy variables with natural variability impact on
    # renewable energy ERDE
    inf_res_var = func.compute_inf_res_var(inf_res_var, inf_res_lag, var.var_erde_loss_rate)
    
    # Compute energy variables with natural variability and climate
    # change impacts on renewable energy ERDE
    inf_res_cc = func.compute_inf_res_cc(inf_res_cc, inf_res_var, ghg_atm, var.ghg_impact_factor)

    ##############################################################
    #####  VI - Non-renewable and renewable coupled scenario  ####
    ##############################################################

    coupled_res = func.compute_coupled_res(coupled_res, fin_res, inf_res_lag)
    
    # Compute all the energy variables (erde, gross, eroi and net)
    coupled_res_var = func.compute_coupled_res(coupled_res_var, fin_res, inf_res_var)
    
    # Compute all the energy variables (erde, gross, eroi and net)
    coupled_res_cc = func.compute_coupled_res(coupled_res_cc, fin_res, inf_res_cc)

    ######################################
    ###########  VI - Plots  #############
    ######################################
    
    # Ploting energy variables time evolution for a finite energy
    func.plot_energy_evol('Figure 2', 'non durable', fin_res, var.current_time)

    # Ploting GHG time evolution for the GHG-emitting
    # finite energy resource
    func.plot_ghg_evol('Figure 3', ghg_atm, var.current_time)

    # Ploting energy variables time evolution for
    # an infinite energy resource
    func.plot_energy_evol('Figure 4', 'durable', inf_res_lag, var.current_time)

    # Ploting energy variables time evolution for
    # an infinite energy resource + natural variability
    func.plot_inf_res_var('Figure 5', 'durable', inf_res_lag, inf_res_var, var.current_time)

    # Ploting energy variables time evolution for an energy system
    # in transition from finite to infinite energy resources
    func.plot_energy_evol('Figure 6', 'en transition', coupled_res, var.current_time)
    
    # Ploting infinite energy eroi with uncertainties (case 1 or 2)
    if var.inf_res_time_lag==100:
        func.plot_case_study('Figure 7', 'case 1', ghg_atm, coupled_res, coupled_res_var, coupled_res_cc, var.current_time)
    else:
        func.plot_case_study('Figure 8', 'case 2', ghg_atm, coupled_res, coupled_res_var, coupled_res_cc, var.current_time)
