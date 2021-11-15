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

    # Creating all the dataframes variables
    # Finite energy resource
    fin_res = pd.DataFrame()

    # Infinite energy resource (ideal)
    inf_res = pd.DataFrame()

    # GHG variables influenced by finite energy resource
    ghg_atm = pd.DataFrame()

    # Infinite energy resource (natural variability)
    inf_res_var = pd.DataFrame()
    
    # Infinite energy resource (natural variability +
    # GHG concentration increase-induced climate change influence)
    inf_res_cc = pd.DataFrame()

    # Coupled (Finite + lagged Infinite) energy resource
    coupled_res = pd.DataFrame()

    # Creating a list with all the variables
    var_list = [fin_res, inf_res, ghg_atm, inf_res_var, inf_res_cc, coupled_res]

    # Initializing all the variables with a time vector
    var_list = func.set_df_time(var_list, var.time_length)
    
    ###############################################################
    #### II - Non-renewable variables derived from variables.py ###
    ###############################################################

    # Compute eroi and erde
    fin_res = func.compute_fin_res(fin_res, var.fin_res_dens, var.fin_res_avail, var.fin_res_peak_time, var.fin_res_min_erde, var.fin_res_abun)

    
    ###############################################################
    ##### III - Renewable variables derived from variables.py #####
    ###############################################################

    # Compute eroi and erde
    inf_res = func.compute_inf_res(inf_res, var.inf_res_max_ener, var.inf_res_avail, var.inf_res_min_erde_sys, var.inf_res_time_best_tech, var.inf_res_infra_deploy, var.inf_res_inter_coop, var.inf_res_time_infra_tech, var.inf_res_time_infra_deploy)

    
    ###############################################################
    ######### IV - GHG variables derived from variables.py ########
    ###############################################################

    # Compute ghg variables
    ghg_atm = func.compute_ghg_atm(ghg_atm, fin_res, var.ghg_max_unit_emis, var.ghg_emis_intens, var.ghg_disint_cst, var.ghg_preind_conc)

    
    ##############################################################
    ######  V - Current renewable energy sources case study  #####
    ##############################################################

    # Compute inf_res_var variables
    inf_res_var = func.compute_inf_res_var(inf_res_var, inf_res, var.var_erde_loss_rate)
    
    # Compute inf_res_cc variables
    inf_res_cc = func.compute_inf_res_cc(inf_res_cc, inf_res_var, ghg_atm)

    
    ######################################
    ###########  VI - Plots  #############
    ######################################
    
    # Ploting energy variables time evolution for a finite energy
    func.plot_energy_evol('Figure 1', 'Finite', fin_res)

    # Ploting energy variables time evolution for an infinite energy
    func.plot_energy_evol('Figure 2', 'Infinite', inf_res)

    # Ploting GHG time evolution for the emitting finite energy
    func.plot_ghg_evol('Figure 3', ghg_atm)

    # Ploting infinite enrergy eroi with uncertainties
    func.plot_case_study('Figure 4', ghg_atm, inf_res, inf_res_var, inf_res_cc)
