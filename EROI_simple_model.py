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
    #### I - Non-renewable variables derived from variables.py ####
    ###############################################################

    # Create an empty dataframe
    fin_res = pd.DataFrame()

    # Set the time dimension
    fin_res = func.set_df_time(fin_res, var.time_length)

    # Compute eroi and erde
    fin_res = func.compute_fin_res(fin_res, var.fin_res_dens, var.fin_res_avail, var.fin_res_peak_time, var.fin_res_min_erde, var.fin_res_abun)

    
    ###############################################################
    ##### II - Renewable variables derived from variables.py ######
    ###############################################################

    # Create an empty dataframe
    inf_res = pd.DataFrame()

    # Set the time dimension
    inf_res = func.set_df_time(inf_res, var.time_length)

    # Compute eroi and erde
    inf_res = func.compute_inf_res(inf_res, var.inf_res_max_ener, var.inf_res_avail, var.inf_res_min_erde_sys, var.inf_res_time_best_tech, var.inf_res_infra_deploy, var.inf_res_inter_coop, var.inf_res_time_infra_tech, var.inf_res_time_infra_deploy)

    
    ###############################################################
    ######## III - GHG variables derived from variables.py ########
    ###############################################################

    # Create an empty dataframe
    ghg_atm = pd.DataFrame()

    # Set the time dimension
    ghg_atm = func.set_df_time(ghg_atm, var.time_length)

    # Compute ghg variables
    ghg_atm = func.compute_ghg_atm(ghg_atm, fin_res, var.ghg_max_unit_emis, var.ghg_emis_intens, var.ghg_disint_cst, var.ghg_preind_conc)

    
    ##############################################################
    #####  IV - Current renewable energy sources case study  #####
    ##############################################################

    # Create an empty dataframe
    inf_res_var = pd.DataFrame()
    inf_res_cc = pd.DataFrame()

    # Set the time dimension
    inf_res_var = func.set_df_time(inf_res_var, var.time_length)
    inf_res_cc = func.set_df_time(inf_res_cc, var.time_length)

    # Compute inf_res_var variables
    inf_res_var = func.compute_inf_res_var(inf_res_var, inf_res, var.var_erde_loss_rate)
    
    # Compute inf_res_cc variables
    inf_res_cc = func.compute_inf_res_cc(inf_res_cc, inf_res_var, ghg_atm)

    
    ######################################
    ############  V - Plots  #############
    ######################################
    
    # Ploting energy variables time evolution for a finite energy
    func.plot_energy_evol('Figure 1', 'Finite', fin_res)

    # Ploting energy variables time evolution for an infinite energy
    func.plot_energy_evol('Figure 2', 'Infinite', inf_res)

    # Ploting GHG time evolution for the emitting finite energy
    func.plot_ghg_evol('Figure 3', ghg_atm)

    # Ploting infinite enrergy eroi with uncertainties
    func.plot_case_study('Figure 4', ghg_atm, inf_res, inf_res_var, inf_res_cc)
