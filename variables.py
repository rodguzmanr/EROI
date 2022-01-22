#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Variables needed by the EROI_simple_model.py script to compute
   and plot the figures of the studied indealized scenarios.

This file contains all the parameters of the EROI simple model
needed to compute the time evolution of all the energy system
variables described by this model for non-renewable and renewable
energy sources.

November 2021, v1 by Rodrigo Guzman

"""

import numpy as np

### Parameters to be set by the user in order to run the   ###
### EROI simple model, which performs calculations and     ###
### plots the main variables for a (I) non-renewable       ###
### energy source emitting (III) Greehouse gases and a     ###
### (II) non-emitting renewable energy source.             ###


# Time extent for the full simulation
time_length = 500

# Current time
current_time = 100

##############################################################
##########   I - Non-renewable energy hypotheses   ###########
##############################################################

#################   a) GEP hypotheses   #####################

# For this type of energy resource, we consider the Gross
# Energy Production (GEP) time evolution from the begining of
# exploiting the resource till its full depletion as a gaussian.

# Resource availability, arbitrary quantity measuring how
# easily available the resource is throughout the world
fin_res_avail = 50

# Intrinsic energy density, arbitrary quantity measuring how
# much energy can be found within one unit of the energy source
fin_res_dens = 3000

# EROI peak, arbitrary moment in time where the EROI curve peaks
fin_res_peak_time = 150

#################   b) ERDE hypotheses   #####################

# For this type of energy resource, we consider the Energy
# Required to Deliver Energy (ERDE) time evolution from the
# begining of exploiting this resource follows an exponential.

# Minimum energy required to start delivering the resource energy
fin_res_min_erde = 1

# Energy abundance, arbitrary quantity measuring how abundant
# the resource is throughout the world
fin_res_abun = 100


##############################################################
###########   II - Renewable energy hypotheses   #############
##############################################################

# Time lag with respect to the beginning of the simulation
inf_res_time_lag = 100 # case 1=100, case 2=150, 200

#################   a) GEP hypotheses   #####################

# For this type of energy resource, we consider the GEP time
# evolution is an asymptotic ascending curve.

# Maximum energy per time unit target, arbitrary quantity
# measuring how much energy should be enough to achieve a
# worldwide sustainable system
inf_res_max_ener = 30

# Time required to fully deploy the worldwide infrastructure
inf_res_time_infra_deploy = 80

#################   b) ERDE hypotheses   #####################

# For this type of energy resource, we consider the ERDE time
# evolution from the begining of exploiting this resource
# follows an inversed sink potential curve.

# Minimum system energy required to keep delivering the resource
# energy once an hypothetical equilibrium state is reached
inf_res_min_erde_sys = 5

# Minimum time required to achieve best technological
# performance, arbitrary quantity measuring the transforming
# performance of a resource energy unit
inf_res_time_best_tech = 50

# International cooperation level, arbitrary quantity measuring
# how easy international relationships are to develop the
# required infrastructure to be deployed
inf_res_inter_coop = 500

# Worldwide deployment required, arbitrary quantity measuring
# how widely the technological deployment is needed before 
# the resource can be exploited at its full potential
inf_res_infra_deploy = 100

# Minimum time required to achieve best technological
# performance for the infrastructure to be deployed
inf_res_time_infra_tech = 50

##############################################################
#########   III - Greenhouse gas (GHG) hypotheses   ##########
##############################################################

#############   a) GHG emissions hypotheses   ################

# We consider the GHG emitted by the non-renewable energy
# source is proportional to the amount of this resource
# consumed, weighted by an emitting factor which exponetially
# decreases with time.

# Maximum GHG emission per finite energy unit consumed
ghg_max_unit_emis = 1

# Intrinsic emission intensity, arbitrary quantity measuring how
# difficult is to make the GHG emissions decrease with time
# due to efficiency, technological development, etc
ghg_emis_intens = 100

############   b) GHG concentration hypotheses   #############

# We consider the GHG concentration in the atmosphere to be a
# recursive function depending on the passed emissions and
# having an exponential decrease driven by a disintegration
# constant.

# GHG disintegration constant, arbitrary quantity measuring
# how slowly the GHG gets caught up by natural sinks 
ghg_disint_cst = 5e5

# Preindustrial GHG concentration, arbitrary quantity defining
# how much GHG was in the atmosphere before starting to consume
# the non-renewable energy resource
ghg_preind_conc = 280

# Constant estimating the impact on the environment caused by
# the GHG concentration increase in the atmosphere
ghg_impact_factor = 0.1

##############################################################
####   IV - Current renewable energy sources case study   ####
##############################################################

# Adding some factors related to energy production uncertainty
# to describe in a more realistic way the potential EROI time
# evolution.

# Efficiency loss in energy production due to the intrisic
# unpredictability (variability) of the renewable energy source
var_erde_loss_rate = 30.


