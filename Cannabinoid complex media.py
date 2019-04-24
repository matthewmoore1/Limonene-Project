#Setting complex media

#import cobra models and required tools
import cobra
import os
from os.path import join
import matplotlib.pyplot as plt
import pandas as pd
from __future__ import print_function
from cobra import Model, Reaction, Metabolite
import numpy as np

#import iJO1366 - can model
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_can_v1.xml"))

#Setting complex media
model.reactions.EX_ca2_e.lower_bound = -0.057
model.reactions.EX_mg2_e.lower_bound = -0.032
model.reactions.EX_k_e.lower_bound = -2.92
model.reactions.EX_na1_e.lower_bound = -10.96
model.reactions.EX_ala__L_e.lower_bound = -4.85
model.reactions.EX_arg__L_e.lower_bound = -2.47
model.reactions.EX_asp__L_e.lower_bound = -6.17
model.reactions.EX_cys__L_e.lower_bound = -0.43
model.reactions.EX_glu__L_e.lower_bound = -16.33
model.reactions.EX_gly_e.lower_bound = -3.23
model.reactions.EX_his__L_e.lower_bound = -2.36
model.reactions.EX_ile__L_e.lower_bound = -4.8
model.reactions.EX_leu__L_e.lower_bound = -7.7
model.reactions.EX_lys__L_e.lower_bound = -5.6
model.reactions.EX_met__L_e.lower_bound = -2.01
model.reactions.EX_phe__L_e.lower_bound = -3.28
model.reactions.EX_pro__L_e.lower_bound = -9.23
model.reactions.EX_ser__L_e.lower_bound = -6.23
model.reactions.EX_thr__L_e.lower_bound = -4.27
model.reactions.EX_trp__L_e.lower_bound = -0.74
model.reactions.EX_tyr__L_e.lower_bound = -1.30
model.reactions.EX_val__L_e.lower_bound = -6.22
model.reactions.EX_glyc_e.lower_bound = -43.4
model.reactions.EX_glc_e.lower_bound = -16.65
model.reactions.EX_o2_e.lower_bound = -1000  #set as theoretical maximum so not rate limiting
model.reactions.EX_pi_e.lower_bound = -1000  #set as theoretical maximum so not rate limiting
model.reactions.EX_thm_e.lower_bound = -0.19
model.reactions.EX_pydxn_e.lower_bound = -0.34
model.reactions.EX_pnto__R_e.lower_bound = -1.15
model.reactions.EX_chol_e.lower_bound = -34.56
model.reactions.EX_btn_e.lower_bound = -0.04

model.reactions.EX_cbl1_e.lower_bound = -1000 #metal ions set as theoretical maximum so not rate limiting
model.reactions.EX_cl_e.lower_bound = -1000
model.reactions.EX_cobalt2_e.lower_bound = -1000
model.reactions.EX_cu2_e.lower_bound = -1000
model.reactions.EX_fe2_e.lower_bound = -1000
model.reactions.EX_fe3_e.lower_bound = -1000
model.reactions.EX_mn2_e.lower_bound = -1000
model.reactions.EX_mobd_e.lower_bound = -1000
model.reactions.EX_ni2_e.lower_bound = -1000
model.reactions.EX_sel_e.lower_bound = -1000
model.reactions.EX_slnt_e.lower_bound = -1000
model.reactions.EX_tungs_e.lower_bound = -1000
model.reactions.EX_zn2_e.lower_bound = -1000

model.reactions.EX_nh4_e.lower_bound = -1.080082 #flux in minimal media - amino acids added so not rate limiting
model.reactions.EX_so4_e.lower_bound = -0.025221 #flux in minimal media - amino acids added so not rate limiting

#model biomass
model.objective = "Ec_biomass_iJO1366_core_53p95M"
model.optimize().f
model.summary()

#model OAC flux
model.objective = "OAC"
model.optimize().f
model.summary()

#model CBGAS flux
model.objective = "CBGAS"
model.optimize().f
model.summary()

#save model in complex media
cobra.io.write_sbml_model(model,filename="C:/Work/ecoli_mod/eco_can_v2.xml")

#Check maximal production for olivetolic acid - Define the maximal production of the OTA for the constraints
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_can_v2.xml"))

model.objective = "OAC"
OACo = model.optimize().f  

#Application
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_can_v2.xml"))
model.objective = "Ec_biomass_iJO1366_core_53p95M"
 
OAC = list(np.linspace(0, OACo, 100))
bm = {}
flux = pd.DataFrame()

b= 0
for a in OAC:
    model.reactions.OAC.upper_bound = a
    model.reactions.OAC.lower_bound = a
    flux[b] =  model.optimize().fluxes
    bm[b] = model.optimize().f
    b=b+1

#BiomassExport = pd.DataFrame(bm)
oac = pd.DataFrame(OAC)
bm=pd.Series(bm)
oac['Biomass']=bm
bm.to_csv("C:/Work/ecoli_mod/OAC_biomass_complex.csv")
oac.to_csv("C:/Work/ecoli_mod/OAC_biomass_complex.csv") #of Lim constraints in loop & Optimization of biomass

#Check maximal production for olivetolic acid - Define the maximal production of the OTA for the constraints
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_can_v2.xml"))

model.objective = "CBGAS"
CBGASo = model.optimize().f 

#Application
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_can_v2.xml"))
model.objective = "Ec_biomass_iJO1366_core_53p95M"
 
CBGAS = list(np.linspace(0, CBGASo, 100))
bm = {}
flux = pd.DataFrame()

b= 0
for a in CBGAS:
    model.reactions.CBGAS.upper_bound = a
    model.reactions.CBGAS.lower_bound = a
    flux[b] =  model.optimize().fluxes
    bm[b] = model.optimize().f
    b=b+1

#BiomassExport = pd.DataFrame(bm)
cbgas = pd.DataFrame(CBGAS)
bm=pd.Series(bm)
cbgas['Biomass']=bm
bm.to_csv("C:/Work/ecoli_mod/CBGAS_biomass_complex.csv")
cbgas.to_csv("C:/Work/ecoli_mod/CBGAS_biomass_complex.csv") #of Lim constraints in loop & Optimization of biomass

