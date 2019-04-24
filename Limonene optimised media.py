#Setting optimised media

#import cobra and required tools
import cobra
import os
from os.path import join
import matplotlib.pyplot as plt
import pandas as pd
from __future__ import print_function
from cobra import Model, Reaction, Metabolite
import numpy as np

#import iJO1366 - Limo complex media model
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_Lim_v2.xml"))

#setting improved media constraints
model.reactions.EX_pro__L_e.lower_bound = -27.69 #flux of amino acid tripled from complex media
model.reactions.EX_ser__L_e.lower_bound = -18.69 #flux of amino acid tripled from complex media
model.reactions.EX_glu__L_e.lower_bound = -48.99 #flux of amino acid tripled from complex media
model.reactions.EX_asp__L_e.lower_bound = -18.57 #flux of amino acid tripled from complex media
model.reactions.EX_ala__L_e.lower_bound = -14.55 #flux of amino acid tripled from complex media
model.reactions.EX_thr__L_e.lower_bound = -12.81 #flux of amino acid tripled from complex media
model.reactions.EX_arg__L_e.lower_bound = -7.41 #flux of amino acid tripled from complex media
model.reactions.EX_trp__L_e.lower_bound = -2.22 #flux of amino acid tripled from complex media
model.reactions.EX_cys__L_e.lower_bound = -1.29 #flux of amino acid tripled from complex media
model.reactions.EX_glyc_e.lower_bound = -130.2 #flux of amino acid tripled from complex media
model.reactions.EX_glc_e.lower_bound = -49.95 #flux of amino acid tripled from complex media
model.reactions.EX_gly_e.lower_bound = -9.69 #flux of amino acid tripled from complex media
model.reactions.EX_mg2_e.lower_bound = -0.04 #flux increased from 0.032 to 0.04

#setting biomass as objective function
model.objective = "Ec_biomass_iJO1366_core_53p95M"
model.optimize().f

#setting limonene flux as objective function
model.objective = "LIMS"
model.optimize().f

#save model in optimized media
cobra.io.write_sbml_model(model,filename="C:/Work/ecoli_mod/eco_Lim_v5.xml")

#Check maximal production for limonene - Define the maximal production of the LIM for the constraints
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_Lim_v5.xml"))

model.objective = "LIMS"
LIMo = model.optimize().f 

#Application
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_Lim_v5.xml"))
model.objective = "Ec_biomass_iJO1366_core_53p95M"
 
LIM = list(np.linspace(0, LIMo, 100))
bm = {}
flux = pd.DataFrame()

b= 0
for a in LIM:
    model.reactions.LIMS.upper_bound = a
    model.reactions.LIMS.lower_bound = a
    flux[b] =  model.optimize().fluxes
    bm[b] = model.optimize().f
    b=b+1

#BiomassExport = pd.DataFrame(bm)
Lim = pd.DataFrame(LIM)
bm=pd.Series(bm)
Lim['Biomass']=bm
bm.to_csv("C:/Work/ecoli_mod/LIM_biomass.csv")
Lim.to_csv("C:/Work/ecoli_mod/LIM_biomass_v5.csv") #of Lim constraints in loop & Optimization of biomass
flux.to_csv('C:/Work/ecoli_mod/flux_biomass_v5.csv') #metabolite flux