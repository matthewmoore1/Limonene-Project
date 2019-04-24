#Setting minimal media

#import cobra and required tools
import cobra
import os
from os.path import join
import matplotlib.pyplot as plt
import pandas as pd
from __future__ import print_function
from cobra import Model, Reaction, Metabolite
import numpy as np
from cobra.medium import minimal_medium

#import iJO1366 - can model
data_dir="C:\Work\ecoli_mod"
model=cobra.io.read_sbml_model(join(data_dir,"eco_can_v1.xml"))
#set biomass as objective function
model.objective="Ec_biomass_iJO1366_core_53p95M"

#minimal medium biomass
cobra.medium.minimal_medium(model)

#minimal medium OAC production
model.objective = 'OAC'
cobra.medium.minimal_medium(model)

#minimal medium CBGAS production
model.objective = 'CBGAS'
cobra.medium.minimal_medium(model)

#set minimal media
data_dir="C:\Work\ecoli_mod"
model=cobra.io.read_sbml_model(join(data_dir,"eco_can_v1.xml"))
model.objective="Ec_biomass_iJO1366_core_53p95M" 

model.reactions.EX_ca2_e.lower_bound = -0.000521
model.reactions.EX_cl_e.lower_bound = -0.000521
model.reactions.EX_cobalt2_e.lower_bound = -0.000002
model.reactions.EX_cu2_e.lower_bound = -0.000071
model.reactions.EX_fe2_e.lower_bound = -0.000825
model.reactions.EX_fe3_e.lower_bound = -0.000781
model.reactions.EX_glc_e.lower_bound = -1.465741
model.reactions.EX_k_e.lower_bound = -0.019519
model.reactions.EX_mg2_e.lower_bound = -0.000867
model.reactions.EX_mn2_e.lower_bound = -0.000069
model.reactions.EX_mobd_e.lower_bound = -0.000013
model.reactions.EX_nh4_e.lower_bound = -1.080082
model.reactions.EX_ni2_e.lower_bound = -0.000032
model.reactions.EX_o2_e.lower_bound = -2.039876
model.reactions.EX_pi_e.lower_bound = -0.096463
model.reactions.EX_so4_e.lower_bound = -0.025221
model.reactions.EX_zn2_e.lower_bound = -0.000034

#save model in minimal media
cobra.io.write_sbml_model(model,filename="C:/Work/ecoli_mod/eco_can_v3.xml")

#Check maximal production for OAC - Define the maximal production of the OAC for the constraints
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_can_v3.xml"))

model.objective = "OAC"
OACo = model.optimize().f 

#Application
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_can_v3.xml"))
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
bm.to_csv("C:/Work/ecoli_mod/OAC_biomass_minimal.csv")
oac.to_csv("C:/Work/ecoli_mod/OAC_biomass_minimal.csv") #of Lim constraints in loop & Optimization of biomass

#Check maximal production for CBGAS - Define the maximal production of the CBGAS for the constraints
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_can_v3.xml"))

model.objective = "CBGAS"
CBGASo = model.optimize().f

#Application
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_can_v3.xml"))
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
bm.to_csv("C:/Work/ecoli_mod/CBGAS_biomass_minimal.csv")
cbgas.to_csv("C:/Work/ecoli_mod/CBGAS_biomass_minimal.csv") #of Lim constraints in loop & Optimization of biomass

