#Limonene pathway introduction

#import cobra and required tools
import cobra
import os
from os.path import join
import matplotlib.pyplot as plt
import pandas as pd
from __future__ import print_function
from cobra import Model, Reaction, Metabolite
import numpy as np

#import iJO1366 model
data_dir="C:\Work\ecoli_mod"
model=cobra.io.read_sbml_model(join(data_dir,"iJO1366.xml"))
#set biomass as objective function
model.objective="Ec_biomass_iJO1366_core_53p95M"
model.optimize().f

#Addition of mevalonate pathway

#Setting reversibility of reaction 1
model.reactions.get_by_id("GPPS").lower_bound=-1000
model.reactions.get_by_id("GPPS").reversibility

#Reaction 2 of MEV pathway - HMGCOAS
reaction = Reaction('HMGCOAS')
reaction.name = 'Hydroxymethyl glutaryl-CoA synthase'
reaction.subsystem = 'Limonene synthesis'
reaction.lower_bound = 0
reaction.upper_bound = 1000

    #Add metabolite HMG Coa
HMGcoa_c = Metabolite(
    'HMGcoa_c',
    formula='C27H40N7O20P3S',
    name='Hydroxymethyl glutaryl-CoA',
    compartment='c')
accoa_c = model.metabolites.accoa_c
aacoa_c = model.metabolites.aacoa_c
coa_c = model.metabolites.coa_c
h2o_c = model.metabolites.h2o_c


    
    #Add reaction
reaction.add_metabolites({
    accoa_c:-1.0,
    aacoa_c:-1.0,
    h2o_c: -1.0,
    HMGcoa_c:1.0,
    coa_c:1.0,
    
})

reaction.gene_reaction_rule='(SA_RS13375)'
reaction.genes

model.add_reactions([reaction])

#Reaction 3 of MEV pathway - HMGCOAR
reaction = Reaction('HMGCOAR')
reaction.name = 'Hydroxymethyl glutaryl-CoA reductase'
reaction.subsystem = 'Limonene synthesis'
reaction.lower_bound = 0
reaction.upper_bound = 1000

    #Add metabolite mevalonate
mev__R_c = Metabolite(
    'mev__R_c',
    formula='C6H11O4',
    name='Mevalonate',
    compartment='c')
HMGcoa_c = model.metabolites.HMGcoa_c
nadph_c = model.metabolites.nadph_c
coa_c = model.metabolites.coa_c
nadp_c = model.metabolites.nadp_c
h_c = model.metabolites.h_c

    #Add reaction
reaction.add_metabolites({
    HMGcoa_c:-1.0,
    nadph_c:-2.0,
    h_c:-1.0,
    mev__R_c:1.0,
    coa_c:1.0,
    nadp_c:2.0,
    
})

reaction.gene_reaction_rule='(SA_RS13370)'
reaction.genes

model.add_reactions([reaction])

#Reaction 4 of MEV pathway - MEVK1
reaction = Reaction('MEVK1')
reaction.name = 'Mevalonate kinase'
reaction.subsystem = 'Limonene synthesis'
reaction.lower_bound = 0
reaction.upper_bound = 1000

    #Add metabolite mevalonate
pmev_c = Metabolite(
    'pmev_c',
    formula='C6H10O7P',
    name='Phosphomevalonate',
    compartment='c')
mev__R_c = model.metabolites.mev__R_c
atp_c = model.metabolites.atp_c
adp_c = model.metabolites.adp_c
h_c = model.metabolites.h_c

    #Add reaction
reaction.add_metabolites({
    mev__R_c:-1.0,
    atp_c:-1.0,
    pmev_c:1.0,
    adp_c:1.0,
    h_c:1.0,
   
})

reaction.gene_reaction_rule='(YMR208W)'
reaction.genes

model.add_reactions([reaction])

#Reaction 5 of MEV pathway - PMEVK
reaction = Reaction('PMEVK')
reaction.name = 'Phosphomevalonate kinase'
reaction.subsystem = 'Limonene synthesis'
reaction.lower_bound = -1000
reaction.upper_bound = 1000

    #Add metabolite MVAPP
dpmev_c = Metabolite(
    'dpmev_c',
    formula='C6H10O10P2',
    name='Diphosphomevalonate',
    compartment='c')
pmev_c = model.metabolites.pmev_c
atp_c = model.metabolites.atp_c
adp_c = model.metabolites.adp_c

    #Add reaction
reaction.add_metabolites({
    pmev_c:-1.0,
    atp_c:-1.0,
    dpmev_c:1.0,
    adp_c:1.0,
   
})

reaction.gene_reaction_rule='(YMR220W)'
reaction.genes

model.add_reactions([reaction])

#Reaction 6 of MEV pathway - DPMVD
reaction = Reaction('DPMVD')
reaction.name = 'Mevalonate decarboxylase'
reaction.subsystem = 'Limonene synthesis'
reaction.lower_bound = 0
reaction.upper_bound = 1000

dpmev_c = model.metabolites.dpmev_c
atp_c = model.metabolites.atp_c
ipdp_c = model.metabolites.ipdp_c
adp_c = model.metabolites.adp_c
pi_c = model.metabolites.pi_c
co2_c = model.metabolites.co2_c

    #Add reaction
reaction.add_metabolites({
    dpmev_c:-1.0,
    atp_c:-1.0,
    ipdp_c:1.0,
    adp_c:1.0,
    pi_c:1.0,
    co2_c:1.0,
   
})

reaction.gene_reaction_rule='(YNR043W)'
reaction.genes

model.add_reactions([reaction])

#Adding limonene synthase reaction
reaction = Reaction('LIMS')
reaction.name = 'Limonene synthase'
reaction.subsystem = 'Limonene synthesis'
reaction.lower_bound = 0
reaction.upper_bound = 1000

    #Add metabolite limonene
lim_c = Metabolite(
    'lim_c',
    formula='C10H16',
    name='Limonene',
    compartment='c')
grdp_c = model.metabolites.grdp_c
ppi_c = model.metabolites.ppi_c

    
    #Add reaction
reaction.add_metabolites({
    grdp_c:-1.0,
    lim_c:1.0,
    ppi_c:1.0,
   
})

reaction.gene_reaction_rule='(lim1)'
reaction.genes

model.add_reactions([reaction])

#Add export reaction
reaction = Reaction('EX_LIM')
reaction.name = 'Exchange Limonene'
reaction.subsystem = 'Exchange'
reaction.lower_bound = 0
reaction.upper_bound = 1000

    #Add metabolite limonene
lim_e = Metabolite(
    'lim_e',
    formula='C10H16',
    name='Limonene',
    compartment='e')
    
    #Add reaction
reaction.add_metabolites({
    lim_e:1.0  
})

model.add_reactions([reaction])

#Add sink reaction
reaction = Reaction('DM_LIM')
reaction.name = 'Sink Limonene'
reaction.subsystem = 'Exchange'
reaction.lower_bound = 0
reaction.upper_bound = 1000

    #Add metabolite limonene
lim_c = model.metabolites.lim_c

    
    #Add reaction
reaction.add_metabolites({
    lim_c:-1.0
})

model.add_reactions([reaction])

#Save adapted model
cobra.io.write_sbml_model(model,filename="C:/Work/ecoli_mod/eco_Lim_v1.xml")

#Check maximal production for limonene - Define the maximal production of the LIM for the constraints
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_Lim_v1.xml"))

model.objective = "LIMS"
LIMo = model.optimize().f  

### Application
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_Lim_v1.xml"))
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
Lim.to_csv("C:/Work/ecoli_mod/LIM_biomass.csv") #of Lim constraints in loop & Optimization of biomass
