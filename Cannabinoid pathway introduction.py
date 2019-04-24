#Cannabinoid pathway introduction

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
reaction.subsystem = 'Cannabinoid production'
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
reaction.subsystem = 'Cannabinoid production'
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
reaction.subsystem = 'Cannabinoid production'
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
reaction.subsystem = 'Cannabinoid production'
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

#Reaction 6 of MEV pathway
reaction = Reaction('DPMVD')
reaction.name = 'Mevalonate decarboxylase'
reaction.subsystem = 'Cannabinoid production'
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

#Addition of fatty acid biosynthesis pathway

#Reaction 1 of fatty acid biosynthesis pathway - TSK reaction 1
reaction = Reaction('TKS1')
reaction.name = 'TKS1 reaction 1'
reaction.subsystem = 'Cannabinoid production'
reaction.lower_bound = 0
reaction.upper_bound = 1000

#metabolites 
malcoa_c = Metabolite(
    'malcoa_c',
    formula='C24H33N7O19P3S',
    name='Malonyl-CoA',
    compartment='c')
hxcoa_c = Metabolite(
    'hxcoa_c',
    formula='C27H42N7O17P3S',
    name='Hexanoyl-CoA',
    compartment='c')
coa_c = Metabolite(
    'coa_c',
    formula='C21H32N7O16P3S',
    name='Coenzyme A',
    compartment='c')
co2_c = Metabolite(
    'co2_c',
    formula='CO2',
    name='CO2',
    compartment='c')
oocoa_c = Metabolite(
    '3oocoa_c',
    formula='C29H44N7O18P3S',
    name='3-Oxooctanoyl-CoA',
    compartment='c')
h_c = Metabolite(
    'h_c',
    formula='H',
    name='H+',
    compartment='c')

    #Add reaction
reaction.add_metabolites({
    malcoa_c:-1.0,
    hxcoa_c:-1.0,
    h_c:-1.0,
    coa_c:1.0,
    co2_c:1.0,
    oocoa_c:1.0,
    
})

reaction.gene_reaction_rule='(TKS_G)'
reaction.genes

model.add_reactions([reaction])

#Reaction 2 of fatty acid biosynthesis pathway - TSK1 reaction 2
reaction = Reaction('TKS2')
reaction.name = 'TKS1 reaction 2'
reaction.subsystem = 'Cannabinoid production'
reaction.lower_bound = 0
reaction.upper_bound = 1000

#metabolites
docoa_c = Metabolite(
    '35docoa_c',
    formula='C31H46N7O19P3S',
    name='3,5-dioxodecanoyl-CoA',
    compartment='c',)
oocoa_c = Metabolite(
    '3oocoa_c',
    formula='C29H44N7O18P3S',
    name='3-Oxooctanoyl-CoA',
    compartment='c')
malcoa_c = Metabolite(
    'malcoa_c',
    formula='C24H33N7O19P3S',
    name='Malonyl-CoA',
    compartment='c')
coa_c = Metabolite(
    'coa_c',
    formula='C21H32N7O16P3S',
    name='Coenzyme A',
    compartment='c')
co2_c = Metabolite(
    'co2_c',
    formula='CO2',
    name='CO2',
    compartment='c')
h_c = Metabolite(
    'h_c',
    formula='H',
    name='H+',
    compartment='c')

    #Add reaction
reaction.add_metabolites({
    malcoa_c:-1.0,
    oocoa_c:-1.0,
    h_c:-1.0,
    coa_c:1.0,
    co2_c:1.0,
    docoa_c:1.0,
   
})

reaction.gene_reaction_rule='(TKS_G)'
reaction.genes

model.add_reactions([reaction])

#Reaction 3 of fatty acid biosynthesis pathway - TSK1 reaction 3
reaction = Reaction('TKS3')
reaction.name = 'TKS1 reaction 3'
reaction.subsystem = 'Cannabinoid production'
reaction.lower_bound = 0
reaction.upper_bound = 1000

    #Add metabolite 35tocoa
tocoa_c = Metabolite(
    '357tocoa_c',
    formula='C33H52N7O20P3S',
    name='3,5,7-trioxododecanoyl-CoA',
    compartment='c')
docoa_c = Metabolite(
    '35docoa_c',
    formula='C31H46N7O19P3S',
    name='3,5-dioxodecanoyl-CoA',
    compartment='c')
malcoa_c = Metabolite(
    'malcoa_c',
    formula='C24H33N7O19P3S',
    name='Malonyl-CoA',
    compartment='c')
coa_c = Metabolite(
    'coa_c',
    formula='C21H32N7O16P3S',
    name='Coenzyme A',
    compartment='c')
co2_c = Metabolite(
    'co2_c',
    formula='CO2',
    name='CO2',
    compartment='c')
h_c = Metabolite(
    'h_c',
    formula='H',
    name='H+',
    compartment='c')

    #Add reaction
reaction.add_metabolites({
    malcoa_c:-1.0,
    docoa_c:-1.0,
    h_c:-5.0,
    coa_c:1.0,
    co2_c:1.0,
    tocoa_c:1.0,
   
})

reaction.gene_reaction_rule='(TKS_G)'
reaction.genes

model.add_reactions([reaction])

#3oocoa hydrolysis reaction
reaction = Reaction('oocoaH')
reaction.name = '3oocoa hydrolysis'
reaction.subsystem = 'Cannabinoid production'
reaction.lower_bound = -1000
reaction.upper_bound = 1000

#metabolites
pdal_c = Metabolite(
    'pdal_c',
    formula='C10H14O3',
    name='Pentyl diacetic acid lactone',
    compartment='c')
oocoa_c = Metabolite(
    '3oocoa_c',
    formula='C29H44N7O18P3S',
    name='3-Oxooctanoyl-CoA',
    compartment='c')
h2o_c = Metabolite(
    'h2o_c',
    formula='H2O',
    name='H2O',
    compartment='c')
coa_c = Metabolite(
    'coa_c',
    formula='C21H32N7O16P3S',
    name='Coenzyme A',
    compartment='c')
co2_c = Metabolite(
    'co2_c',
    formula='CO2',
    name='CO2',
    compartment='c')
o2_c = Metabolite(
    'o2_c',
    formula='O2',
    name='O2',
    compartment='c')

    #Add reaction
reaction.add_metabolites({
    oocoa_c:-1.0,
    h2o_c:-1.0,
    co2_c:-2.0,
    pdal_c:1.0,
    coa_c:1.0, 
    o2_c:2.0,
   
})

model.add_reactions([reaction])

#35docoa hydrolysis reaction
reaction = Reaction('docoaH')
reaction.name = '35docoa hydrolysis'
reaction.subsystem = 'Cannabinoid production'
reaction.lower_bound = -1000
reaction.upper_bound = 1000

#Add metabolites
htal_c = Metabolite(
    'htal_c',
    formula='C12H16O4',
    name='Hexanoyl triacetic acid lactone',
    compartment='c')
docoa_c = Metabolite(
    '35docoa_c',
    formula='C31H46N7O19P3S',
    name='3,5-dioxodecanoyl-CoA',
    compartment='c')
h2o_c = Metabolite(
    'h2o_c',
    formula='H2O',
    name='H2O',
    compartment='c')
co2_c = Metabolite(
    'co2_c',
    formula='CO2',
    name='CO2',
    compartment='c')
coa_c = Metabolite(
    'coa_c',
    formula='C21H32N7O16P3S',
    name='Coenzyme A',
    compartment='c')
o2_c = Metabolite(
    'o2_c',
    formula='O2',
    name='O2',
    compartment='c')

    #Add reaction
reaction.add_metabolites({
    docoa_c:-1.0,
    h2o_c:-1.0,
    co2_c: -2.0,
    htal_c:1.0,
    coa_c:1.0,
    o2_c:2.0,
      
})

model.add_reactions([reaction])

#Olivetol synthesis reaction
reaction = Reaction('oliv1')
reaction.name = 'Olivetol synthesis'
reaction.subsystem = 'Cannabinoid production'
reaction.lower_bound = -1000
reaction.upper_bound = 1000

    #Add metabolite Olivetol
oliv_c = Metabolite(
    'oliv_c',
    formula='C11H16O2',
    name='Olivetol',
    compartment='c')
tocoa_c = Metabolite(
    '357tocoa_c',
    formula='C33H52N7O20P3S',
    name='3,5,7-trioxododecanoyl-CoA',
    compartment='c')
coa_c = Metabolite(
    'coa_c',
    formula='C21H32N7O16P3S',
    name='Coenzyme A',
    compartment='c')
co2_c = Metabolite(
    'co2_c',
    formula='CO2',
    name='CO2',
    compartment='c')
h_c = Metabolite(
    'h_c',
    formula='H',
    name='H+',
    compartment='c')

    #Add reaction
reaction.add_metabolites({
    tocoa_c:-1.0,
    oliv_c:1.0,
    coa_c:1.0,
    co2_c:1.0,
    h_c:4.0,
         
})

reaction.gene_reaction_rule='(TKS_G)'
reaction.genes

model.add_reactions([reaction])

#Olivetolic acid cyclase reaction
reaction = Reaction('OAC')
reaction.name = 'Olivetolic acid cyclase'
reaction.subsystem = 'Cannabinoid production'
reaction.lower_bound = 0
reaction.upper_bound = 1000

    #Add metabolites
ota_c = Metabolite(
    'ota_c',
    formula='C12H16O4',
    name='Olivetolic acid',
    compartment='c')
tocoa_c = Metabolite(
    '357tocoa_c',
    formula='C33H52N7O20P3S',
    name='3,5,7-trioxododecanoyl-CoA',
    compartment='c')
coa_c = Metabolite(
    'coa_c',
    formula='C21H32N7O16P3S',
    name='Coenzyme A',
    compartment='c')
h_c = Metabolite(
    'h_c',
    formula='H',
    name='H+',
    compartment='c')

     #Add reaction
reaction.add_metabolites({
    tocoa_c:-1.0,
    ota_c:1.0,
    coa_c:1.0,
    h_c:4.0,
       
})

reaction.gene_reaction_rule='(OAC_G)'
reaction.genes

model.add_reactions([reaction])

#Cannabigerolic acid synthesis production
reaction = Reaction('CBGAS')
reaction.name = 'Cannabigerolic acid synthase'
reaction.subsystem = 'Cannabinoid production'
reaction.lower_bound = 0
reaction.upper_bound = 1000

    #Add metabolites
cbga_c = Metabolite(
    'cbga_c',
    formula='C22H32O4',
    name='Cannabinogerolic acid',
    compartment='c')
ota_c = Metabolite(
    'ota_c',
    formula='C12H16O4',
    name='Olivetolic acid',
    compartment='c')
grdp_c = Metabolite(
    'grdp_c',
    formula='C10H17O7P2',
    name='Geranyl diphosphate',
    compartment='c')
ppi_c = Metabolite(
    'ppi_c',
    formula='HO7P2',
    name='Diphosphate',
    compartment='c')

    #Add reaction
reaction.add_metabolites({
    ota_c:-1.0,
    grdp_c:-1.0,
    cbga_c:1.0,
    ppi_c:1.0,
   
})

reaction.gene_reaction_rule='(CBGAS_G)'
reaction.genes

model.add_reactions([reaction])

#THCA synthase reaction
reaction = Reaction('THCAS')
reaction.name = 'Tetrahydrocannabinolic acid synthase'
reaction.subsystem = 'Cannabinoid production'
reaction.lower_bound = 0
reaction.upper_bound = 1000

    #Add metabolite tetrahydrocannabinolic acid
thca_c = Metabolite(
    'thca_c',
    formula='C22H30O4',
    name='Tetrahydrocannabinolic acid',
    compartment='c')
cbga_c = Metabolite(
    'cbga_c',
    formula='C22H32O4',
    name='Cannabinogerolic acid',
    compartment='c')
o2_c = Metabolite(
    'o2_c',
    formula='O2',
    name='O2',
    compartment='c')
h2o2_c = Metabolite(
    'h2o2_c',
    formula='H2O2',
    name='Hydrogen peroxide',
    compartment='c')

    #Add reaction
reaction.add_metabolites({
    cbga_c:-1.0,
    o2_c:-1.0,
    thca_c:1.0,
    h2o2_c:1.0,
      
})

reaction.gene_reaction_rule='(THCAS_G)'
reaction.genes

model.add_reactions([reaction])

#Cannabidiolic acid synthase reaction
reaction = Reaction('CBDAS')
reaction.name = 'Cannabidiolic acid synthase'
reaction.subsystem = 'Cannabinoid production'
reaction.lower_bound = 0
reaction.upper_bound = 1000

    #Add metabolite cannabidiolic acid
cbda_c = Metabolite(
    'cbda_c',
    formula='C22H30O4',
    name='Cannabidiolic acid',
    compartment='c')
o2_c = Metabolite(
    'o2_c',
    formula='O2',
    name='O2',
    compartment='c')
cbga_c = Metabolite(
    'cbga_c',
    formula='C22H32O4',
    name='Cannabinogerolic acid',
    compartment='c')
h2o2_c = Metabolite(
    'h2o2_c',
    formula='H2O2',
    name='Hydrogen peroxide',
    compartment='c')

    #Add reaction
reaction.add_metabolites({
    cbga_c:-1.0,
    o2_c:-1.0,
    cbda_c:1.0,
    h2o2_c:1.0,
   
})

reaction.gene_reaction_rule='(CBDAS_G)'
reaction.genes

model.add_reactions([reaction])

#Reaction export THCA

reaction = Reaction('EX_THCA')
reaction.name = 'Exchange Tetrahydrocannabinolic acid'
reaction.subsystem = 'Exchange'
reaction.lower_bound = 0
reaction.upper_bound = 1000

    #Add metabolite 
thca_e = Metabolite(
    'thca_e',
    formula='C22H30O4',
    name='Tetrahydrocannabinolic acid',
    compartment='e')
    
    #Add reaction
reaction.add_metabolites({
    thca_e:1.0  
})

model.add_reactions([reaction])

#Reaction sink THCA

reaction = Reaction('DM_THCA')
reaction.name = 'Sink Tetrahydrocannabinolic acid'
reaction.subsystem = 'Exchange'
reaction.lower_bound = 0
reaction.upper_bound = 1000

    #Add metabolite 
thca_c = model.metabolites.thca_c

    
    #Add reaction
reaction.add_metabolites({
    thca_c:-1.0
})

model.add_reactions([reaction])

#Reaction export CBDA

reaction = Reaction('EX_CBDA')
reaction.name = 'Exchange Cannabidiolic acid'
reaction.subsystem = 'Exchange'
reaction.lower_bound = 0
reaction.upper_bound = 1000

    #Add metabolite 
cbda_e = Metabolite(
    'cbda_e',
    formula='C22H30O4',
    name='Cannabidiolic acid',
    compartment='e')
    
    #Add reaction
reaction.add_metabolites({
    cbda_e:1.0  
})

model.add_reactions([reaction])

#Reaction sink CBDA

reaction = Reaction('DM_CBDA')
reaction.name = 'Sink Cannabidiolic acid'
reaction.subsystem = 'Exchange'
reaction.lower_bound = 0
reaction.upper_bound = 1000

    #Add metabolite 
cbda_c = model.metabolites.cbda_c

    
    #Add reaction
reaction.add_metabolites({
    cbda_c:-1.0
})

model.add_reactions([reaction])

#save adapted model
cobra.io.write_sbml_model(model,filename="C:/Work/ecoli_mod/eco_can_v1.xml")

# Check maximal production for olivetolic acid - Define the maximal production of OTA for the constraints
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_can_v1.xml"))

model.objective = "OAC"
OACo = model.optimize().f

#Application
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_can_v1.xml"))
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
bm.to_csv("C:/Work/ecoli_mod/OAC_biomass.csv")
oac.to_csv("C:/Work/ecoli_mod/OAC_biomass.csv") #of Lim constraints in loop & Optimization of biomass

# Check maximal production for CBGAS - Define the maximal production of CBGAS for the constraints
model.objective = "CBGAS"
CBGASo = model.optimize().f

#Application
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_can_v1.xml"))
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
bm.to_csv("C:/Work/ecoli_mod/CBGAS_biomass.csv")
cbgas.to_csv("C:/Work/ecoli_mod/CBGAS_biomass.csv") #of Lim constraints in loop & Optimization of biomass


