#Single and double gene deletions

#import cobra and required tools
import cobra
import os
from os.path import join
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option('display.max_rows', 2600)
from __future__ import print_function
from cobra import Model, Reaction, Metabolite
import numpy as np
import cobra.test
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)
	
#import iJO1366 - Limo model - complex media
data_dir="C:\Work\ecoli_mod"
model=cobra.io.read_sbml_model(join(data_dir,"eco_Lim_v2.xml"))

#set biomass as objective function
model.objective="Ec_biomass_iJO1366_core_53p95M"

#single gene deletion
G=cobra.flux_analysis.single_gene_deletion(model)
G

#set limonene flux as objective function
model.objective = 'LIMS'

#single gene deletion
model.reactions.Ec_biomass_iJO1366_core_53p95M.upper_bound = 0.1 #constrain biomass
G=cobra.flux_analysis.single_gene_deletion(model)
G

#import iJO1366 - Limo model - optimised media
data_dir="C:\Work\ecoli_mod"
model=cobra.io.read_sbml_model(join(data_dir,"eco_Lim_v5.xml"))

#set biomass as objective function
model.objective="Ec_biomass_iJO1366_core_53p95M"

#single gene deletion
G=cobra.flux_analysis.single_gene_deletion(model)
G

#set limonene flux as objective function
model.objective = 'LIMS'

#single gene deletion
model.reactions.Ec_biomass_iJO1366_core_53p95M.upper_bound = 0.1 #constrain biomass
G=cobra.flux_analysis.single_gene_deletion(model)
G

#import iJO1366 - Limo model - minimal media
data_dir="C:\Work\ecoli_mod"
model=cobra.io.read_sbml_model(join(data_dir,"eco_Lim_v7.xml"))

#set biomass as objective function
model.objective="Ec_biomass_iJO1366_core_53p95M"

#single gene deletion
G=cobra.flux_analysis.single_gene_deletion(model)
G

#set limonene flux as objective function
model.objective = 'LIMS'

#single gene deletion
model.reactions.Ec_biomass_iJO1366_core_53p95M.upper_bound = 0.1 #constrain biomass
G=cobra.flux_analysis.single_gene_deletion(model)
G

#double gene deletion script - for limonene production - minimal media
data_dir="C:\Work\ecoli_mod"
model=cobra.io.read_sbml_model(join(data_dir,"eco_Lim_v7.xml"))
model.reactions.Ec_biomass_iJO1366_core_53p95M.upper_bound = 0.1 #constrain biomass
model.objective = 'LIMS'
G=cobra.flux_analysis.double_gene_deletion(model)
G.to_csv("C:/Work/ecoli_mod/double_gene_deletions_min.csv")

#double gene deletion script - for limonene production - complex media
data_dir="C:\Work\ecoli_mod"
model=cobra.io.read_sbml_model(join(data_dir,"eco_Lim_v2.xml"))
model.reactions.Ec_biomass_iJO1366_core_53p95M.upper_bound = 0.1 #constrain biomass
model.objective = 'LIMS'
G=cobra.flux_analysis.double_gene_deletion(model)
G.to_csv("C:/Work/ecoli_mod/double_gene_deletions_min.csv")

#double gene deletion script - for limonene production - optimised media
data_dir="C:\Work\ecoli_mod"
model=cobra.io.read_sbml_model(join(data_dir,"eco_Lim_v5.xml"))
model.reactions.Ec_biomass_iJO1366_core_53p95M.upper_bound = 0.1 #constrain biomass
model.objective = 'LIMS'
G=cobra.flux_analysis.double_gene_deletion(model)
G.to_csv("C:/Work/ecoli_mod/double_gene_deletions_min.csv")