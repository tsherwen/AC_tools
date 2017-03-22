#!/usr/bin/python  
"""
Example processing of KPP prod/loss tags for GEOS-Chem diganotic (ND65). A variety
of functions for working with KPP/SMVGEAR tags are also in AC_tools  
(funcs4GESOSC/funcs_vars).

NOTES
---
 - This code is for working with smvgear diagnostics. for KPP P/L see 
 process_SMVGEAR_prod_loss_tags.py
 - details on the GEOS-Chem diagnostic are in the GEOS-Chem manual 
 (http://acmg.seas.harvard.edu/geos/doc/man/chapter_13.html)
"""
# import modules
import sys
import AC_tools as AC


def main( fam='LOx'):
    """
    Demonstate some processing of KPP tags using functions in AC_tools
    """
      
    
    # ---------  Processing to setup input files for KPP ------

    # --- Prep input for gckpp.kpp
    # Takes a list of reaction numbers or (fam) and prints out 
    # prod/loss in a form that can be copy/pasted into KPP (gckpp.kpp)
    #
    # prt_families4rxns_to_input_to_PROD_LOSS()
        
    # --- Prep input for globchem.spc
    #


    # prt_families4rxns_to_input_to_PROD_LOSS_globchem_spec()


    # ---------  Processing to look at output from KPP ------
    # --- Get a dictionary of all reactions in mechanism by number. 
    # 
#    get_dict_of_KPP_mech()

    # --- Get a dictionary of all tagged reactions. 
    #
    # get_dictionary_of_tagged_reactions()


    # --- Get tag reactions for family. 
    #
    # get_KPP_tagged_rxns(fam=fam)    
    # or function below ... double up?
    # get_tags4_family()

    # --- Get stoichiometry for tags
    #
    # get_stioch_for_family_reactions()
    
    # --- Assign OX family based on reaction component. 
    #
    # 
    # get_Ox_family_tag_based_on_reactants()

    # ---  Get tags in given lists of reaction numbers. 
    #
    #
    # get_tags_in_rxn_numbers()

    
    # ---------  Processing to extract data and process model p/l data ------    
    # --- Extract some data from ND65 [molec/cm3/s] - 	PORL-L=$
    #
    #
    # From manual: "Chemical family P/L rates. 
    # See the FlexChem wiki page for more information."
    #
    # get dictionary of deafuly variable settings
#    Var_rc = AC.get_default_variable_dict()

    # tags?
#    PL_tags =['LOx', 'POx']
#    PL_tags = ['P'+i for i in tags ]
    # extract
#    ars = AC.get_GC_output( wd=Var_rc['wd'], r_list=True, \
#        vars=['PORL_L_S__'+i for i in PL_tags ], trop_limit=Var_rc['trop_limit'])

    # --- Convert prod/loss [molec/cm3/s] to mass per grid box [ g(X)/s ]
    #    
    # AC.process_to_X_per_s()
    

if __name__ == "__main__":
    main( )
    