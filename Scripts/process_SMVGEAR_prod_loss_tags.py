#!/usr/bin/python
import sys
from . import AC_tools as AC
"""
Example processing of SMVGEAR prod/loss tags for GEOS-Chem diganotic (ND65). A variety
of functions for working with KPP/SMVGEAR tags are in AC_tools  
(funcs4GESOSC/variables).

NOTES
---
 - This code is for working with smvgear diagnostics. for KPP P/L see 
 Tagged_GC_KPP_Mechanism4family.py
 - details on the GEOS-Chem diagnostic are in the GEOS-Chem manual 
 (http://acmg.seas.harvard.edu/geos/doc/man/chapter_13.html)
"""

# --- Master debug setting
DEBUG = True


def main(trop_limit=True, res='4x5',  debug=False):
    """ 
    Get prod loss output for a family and print this to screen  
    """
    # --- Get family from Command line (and other vars)
    wd = sys.argv[1]
    spec = sys.argv[2]
    # version?
    ver = AC.iGEOSChem_ver(wd)

    # --- Get all tags for this family (through dictionary route)
    # ( e.g. 'PIOx', 'LIOx',  'P_Iy', 'L_Iy' )
    nums, rxns, tags, Coe = AC.prod_loss_4_spec(wd, spec, ver=ver)
    # beatify reaction strings
    rxnstr_l = [''.join(i[4:]) for i in rxns]
    # one consider one tag per reaction and tagged reactions
    try:
        tags = [i[0] for i in tags]  # just consider first tag
    except:
        print('WARNING! - attempting to process just tagged reactions')
        detail_zip = list(zip(rxnstr_l, list(zip(nums, tags))))
        untagged = [n for n, i in enumerate(tags) if (len(i) < 1)]
        print('Untagged reactions: ', [detail_zip[i] for i in untagged])
        tags = [i for n, i in enumerate(tags) if (n not in untagged)]
        tags = [i[0] for i in tags]  # just consider first tag
#        tags.pop( tags.index('LR71') )  # rm tag for ClOO loss...

    # --- Extract prod loss for these tracers
    # get prod loss IDs
    PDs = [AC.PLO3_to_PD(i, ver=ver, wd=wd, fp=True) for i in tags]
    # extract en mass
    fam_loss = AC.get_GC_output(wd, vars=['PORL_L_S__'+i for i in PDs],
                                trop_limit=trop_limit, r_list=True)
#    print [ ( i.shape, i.sum() ) for i in fam_loss ]
    # Get reference species for family ( e.g. so output is in X g of Y )
    ref_spec = AC.get_ref_spec(spec)
    # get shared variable arrrays
    s_area = get_surface_area(res=res)[..., 0]  # m2 land map
    # convert to mass terms  ( in g X )
    fam_loss = convert_molec_cm3_s_2_g_X_s(ars=fam_loss,
                                           ref_spec=ref_spec, wd=wd, conbine_ars=False,
                                           rm_strat=True, month_eq=True)
    print([i.shape for i in fam_loss])

    # sum and convert to Gg
    p_l = [i.sum() / 1E9 for i in fam_loss]

    # --- print output as: reaction, magnitude, percent of family
    pcent = [np.sum(i)/np.sum(p_l)*100 for i in p_l]
    d = dict(list(zip(tags, list(zip(rxnstr_l, p_l, pcent)))))
    df = pd.DataFrame(d).T
    df.columns = ['rxn', 'Gg X', '% of total']
    # sort
    df = df.sort_values(['% of total'], ascending=False)
    print(df)


if __name__ == "__main__":
    main(debug=DEBUG)
