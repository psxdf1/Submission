'''
Script to create a base community model using AG and GD that can then be modified further.
Uses the micom package's community class that requires construction from a pandas dataframe.
'''
import micom as mc
import pandas as pd

def main():
 
    # AG_uncoupled_filepath = '../../data/AG.xml'
    AG_uncoupled_filepath = 'project/data/AG.xml'
    GD_uncoupled_filepath = 'project/data/GD.xml'

    # Id Labels
    ids  = ["AG_uc","GD_uc"] # _uc to indicate uncoupled 'wt' strains
    taxa = pd.DataFrame({"id":ids}) 
    taxa["genus"] = "ecoli"
    taxa["species"] = ["Ecoli AG", "Ecoli GD"]
    taxa ["file"] = [AG_uncoupled_filepath,GD_uncoupled_filepath]

    com = mc.Community(taxa)

    # Adjustments 
    # Remove from media
    com.reactions.EX_Tyrosol_m.lower_bound = 0 
    com.reactions.EX_phe__L_m.lower_bound = 0 
    com.reactions.EX_tyr__L_m.lower_bound = 0

    # Additionally cut off glucose from AG, these reactions are reset wheb the community is created.?
    com.reactions.GLCtex_copy1__AG_uc.bounds = (0,0)
    com.reactions.GLCtex_copy2__AG_uc.bounds = (0,0)

    com.to_pickle("project/data/community.pickle")
    return com

def display_community_uptake_bounds(com):
    # Quick script to display the uptake bounds of the community as well as the internal constituents.
    from pandas import DataFrame
    com.medium

    tab = []
    for rxn_ in com.medium.keys():
        # print(rxn[:-1])
        rxn = rxn_[:-1]
        c = 0
        ag_bound, gd_bound = 9999, 9999 # Izf it comes out as 999 then there is a problem
        for i in com.reactions:

            if rxn == i.id[:-8]:
                if i.id[-5:-3] == "AG":
                    ag_bound = i.lower_bound
                elif i.id[-5:-3] == "GD":
                    gd_bound = i.lower_bound
                c+=1
                if c ==2:
                    tab.append(
                        {
                            'Reaction' : rxn[:-1],
                            'Community' : -com.medium[rxn_],
                            'AG' : ag_bound,
                            'GD' : gd_bound
                        }
                    )

    return DataFrame(tab)

if __name__ == "__main__":
    com = main()
    print(display_community_uptake_bounds(com))

    # just confirm no uptake of xylose - the uptake bounds are misleading
    # for i in com.metabolites.xyl__D_p__GD_uc.reactions:
    #     print(i, i.bounds)
