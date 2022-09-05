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


if __name__ == "__main__":
    main()

