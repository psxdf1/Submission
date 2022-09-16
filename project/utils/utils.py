from ast import Raise


def foo():
    return "World, Hello"


def set_bound_modifications(model,file):
    '''
    Takes a csv file containing a list of rxn bounds and applies to the model.
    
    Check that all rxns are in the model & that it can optimze afterward.

    model = set_bound_modifications(model, file_path)
    '''
    # use context to make changes and checks
    from csv import reader

    with open(file, "r") as my_file:
        # pass the file object to reader()
        file_reader = reader(my_file)
        # do this for all the rows
        for line in file_reader:
            rxn,lb,ub = line[0], line[1], line[2]
            try:
                rxn = [i for i in model.reactions if i.id == rxn][0]
                if len(lb) > 0:
                    rxn.lower_bound = float(lb) 
                if len(ub) > 0:
                    rxn.upper_bound = float(ub)
            except Exception as e:
                print(e)
                print("Case: ",file,"Model might not contain ", line[0], lb, ub, rxn)
                print(line, len(line))
                print(len(lb),not lb)
                return False

    return model

def runcases(model,path_to_cases, tgt, num_of_files=3):
    # Aimed at testing different internal flux bound cases in individual models
    i = 0
    m_gr = []
    m_tmy = []

    for i in range(num_of_files):
        fn = str(i) + ".csv"
        with model as m: 
            m = set_bound_modifications(m,f"{path_to_cases}/{fn}")
            opt_bm = m.slim_optimize()
            m_gr.append(opt_bm)

            # Constraint the biomass flux to the computed maximum
            # insert AG biomass flux constraint here
            try:
                model.reactions.BIOMASS_Ec_iHK1487_core.lower_bound = opt_bm
            except:
                try:
                    model.reactions.BIOMASS_Ec_iJO1366_WT_53p95M = opt_bm
                except AttributeError:
                    raise 
            # insert GD biomass flux constriant here
            m.objective = tgt
            m_tmy.append(m.slim_optimize())

    return m_gr, m_tmy



def runcases_community(model,path_to_cases, tgt, num_of_files=3):
    # Aimed at testing different internal flux bound cases in community  models
    i = 0
    m_gr = []
    m_tmy = []

    for i in range(num_of_files):
        fn = str(i) + ".csv"
        with model as m: 
            m = set_bound_modifications(m,f"{path_to_cases}/{fn}")

            temp = model.optimize_all()
            opt_bm_AG = temp.AG_uc
            opt_bm_GD = temp.GD_uc
            print("AG BM optimal",opt_bm_AG)
            print("GD BM optimal",opt_bm_GD)
            # # opt_bm = m.slim_optimize()
            # m_gr.append("x")

            model.reactions.BIOMASS_Ec_iJO1366_core_53p95M__AG_uc.lower_bound = 0.818
            model.reactions.BIOMASS_Ec_iHK1487_core__GD_uc.lower_bound = 0.829
            # # Constrain the biomass flux to the computed maximum
            # try:
                # model.reactions.BIOMASS_Ec_iJO1366_core_53p95M__AG_uc.lower_bound = 0.5 * opt_bm_AG
                # model.reactions.BIOMASS_Ec_iHK1487_core__GD_uc.lower_bound = 0.5 * opt_bm_GD

            # except AttributeError:
            #     raise

            m.objective = tgt
            m_tmy.append(m.slim_optimize())

    return m_gr, m_tmy

def pe_data(model_,tgt_str):
    '''
    Function to plot production envelopes
    Assumes that the model's objective function is biomass
    '''
    import numpy as np

    with model_ as model:
        bm = [i for i in model.reactions if i.objective_coefficient > 0][0]
        # print("Original objective", bm.id)

        maximum_growth = model.slim_optimize()
        model.objective = tgt_str
        # print("new objective: ", model.objective.expression)
        pin, res_max, res_min = [],[],[]

        for i in np.arange(0.0,1.01,0.05):
            bm.lower_bound = maximum_growth*i
            # pin.append(i) # For % of maximum
            pin.append(maximum_growth*i) # For growth rate [h-1] 
            model.objective_direction = 'min'
            res_min.append((model.slim_optimize()))
            model.objective_direction = 'max'
            res_max.append((model.slim_optimize()))
        
        # model.objective.lower_bound = 0 # tidy up - not working
        # bm.lower_bound = 0 # sets the lower bound of the biomass rxn back to 0 
        # model.objective = bm.id # Restores biomass rxn back to the objective

    return pin, res_max, res_min, tgt_str


def pe_plot(pin, res_max, res_min, tgt_rxn):
    '''
    Plot production envelope
    '''
    import matplotlib
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots() 
    ax.plot(pin, res_max)
    ax.plot(pin, res_min)
    ax.set_ylabel(f"{tgt_rxn} Production Rate; mmol/gDW.hr")
    ax.set_xlabel("Growth Rate [h¯¹]") # (% of maximum)
    ax.set_title(f"{tgt_rxn} Production Envelope")

def pe(model_,tgt_str):
    a,b,c,d = pe_data(model_, tgt_str)
    pe_plot(a,b,c,d)

def getInternalExchangeFluxes(com, sol):
    '''
    Takes a communty model and a solution object derived from that model and
    returns a trimmed dataframe of the non-zero internal exchange metabolite fluxes. 
    This function is not generic - only use with AG_GD derived communities.  
    '''
    int_ex = [i.id for i in com.internal_exchanges]
    int_ex = [i[:-7] for i in int_ex]

    x = sol.fluxes.T

    x = x[x.index.isin(int_ex)].sort_values(by=['GD_uc'], ascending = True)
    
    x = x.drop(columns = ["medium"]) # Remove superflous medium column
    x = x.fillna(0)
    x = x[((x.AG_uc != 0) | (x.GD_uc != 0))] # reduce to a meaningul subset of exchange reactions by removing NAN's & zeros etc
    # x = x[((x.AG_uc != 0) | (x.GD_uc != 0)) & (abs(x.AG_uc) + abs(x.GD_uc) > 0)] # reduce to a meaningul subset of exchange reactions by removing NAN's & zeros etc
    
    return x


def diff_internalExFluxes(com,A,B):
    '''
    Compares the fluxes of two cooperative tradeoff solutions from MICOM and outputs the difference between 
    internal exhange reaction fluxes in the two input scenarios. 

    Note:  Has been hardcoded to work with AG/GD models, rework if use any other. 

    Args:
        com: The micom/cobrapy community model underlying the scenarios, required for the exchange fluxe id's
        A: Micom cooperative tradeoff solution object.
        B: Micom cooperative tradeoff solution object.

    Returns:
        x: Dataframe containing both the flux values for the internal exchanges of both scenarios for each strain, 
        but also the difference between scenarios for each strain.

    Raises:
        None
    '''
    x = A.fluxes.T.compare(B.fluxes.T) # Handy function this
    x['AG_diff'] = x.AG_uc.other - x.AG_uc.self # positive values mean t4 has x amount greatr flux than t1 etc
    x['GD_diff'] = x.GD_uc.other - x.GD_uc.self # positive values mean t4 has x amount greatr flux than t1 etc
    x['medium_diff'] = x.medium.other - x.medium.self # positive values mean t4 has x amount greater flux than t1 etc # not needed.
    
    # Grab names of internal exchange fluxes and trim to make compatible with flux output. 
    int_ex = [i.id for i in com.internal_exchanges]
    int_ex = [i[:-7] for i in int_ex]
    # int_ex

    # Narrow focus to exchange fluxes
    x = x[x.index.isin(int_ex)].rename(columns={'self': 'testA', 'other': 'testB'}, level=-1).sort_values(by=['GD_diff'], ascending = True)
    # x.drop(["medium","medium_diff"])
    del x["medium"]
    del x["medium_diff"]
    
    return x



def main():
    pass
    # aaSensitivity()

if __name__ == "__main__":
    main()

