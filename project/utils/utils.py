
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
    # x = m.slim_optimize()
    # if x > 0: # this is wonky
    #     return m
    # else:
    #     print(x)
    #     print("Not able to optimize after bound modifications")
    #     return False

    return model


def runcases(model,path_to_cases, tgt):
    # Aimed at testing different internal flux bound cases in individual models
    i = 0
    m_gr = []
    m_tmy = []

    for i in range(3):
        fn = str(i) + ".csv"
        with model as m: 
            m = set_bound_modifications(m,f"{path_to_cases}/{fn}")
            m_gr.append(m.slim_optimize())
            m.objective = tgt
            m_tmy.append(m.slim_optimize())

    return m_gr, m_tmy




def set_media(model, file, display= False):
    '''
    Gets a media representation from a file and ..
    '''

    pass


def indMod_case(model, modifications):
    gr = 0
    tmy = 0

    return gr, tmy # growth rate, theoretical maximum yeiled


def main():
    pass

if __name__ == "__main__":
    mod = [1]
    set_bound_modifications(mod,"project/data/ind_mod_cases/AG_case1.csv")

