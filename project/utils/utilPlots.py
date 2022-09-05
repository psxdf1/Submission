
def generate_uptake_sensitiviy_data(model, exRxn ,uptakeMaximum, step):
    '''
    Intended to calculate the effect of variance of a single aa availiblity on the growth rate

    Returns: 
        Should return the growth rate corresponding to uptake constraints.
    '''
    with model as m:
        x = uptakeMaximum
        bound = []
        exRxn_Sensitivity = [] # List to contain growth rate in response to aa availibility

        rxn = [i for i in model.reactions if i.id == exRxn][0]

        while x <= 0:
            rxn.lower_bound = x
            exRxn_Sensitivity.append(m.slim_optimize())
            bound.append(abs(x))
            x += step

    # Just reverse both lists for ridyness
    bound.reverse()
    exRxn_Sensitivity.reverse()

    return bound, exRxn_Sensitivity


def plotIntExchFlux(com, exchngFluxes, _title = "Internal Exchange Flux Values - Scenario X",save_output=False):
   
    import matplotlib.pyplot as plt
    import numpy as np

    # int_ex = [i.id for i in com.internal_exchanges]
    # int_ex = [i[:-7] for i in int_ex]

    # x = sol.fluxes.T

    # x = x[x.index.isin(int_ex)].rename(columns={'self': 'testA', 'other': 'testB'}, level=-1).sort_values(by=['GD_uc'], ascending = False)
    
    # Start plotting
    fig, ax = plt.subplots(figsize=(20,10))
    x_vals = np.arange(len(exchngFluxes.index))

    y1 = [i for i in exchngFluxes.AG_uc]
    y2 = [i for i in exchngFluxes.GD_uc]

    width = 0.4
    offset = width//2

    ax.bar(x_vals + 0.1, y1, width=width, label="AG")
    ax.bar(x_vals - 0.1, y2, width=width, label="GD")

    ax.tick_params("x", rotation=90, labelsize=15)
    ax.set_xticks(x_vals,[i for i in exchngFluxes.index])
    ax.tick_params(axis='y', labelsize=20)

    # print(ax.get_ylim())
    # ax.set_yscale('log')
    ax.set_yscale('symlog')

    ax.set_xlabel('Internal Exhange Reaction',fontsize=20)
    ax.set_ylabel('Flux Value Difference',fontsize=20)

    ax.legend()
    ax.set_title(f"{_title}",fontsize=30)
    # ax.title.set_size(40)

    ax.grid()

    if save_output:
        fig.savefig(f"C:/Users/domin/Masters-Thesis---Draft-1/images/Micom_Tests/{_title}.pdf", bbox_inches = 'tight',pad_inches=0.2)
    






def plot_IntExchDiff(x,_title="Differences Between Internal Exchange Reactions; Scenario x vs Scenario x",save_output=False):
    from matplotlib import pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(figsize=(20,10))
    x_vals = np.arange(len(x.index))

    y1 = [i for i in x.AG_diff]
    y2 = [i for i in x.GD_diff]

    width = 0.4
    offset = width//2

    ax.bar(x_vals + 0.1, y1, width=width, label="AG")
    ax.bar(x_vals - 0.1, y2, width=width, label="GD")

    ax.tick_params("x", rotation=90)
    ax.set_xticks(x_vals,[i for i in x.index])
    ax.tick_params(axis='y', labelsize=20)

    ax.set_yscale('symlog')

    ax.set_xlabel('Internal Exhange Reaction',fontsize=20)
    ax.set_ylabel('Flux Value Difference',fontsize=20)

    ax.legend()
    ax.set_title(f"{_title}",fontsize=30)
    # ax.title.set_size(40)

    ax.grid()
    if save_output:
        fig.savefig(f"C:/Users/domin/Masters-Thesis---Draft-1/images/Micom_Tests/{_title}.pdf", bbox_inches = 'tight',pad_inches=0.2)
    