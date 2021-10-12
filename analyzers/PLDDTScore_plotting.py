import sys
sys.path.append("../alphafold_analysis_for_mutation")

import matplotlib.pyplot as plt
import pandas as pd

configs = [{
    "statistics_name": "min",
    "ylabel": "Minimum",
    "mutation_site_vs_plddt": 'mutation_site_vs_minimum_plddt_confidence',
    "wt_versus_plddt": "minimum_plddt_confidence_for_each_wt"
    },
    {
    "statistics_name": "max",
    "ylabel": "Maximum",
    "mutation_site_vs_plddt": 'mutation_site_vs_maximum_plddt_confidence',
    "wt_versus_plddt": "maximum_plddt_confidence_for_each_wt"
    },
    {
    "statistics_name": "avg",
    "ylabel": "Average",
    "mutation_site_vs_plddt": 'mutation_site_vs_average_plddt_confidence',
    "wt_versus_plddt": "average_plddt_confidence_for_each_wt"
    },
    {
    "statistics_name": "median",
    "ylabel": "Median",
    "mutation_site_vs_plddt": 'mutation_site_vs_median_plddt_confidence',
    "wt_versus_plddt": "median_plddt_confidence_for_each_wt"
    },
    {
    "statistics_name": "std",
    "ylabel": "Standard deviation",
    "mutation_site_vs_plddt": 'mutation_site_vs_std_plddt_confidence',
    "wt_versus_plddt": "std_plddt_confidence_for_each_wt"
    },
]


def get_mutation_site_vs_statistics(wt_plddt_statistics_df, 
                                    mt_plddt_statistics_df,
                                    statistics_name="median", 
                                    ylabel="Median", 
                                    filename="mutation_site_vs_median_plddt_confidence",
                                    do_plot=False,
                                    neighbor=0):
    """[summary]

    Args:
        wt_plddt_statistics_df ([type]): [description]
        mt_plddt_statistics_df ([type]): [description]
        statistics_name (str, optional): [description]. Defaults to "median".
        ylabel (str, optional): [description]. Defaults to "Median".
        filename (str, optional): [description]. Defaults to "mutation_site_vs_median_plddt_confidence".
        do_plot (bool, optional): [description]. Defaults to False.
        neighbor (int, optional): [description]. Defaults to 0.

    Returns:
        [type]: [description]
    """
    # wt_plddt_statistics_df = pd.read_excel("outputs/wt_mutation_site_plddt_statistics.xlsx")
    # mt_plddt_statistics_df = pd.read_excel("outputs/mt_mutation_site_plddt_statistics.xlsx")
    x_wt_mutation_site = wt_plddt_statistics_df["mutation_site"].astype(int)
    y_wt_plddt_stat = wt_plddt_statistics_df[statistics_name].astype(float)
    x_mt_mutation_site = mt_plddt_statistics_df["mutation_site"].astype(int)
    y_mt_plddt_stat = mt_plddt_statistics_df[statistics_name].astype(float)
    
    if do_plot:
        fig = plt.figure()
        gs = fig.add_gridspec(2, hspace=0)
        axs = gs.subplots(sharex=True, sharey=True)
        axs[0].scatter(x_wt_mutation_site, y_wt_plddt_stat, color='green', label="Wild-type")
        axs[0].legend()
        axs[1].scatter(x_mt_mutation_site, y_mt_plddt_stat, color='red', label="Variant")
        axs[1].legend()
        for ax in axs:
            ax.label_outer()
        plt.legend(loc="best")
        plt.xlabel("Mutation site")
        plt.ylabel("Confidence PLDDT ({})".format(ylabel.lower()))
        # plt.show()
        plt.savefig("output_images/mutation_site_specific_confidence_plots/{}_{}_neighbor.pdf".format(filename, neighbor), dpi=300, format="pdf")
        plt.close()
    
    return x_wt_mutation_site, y_wt_plddt_stat, x_mt_mutation_site, y_mt_plddt_stat


def plot_mutation_site_vs_statistics_1():
    """For five plots with alignment.
    """
    fig = plt.figure(figsize=(10, 6))
    fig.text(0.5, 0.01, "Mutation site", ha='center')
    fig.text(0.02, 0.5, "Confidence (PLDDT)", va='center', rotation='vertical')
    gs = fig.add_gridspec(5,6, hspace=0.0)
    # gs.subplots(sharex=True, sharey=False)
    x_wt_mutation_site, y_wt_plddt_stat, x_mt_mutation_site, y_mt_plddt_stat = get_mutation_site_vs_statistics(statistics_name=configs[0]["statistics_name"], ylabel=configs[0]["ylabel"], filename=configs[0]["mutation_site_vs_plddt"], do_plot=False)
    ax1 = fig.add_subplot(gs[0, 0:2])
    ax2 = fig.add_subplot(gs[1, 0:2])
    ax1.scatter(x_wt_mutation_site, y_wt_plddt_stat, color='green', label="Wild-type")
    ax2.scatter(x_mt_mutation_site, y_mt_plddt_stat, color='red', label="Variant")
    ax2.set_ylabel(configs[0]["ylabel"])
    
    x_wt_mutation_site, y_wt_plddt_stat, x_mt_mutation_site, y_mt_plddt_stat = get_mutation_site_vs_statistics(statistics_name=configs[1]["statistics_name"], ylabel=configs[1]["ylabel"], filename=configs[1]["mutation_site_vs_plddt"], do_plot=False)
    ax3 = fig.add_subplot(gs[0, 2:4])
    ax4 = fig.add_subplot(gs[1, 2:4])
    ax3.scatter(x_wt_mutation_site, y_wt_plddt_stat, color='green', label="Wild-type")
    ax4.scatter(x_mt_mutation_site, y_mt_plddt_stat, color='red', label="Variant")
    ax4.set_ylabel(configs[1]["ylabel"])
    
    x_wt_mutation_site, y_wt_plddt_stat, x_mt_mutation_site, y_mt_plddt_stat = get_mutation_site_vs_statistics(statistics_name=configs[2]["statistics_name"], ylabel=configs[2]["ylabel"], filename=configs[2]["mutation_site_vs_plddt"], do_plot=False)
    ax5 = fig.add_subplot(gs[0, 4:6])
    ax6 = fig.add_subplot(gs[1, 4:6])
    ax5.scatter(x_wt_mutation_site, y_wt_plddt_stat, color='green', label="Wild-type")
    ax6.scatter(x_mt_mutation_site, y_mt_plddt_stat, color='red', label="Variant")
    ax6.set_ylabel(configs[2]["ylabel"])
    
    x_wt_mutation_site, y_wt_plddt_stat, x_mt_mutation_site, y_mt_plddt_stat = get_mutation_site_vs_statistics(statistics_name=configs[3]["statistics_name"], ylabel=configs[3]["ylabel"], filename=configs[3]["mutation_site_vs_plddt"], do_plot=False)
    ax7 = fig.add_subplot(gs[3, 1:3])
    ax8 = fig.add_subplot(gs[4, 1:3])
    ax7.scatter(x_wt_mutation_site, y_wt_plddt_stat, color='green', label="Wild-type")
    ax8.scatter(x_mt_mutation_site, y_mt_plddt_stat, color='red', label="Variant")
    ax8.set_ylabel(configs[3]["ylabel"])
    
    x_wt_mutation_site, y_wt_plddt_stat, x_mt_mutation_site, y_mt_plddt_stat = get_mutation_site_vs_statistics(statistics_name=configs[4]["statistics_name"], ylabel=configs[4]["ylabel"], filename=configs[4]["mutation_site_vs_plddt"], do_plot=False)
    ax9 = fig.add_subplot(gs[3, 3:5])
    ax10 = fig.add_subplot(gs[4, 3:5])
    ax9.scatter(x_wt_mutation_site, y_wt_plddt_stat, color='green', label="Wild-type")
    ax10.scatter(x_mt_mutation_site, y_mt_plddt_stat, color='red', label="Variant")
    ax10.set_ylabel(configs[4]["ylabel"])
    
    fig.tight_layout()
    # plt.show()
    plt.savefig("output_images/plddt_confidence_at_mutation_site.pdf", dpi=300, format="pdf")
    plt.close()

    
def plot_mutation_site_vs_statistics_0():
    """No alighment of the figures. But will work.
    """
    fig = plt.figure()
    fig.text(0.5, 0.01, "Mutation site", ha='center')
    fig.text(0.02, 0.5, "Confidence (PLDDT)", va='center', rotation='vertical')
    gs = fig.add_gridspec(5,3, hspace=0.0)
    axs = gs.subplots(sharex=True, sharey=False)
    x_wt_mutation_site, y_wt_plddt_stat, x_mt_mutation_site, y_mt_plddt_stat = get_mutation_site_vs_statistics(statistics_name=configs[0]["statistics_name"], ylabel=configs[0]["ylabel"], filename=configs[0]["mutation_site_vs_plddt"], do_plot=False)
    axs[0,0].scatter(x_wt_mutation_site, y_wt_plddt_stat, color='green', label="Wild-type")
    axs[1,0].scatter(x_mt_mutation_site, y_mt_plddt_stat, color='red', label="Variant")
    axs[1,0].set_ylabel(configs[0]["ylabel"])
    
    x_wt_mutation_site, y_wt_plddt_stat, x_mt_mutation_site, y_mt_plddt_stat = get_mutation_site_vs_statistics(statistics_name=configs[1]["statistics_name"], ylabel=configs[1]["ylabel"], filename=configs[1]["mutation_site_vs_plddt"], do_plot=False)
    axs[0,1].scatter(x_wt_mutation_site, y_wt_plddt_stat, color='green', label="Wild-type")
    axs[1,1].scatter(x_mt_mutation_site, y_mt_plddt_stat, color='red', label="Variant")
    axs[1,1].set_ylabel(configs[1]["ylabel"])
    
    x_wt_mutation_site, y_wt_plddt_stat, x_mt_mutation_site, y_mt_plddt_stat = get_mutation_site_vs_statistics(statistics_name=configs[2]["statistics_name"], ylabel=configs[2]["ylabel"], filename=configs[2]["mutation_site_vs_plddt"], do_plot=False)
    axs[0,2].scatter(x_wt_mutation_site, y_wt_plddt_stat, color='green', label="Wild-type")
    axs[1,2].scatter(x_mt_mutation_site, y_mt_plddt_stat, color='red', label="Variant")
    axs[1,2].set_ylabel(configs[2]["ylabel"])
    
    axs[2,0].remove()
    axs[2,1].remove()
    axs[2,2].remove()
    
    x_wt_mutation_site, y_wt_plddt_stat, x_mt_mutation_site, y_mt_plddt_stat = get_mutation_site_vs_statistics(statistics_name=configs[3]["statistics_name"], ylabel=configs[3]["ylabel"], filename=configs[3]["mutation_site_vs_plddt"], do_plot=False)
    axs[3,0].scatter(x_wt_mutation_site, y_wt_plddt_stat, color='green', label="Wild-type")
    axs[4,0].scatter(x_mt_mutation_site, y_mt_plddt_stat, color='red', label="Variant")
    axs[4,0].set_ylabel(configs[3]["ylabel"])
    
    x_wt_mutation_site, y_wt_plddt_stat, x_mt_mutation_site, y_mt_plddt_stat = get_mutation_site_vs_statistics(statistics_name=configs[4]["statistics_name"], ylabel=configs[4]["ylabel"], filename=configs[4]["mutation_site_vs_plddt"], do_plot=False)
    axs[3,1].scatter(x_wt_mutation_site, y_wt_plddt_stat, color='green', label="Wild-type")
    axs[4,1].scatter(x_mt_mutation_site, y_mt_plddt_stat, color='red', label="Variant")
    axs[4,1].set_ylabel(configs[4]["ylabel"])
    
    axs[3,2].remove()
    axs[4,2].remove()
    
    fig.tight_layout(pad=1)
    plt.show()
# plot_mutation_site_vs_statistics_0() 

def get_plddt_confidence_score_for_all_WT_proteins(mt_plddt_statistics_df, 
                                                   statistics_name="median", 
                                                   ylabel="Median", 
                                                   filename="minimum_plddt_confidence_for_each_wt",
                                                   do_plot=False,
                                                   neighbor=0):
    """It returns the x, y values for plotting confidence vs. WT 
    as well as single-plots if do_plot is true

    Args:
        mt_plddt_statistics_df ([type]): [description]
        statistics_name (str, optional): [description]. Defaults to "median".
        ylabel (str, optional): [description]. Defaults to "Median".
        filename (str, optional): [description]. Defaults to "minimum_plddt_confidence_for_each_wt".
        do_plot (bool, optional): [description]. Defaults to False.
        neighbor (int, optional): [description]. Defaults to 0.

    Returns:
        [type]: [description]
    """
    ssym_classified_full_df = pd.read_excel("data/ssym_classified_full.xlsx")
    # mt_plddt_statistics_df = pd.read_excel("outputs/mt_mutation_site_plddt_statistics.xlsx")
    mt_plddt_statistics_df = mt_plddt_statistics_df.rename(columns={"pdb_id": "inverse_pdb_id"})
    mt_plddt_statistics_with_wtpdbid_df = mt_plddt_statistics_df.merge(right=ssym_classified_full_df[["pdb_id", "inverse_pdb_id"]],
                                                                   how="left", left_on="inverse_pdb_id", right_on="inverse_pdb_id", 
                                                                   left_index=False, right_index=False,)  
    wt_specific_plddt =[]
    wt_pdb_ids = mt_plddt_statistics_with_wtpdbid_df["pdb_id"].unique()
    for wt_pdb_id in wt_pdb_ids:
        x = pd.DataFrame()
        x[wt_pdb_id] = mt_plddt_statistics_with_wtpdbid_df[mt_plddt_statistics_with_wtpdbid_df["pdb_id"]==wt_pdb_id][statistics_name]
        wt_specific_plddt.append(x.reset_index(drop=True))

    temp_df = pd.DataFrame()
    temp_df = pd.concat(wt_specific_plddt, axis=1)
    temp_df.columns = [''] * len(temp_df.columns) # removed the column names as they were pdb id
    # print(temp_df.head())
    if do_plot:
        temp_df.boxplot(grid=True) 
        plt.legend(loc="best")
        plt.xlabel("Dataset (Wild-type)")
        plt.ylabel("Confidence PLDDT ({})".format(ylabel.lower()))
        # plt.show()
        plt.savefig("output_images/wt_specific_plddt_confidence_plots/{}_{}_neighbor.pdf".format(filename, neighbor), dpi=300, format="pdf")
        plt.close()
    return temp_df
  
    

def plot_plddt_confidence_score_for_all_WT_proteins():
    fig, axs = plt.subplots(nrows=2, ncols=3, sharex=True, figsize=(10, 6))
    fig.text(0.5, 0.01, "Dataset (Wild-type)", ha='center')
    fig.text(0.02, 0.5, "Confidence (PLDDT)", va='center', rotation='vertical')
    k=-1
    for i, x in enumerate(configs):
        temp_df = get_plddt_confidence_score_for_all_WT_proteins(statistics_name=x["statistics_name"], 
                                                ylabel=x["ylabel"], filename=x["wt_versus_plddt"])
        if i%3 == 0: k=k+1
        # print(k, i%3)
        temp_df.boxplot(grid=True, ax=axs[k, i%3])
        axs[k, i%3].set_ylabel(x["ylabel"])

    axs[1, 2].remove()
    fig.tight_layout()
    plt.show()


def plot_plddt_confidence_score_for_all_WT_proteins_1():
    """Better image with colspan.
    """
    fig = plt.figure(0)
    fig.text(0.5, 0.01, "Dataset (Wild-type)", ha='center')
    fig.text(0.02, 0.5, "Confidence (PLDDT)", va='center', rotation='vertical')
    axs = []
    axs.append(plt.subplot2grid((2,6), (0,0), colspan=2))
    axs.append(plt.subplot2grid((2,6), (0,2), colspan=2))
    axs.append(plt.subplot2grid((2,6), (0,4), colspan=2))
    axs.append(plt.subplot2grid((2,6), (1,1), colspan=2))
    axs.append(plt.subplot2grid((2,6), (1,3), colspan=2))

    for i, x in enumerate(configs):
        temp_df = get_plddt_confidence_score_for_all_WT_proteins(statistics_name=x["statistics_name"], 
                                                                 ylabel=x["ylabel"], 
                                                                 filename=x["wt_versus_plddt"],
                                                                 do_plot=False)
        plt.cla()
        temp_df.boxplot(grid=True, ax=axs[i])
        axs[i].set_ylabel(x["ylabel"])
        

    fig.tight_layout(pad=1)
    # plt.show()
    plt.savefig("output_images/plddt_confidence_vs_dataset.pdf", dpi=300, format="pdf")
    plt.close()

# plot_plddt_confidence_score_for_all_WT_proteins_1()
# plot_mutation_site_vs_statistics_1()


for config in configs:
    for neighbor in range(3):
        wt_plddt_statistics_df = pd.read_excel("outputs/wt_mutation_plddt_statistics_{}_neighbor.xlsx".format(neighbor))
        mt_plddt_statistics_df = pd.read_excel("outputs/mt_mutation_plddt_statistics_{}_neighbor.xlsx".format(neighbor))
        statistics_name = config["statistics_name"]
        ylabel = config["ylabel"]
        filename = config["mutation_site_vs_plddt"]
        do_plot = True
        get_mutation_site_vs_statistics(wt_plddt_statistics_df, mt_plddt_statistics_df, 
                                        statistics_name, ylabel, filename, do_plot, neighbor)
    # break
    
for config in configs:
    for neighbor in range(3):
        mt_plddt_statistics_df = pd.read_excel("outputs/mt_mutation_plddt_statistics_{}_neighbor.xlsx".format(neighbor))
        statistics_name = config["statistics_name"]
        ylabel = config["ylabel"]
        filename = config["mutation_site_vs_plddt"]
        do_plot = True
        get_plddt_confidence_score_for_all_WT_proteins(mt_plddt_statistics_df, 
                                        statistics_name, ylabel, filename, do_plot, neighbor)