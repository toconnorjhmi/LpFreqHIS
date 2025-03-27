import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# set up seaborn parameters
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

def h_histogram(h_freq_df, file, bins=50, label_offset = 1.5):
    h_freq_df['Freq'] = h_freq_df['Freq'] * 100
    lpg1596 = 2.38 # percent his in 1596

    fig, ax = plt.subplots(figsize=(8, 5))
    sns.histplot(h_freq_df, bins = bins, x = "Freq")

    plt.xlabel("Percent Histidines")
    plt.ylabel("Frequency")
    
    # Add vertical line at lpg1596
    ax.axvline(lpg1596, color='black', linestyle='dashed', linewidth=1, label=f'x = {lpg1596}')
    ax.text(lpg1596 + label_offset, 
            ax.get_ylim()[1]*0.9, 'Lpg1596', color='black', ha='center', 
            fontsize=12, fontstyle='italic',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

    ax.set_xlim(0, 18)
    ax.set_xticks(range(0, 19, 2))
    plt.savefig(file)
    plt.show()

def overlaid_plot(k12_df, lp02_df, file, bins=50):
    k12_df['Freq'] = k12_df['Freq'] * 100
    lp02_df['Freq'] = lp02_df['Freq'] * 100

    fig, ax = plt.subplots(figsize=(8, 5))
    sns.histplot(k12_df, bins = bins, x = "Freq", alpha=0.4, kde=True, label='$\it{E. Coli} \ K12$')
    sns.histplot(lp02_df, bins=bins, x="Freq", alpha=0.4, kde=True, label='$\it{L. pneumophila}$')

    plt.xlabel("Percent Histidines")
    plt.legend()
    plt.ylabel("Frequency")

    ax.set_xlim(0, 18)
    ax.set_xticks(range(0, 19, 2))
    ax.set_yticks(range(0, 401, 50))
    plt.savefig(file)
    plt.show()    

def lpg1596_sliding_window(file):
    h_freq_df = pd.read_csv("../out/sliding_window_1596.csv")
    
    fig, ax = plt.subplots(figsize=(8,5))
    ax.set_xlim(0, 672)
    ax.set_ylim(0, 2)
    ax.yaxis.get_major_locator().set_params(integer=True)
    sns.lineplot(h_freq_df, x='Position', y='Histidine_Count')
    

    plt.xlabel("Position")
    plt.ylabel("# Histidines")
    plt.savefig(file)
    plt.show()

def ORF_sliding_window(window_size = 24):
    df = pd.read_csv('H_run_Lp02ORFs.tsv', sep='\t')  
    window_size = 12
    rolling_max = df['longest_H_run'].rolling(window=window_size).mean()

    result_df = pd.DataFrame({
        'original_index': range(len(df)),
        'max_H_run_length': rolling_max
    })

    # Remove NaN values from the beginning (first window_size-1 rows)
    result_df = result_df.dropna()

    # include the gene name from the original data
    result_df['gene_name'] = df['name'].iloc[result_df['original_index'].values].values

    return result_df

def orf_window_plot(df, window_size, file):
    fig, ax = plt.subplots(figsize=(8,5))
    
    sns.lineplot(df, x='original_index', y='max_H_run_length')
    plt.ylabel(f"Mean Run Of Histidines/{window_size} genes")
    plt.xlabel(None)
    ax.set_xlim(0, 3087)
    ax.set_ylim(0.7, 1.71)

    plt.savefig(file)
    plt.show()


if __name__ == "__main__":
    h_freq_df_lp = pd.read_csv("./out/NCBI_Lp_histFreq.csv")
    h_freq_k12 = pd.read_csv("./out/NCBI_K12_histFreq.csv")

    h_histogram(h_freq_k12, file ='./plots/k12_h_histogram.pdf')
    h_histogram(h_freq_df_lp, file = "./plots/lp_h_histogram.pdf")
    overlaid_plot(h_freq_k12, h_freq_df_lp, file = "./plots/k12_lp02_histoverlay_marker.pdf")

    lpg1596_sliding_window(file='./plots/1596_sliding_window.pdf')
    orf_window_df = ORF_sliding_window()
    orf_window_plot(orf_window_df, window_size=24, file = 'ORF_H_rollingmean.pdf')

    