from utils import *

def load_model(model_name):
    data = []
    with open(f'target/dms/dms_{model_name}.log') as f:
        for line in f:
            line = line.rstrip().split(' | ')[-1]

            if line.startswith('Results for '):
                curr_prot = line.split('/')[-1].split('.')[0]
                continue

            if line.startswith('\tDMS'):
                line = line.rstrip(':').strip()
                # Split on the first dash after "DMS_"
                dms_part, model_part = line.split('-', 1)
                dms_name = dms_part.strip()
                model_name = model_part.strip()
                continue

            if line.startswith('\t\tSpearman r = '):
                corr = abs(float(line.split(',')[0].split()[-1]))
                p = float(line.split()[-1])
                data.append([curr_prot, dms_name, model_name, corr, p])

    return data

"""
def calculate_win_rates(df):
    # Uses iindividual not avearge 
    # Group by protein and find the model with the highest correlation for each
    winners = df.loc[df.groupby('protein')['corr'].idxmax()]
    total_proteins = len(winners)

    # Calculate overall win rates
    overall_win_rates = winners['model_name'].value_counts() / total_proteins

    # Calculate head-to-head win rate for ESM2 vs ESM1b
    head_to_head = df[df['model_name'].isin(['esm2', 'esm1b'])].copy()
    head_to_head['rank'] = head_to_head.groupby('protein')['corr'].rank(ascending=False, method='min')
    esm2_wins = (head_to_head[head_to_head['model_name'] == 'esm2']['rank'] == 1).sum()
    total_comparisons = len(head_to_head['protein'].unique())
    esm2_win_rate = esm2_wins / total_comparisons

    return overall_win_rates, esm2_win_rate
"""

def calculate_win_rates(df):
    # Calculate average correlation for each protein and model
    avg_corr = df.groupby(['protein', 'model_name'])['corr'].mean().reset_index()
    
    # Find the winner for each protein
    winners = avg_corr.loc[avg_corr.groupby('protein')['corr'].idxmax()]
    total_proteins = len(winners['protein'].unique())

    # Calculate overall win rates
    overall_win_rates = winners['model_name'].value_counts() / total_proteins

    # Calculate head-to-head win rate for ESM2 vs ESM1b
    head_to_head = avg_corr[avg_corr['model_name'].isin(['esm2', 'esm1b'])].copy()
    head_to_head['rank'] = head_to_head.groupby('protein')['corr'].rank(ascending=False, method='min')
    esm2_wins = (head_to_head[head_to_head['model_name'] == 'esm2']['rank'] == 1).sum()
    total_comparisons = len(head_to_head['protein'].unique())
    esm2_win_rate = esm2_wins / total_comparisons

    return overall_win_rates, esm2_win_rate
        
if __name__ == '__main__':
    data = []
    data += load_model('esm2-3B')  # New model
    data += load_model('esm2')
    data += load_model('esm1b')
    data += load_model('esm1-670D')  # New model
    data += load_model('tape')
    

    df = pd.DataFrame(data, columns=[
        'protein',
        'dms',
        'model_name',
        'corr',
        'pval',
    ])

    # Remove duplicates based on 'protein', 'dms', and 'model_name'
    df_deduped = df.drop_duplicates(subset=['protein', 'dms', 'model_name']).reset_index()

    # Print information about removed duplicates
    num_duplicates = len(df) - len(df_deduped)
    print(f"Removed {num_duplicates} duplicate entries")

    # Use the deduplicated DataFrame for further processing
    df = df_deduped

    # Calculate win rates
    overall_win_rates, esm2_vs_esm1b_win_rate = calculate_win_rates(df)

    # Print win rates
    print("Overall Win Rates:")
    print(overall_win_rates)
    print(f"\nESM2 vs ESM1b Win Rate: {esm2_vs_esm1b_win_rate:.2f}")

    # Define the color scheme and order
    color_dict = {
        'esm2-3B': '#ff0000',    # Red
        'esm2': '#2ca02c',        # Green
        'esm1b': '#1f77b4',       # Dark Blue
        'DeepSequence': '#ff7f0e',  # Orange
        'tape': '#9ecae1',        # Light Blue
        'esm1-670D': '#808080'    # Gray
    }
    model_order = ['esm2-3B', 'esm2', 'esm1b', 'DeepSequence', 'esm1-670D', 'tape']

    # Set the style for seaborn
    sns.set_style("whitegrid")

    # Ensure the DataFrame has 'model_name' as a categorical type with the desired order
    df['model_name'] = pd.Categorical(df['model_name'], categories=model_order, ordered=True)

    plt.figure(figsize=(14, 5))  # Increased figure width to accommodate new models

    # Create the bar plot
    ax = sns.barplot(
        data=df,
        x='protein',
        y='corr',
        hue='model_name',
        hue_order=model_order,
        palette=color_dict,
        ci=None,
    )

    # Add the strip plot
    sns.stripplot(
        data=df,
        x='protein',
        y='corr',
        hue='model_name',
        hue_order=model_order,
        palette=color_dict,
        dodge=True,
        jitter=True,
        size=4,
        edgecolor='black',
        linewidth=0.5,
    )

    # Customize the plot
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('|Spearman r|')
    plt.xlabel('Protein')
    plt.title('DMS Plot', fontsize=16)

    # Remove top and right spines
    sns.despine()

    # Adjust legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:6], labels[:6], title='Model Name', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Adjust layout and save
    plt.tight_layout()
    plt.savefig('figures/plot_dms.svg', bbox_inches='tight')
    plt.close()

    print(df.to_csv())