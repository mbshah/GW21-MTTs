import math
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from scipy.stats import percentileofscore
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import ast
import joblib

# Data preparation part
working_dir = "/home/hb0358/Disk2/Disk2/Data/MTT_mbs/Statistics_GW21/"
pathway_ko_map_file = working_dir + "pathway_pedia_only_protist.tsv"
pathway_mappingdf = pd.read_csv(pathway_ko_map_file, sep="\t", header=0)
pathway_mappingdf.set_index(pathway_mappingdf.columns[0], inplace=True)
#pathway_mappingdf.drop(pathway_mappingdf.columns[0], axis=1, inplace=True)
pathway_mappingdf['KEGG_Orthologs'] = pathway_mappingdf['KEGG_Orthologs'].str.split('; ').apply(tuple)

# params
n_clusters = 6
time_pointsl = [1, 12, 24, 48, 96, 168, 240]
col_name_to_cluster='Log2FC_over_time'
cluster_level="ko" #change to litrally anything else to use pathways keep "ko" for orthologs
filter_data=True ## currently only works with pathways not KOs

# output
final_folder=f"ko_{n_clusters}_protist" if cluster_level=="ko" else f"pathway_{n_clusters}_protist"
if filter_data:final_folder=final_folder+"_filtered"
cluster_outdir=working_dir + f'graphics/pathway_clustering/{final_folder}/'
if not os.path.exists(cluster_outdir): os.makedirs(cluster_outdir)
output_file_pre = f"{cluster_outdir}/ortholog_cluster"if cluster_level=="ko" else f"{cluster_outdir}/pathway_cluster"
clusterFile1=f"{cluster_outdir}/cluster_info.tsv"
cluster_pickle=f"{cluster_outdir}/cluster_info.pkl"

#declare input_files
if cluster_level=="ko":
    p_value_file=f'{working_dir}KEGG_DESEQ_P_vals_PROTISTS.tsv'
else:
    p_value_file = f'{working_dir}KEGG_Pathways_DESEQ_P_vals_PROTISTS.tsv'

profile_to_use=f"{working_dir}KO_counts_HMM_version_trs_pred_0.001_PROTIST_Normalised.tsv"
ko_clustering_input=profile_to_use.replace(".tsv","_log2f_inputs.tsv")#providing this file will significantly reduce the time requred (has precalculated log2fold change values)

print(f"Current Configuration:")
print(f"Cluster Level: {cluster_level}")
print(f"Profile file being used: {profile_to_use}")
print(f"Filteration:{filter_data}")
if filter_data:print(f"P_Val file: {p_value_file}")
print(f"Number of Clusters: {n_clusters}")
if cluster_level=="ko": print(f"Does log_2fc_calculations_exist: {os.path.exists(ko_clustering_input)}")
print(f"Was this configuration clustering done before: {os.path.exists(cluster_pickle)}")


def filter_data_p_val(counts, p_val_file=p_value_file,col_name="treatment_padj"):
    p_val_t=pd.read_csv(p_val_file, sep="\t", header=0)
    #print(counts.columns)
    col_name_for_subject="KO" if cluster_level=="ko" else "Pathway"
    significant_pathways=p_val_t[p_val_t[col_name] < 0.05]["gene"].tolist()
    print(f"only {len(significant_pathways)} used for further processing\n{significant_pathways}")
    if cluster_level=="ko":
        filtered_counts=counts[counts.index.isin(significant_pathways)]
    else:
        filtered_counts=counts[counts[col_name_for_subject].isin(significant_pathways)]
    #filtered_counts = filtered_counts[~filtered_counts['Subsystems'].str.contains("Organismal Systems|Human Diseases", na=False)]
    return filtered_counts, significant_pathways


def create_input(): #if we want to cluster at pathway level we have to provide KO_profile and pathway profile in addition to the mapping file provided at the top
    KO_file = profile_to_use
    #pathway_file = working_dir + "KEGG_Pathway_profile_HMM_trs_pred_0.001_edgeR_normalised.tsv"
    my_ko_table_full = pd.read_csv(KO_file, sep="\t", header=0, index_col=0)
    my_pway_table_full = pathway_mappingdf #pd.read_csv(pathway_file, sep="\t", header=0, index_col=0)

    extracted_data = []

    for pway in pathway_mappingdf.index:
        if pway in my_pway_table_full.index:  # and pway == "map00073":
            interest_df = my_ko_table_full[my_ko_table_full.index.isin(pathway_mappingdf.loc[pway, "KEGG_Orthologs"])]
            time_points = sorted(set(sample.split("_")[1] for sample in interest_df.columns))

            # Create a new dataframe for sum of KOs
            sum_df = pd.DataFrame(index=["C1", "C2", "C3", "T1", "T2", "T3", "Avg_C", "Avg_T"], columns=time_points)

            # Calculate sum for each time point
            for time_point in time_points:
                control_samples = [f"C{i}_{time_point}" for i in range(1, 4)]
                test_samples = [f"T{i}_{time_point}" for i in range(1, 4)]

                # Sum of KOs for control samples
                sum_df.loc["C1", time_point] = interest_df[control_samples[0]].sum() + 1
                sum_df.loc["C2", time_point] = interest_df[control_samples[1]].sum() + 1
                sum_df.loc["C3", time_point] = interest_df[control_samples[2]].sum() + 1
                sum_df.loc["Avg_C", time_point] = sum_df.loc[["C1", "C2", "C3"], time_point].mean()

                # Sum of KOs for test samples
                sum_df.loc["T1", time_point] = interest_df[test_samples[0]].sum() + 1
                sum_df.loc["T2", time_point] = interest_df[test_samples[1]].sum() + 1
                sum_df.loc["T3", time_point] = interest_df[test_samples[2]].sum() + 1
                sum_df.loc["Avg_T", time_point] = sum_df.loc[["T1", "T2", "T3"], time_point].mean()
            #calculating log2FoldChange
            sum_df.loc["D"] = np.log2(sum_df.loc["Avg_T"].astype(float) / sum_df.loc["Avg_C"].astype(float))
            sum_df = sum_df.T
            to_append=(pway, sum_df['D'].values, sum_df['Avg_C'].values,sum_df['Avg_T'].values)
            extracted_data.append(to_append)
    data_table = pd.DataFrame([{
        "Pathway": item[0],
        "Log2FC_over_time": ";".join(str(item[1]).replace('[', '').replace(']', '').split()),
        #"Avg_C":";".join(str(item[2]).replace('[', '').replace(']', '').split()).split(";"),
        #"Avg_T":";".join(str(item[3]).replace('[', '').replace(']', '').split()).split(";")
        "Avg_C": item[2],
        "Avg_T": item[3]
    } for item in extracted_data])
    return_table = data_table.merge(pathway_mappingdf, left_on="Pathway", right_on=pathway_mappingdf.index, how="left")
    return_table.drop('KEGG_Orthologs', axis=1, inplace=True)
    return_table["index"] = return_table["Pathway"]
    return_table.set_index("index", inplace=True)
    return return_table


def create_ko_input():
    print("creating cluster input lines from KO table")
    KO_file = profile_to_use
    my_ko_table_full = pd.read_csv(KO_file, sep="\t", header=0, index_col=0)

    extracted_data = []

    for ko in my_ko_table_full.index:
        interest_df = my_ko_table_full.loc[[ko]]
        time_points = sorted(set(sample.split("_")[1] for sample in interest_df.columns))

        # Create a new dataframe for sum of KOs
        sum_df = pd.DataFrame(index=["C1", "C2", "C3", "T1", "T2", "T3", "Avg_C", "Avg_T"], columns=time_points)

        # Calculate sum for each time point
        for time_point in time_points:
            control_samples = [f"C{i}_{time_point}" for i in range(1, 4)]
            test_samples = [f"T{i}_{time_point}" for i in range(1, 4)]
            #print(test_samples)

            # Sum of KOs for control samples
            sum_df.loc["C1", time_point] = interest_df[control_samples[0]].sum() + 1
            sum_df.loc["C2", time_point] = interest_df[control_samples[1]].sum() + 1
            sum_df.loc["C3", time_point] = interest_df[control_samples[2]].sum() + 1
            sum_df.loc["Avg_C", time_point] = sum_df.loc[["C1", "C2", "C3"], time_point].mean()

            # Sum of KOs for test samples
            sum_df.loc["T1", time_point] = interest_df[test_samples[0]].sum() + 1
            sum_df.loc["T2", time_point] = interest_df[test_samples[1]].sum() + 1
            sum_df.loc["T3", time_point] = interest_df[test_samples[2]].sum() + 1
            sum_df.loc["Avg_T", time_point] = sum_df.loc[["T1", "T2", "T3"], time_point].mean()

        #calculating log2FoldChange
        sum_df.loc["D"] = np.log2(sum_df.loc["Avg_T"].astype(float) / sum_df.loc["Avg_C"].astype(float))
        sum_df = sum_df.T
        to_append = (ko, sum_df['D'].values, sum_df['Avg_C'].values, sum_df['Avg_T'].values)
        extracted_data.append(to_append)

    data_table = pd.DataFrame([{
        "KO": item[0],
        "Log2FC_over_time": ";".join(str(item[1]).replace('[', '').replace(']', '').split()),
        "Avg_C": item[2],
        "Avg_T": item[3]
    } for item in extracted_data])
    data_table.to_csv(ko_clustering_input, sep="\t", index=False)
    return data_table


def normalize_list(lst):
    min_val = min(lst)
    max_val = max(lst)
    #ret_val=[(x - min_val) / (max_val - min_val) if max_val != min_val else x for x in lst]
    ret_val=[math.log2(x) for x in lst]
    return ret_val


# Function to calculate z-scores
def calculate_z_scores_per_timepoint(values, means, stds):
    ret=[(x - mean) / std if std != 0 else 0 for x, mean, std in zip(values, means, stds)]
    return ret


def make_list(in_str):
    in_str=str(in_str)
    if in_str.__contains__(","):
        ret_lst=ast.literal_eval(in_str)
    else:
        ret_lst=ast.literal_eval(str(in_str).replace(" ",","))
    return ret_lst


def input_from_file(file):
    print(f"Reading {file} since it exists")
    return_table=pd.read_csv(file, sep="\t", header=0, index_col=0)

    return_table['Avg_C'] = return_table['Avg_C'].apply(lambda x: make_list(x))
    return_table['Avg_T'] = return_table['Avg_T'].apply(lambda x: make_list(x))
    return return_table


def cluster_and_pval(tablei,col_t_c=col_name_to_cluster,no_clusters=n_clusters):
    print("Clustering Now")
    # Prepare data for clustering
    extracted_data = [(index, [float(x) for x in row[col_t_c].split(';')]) for index, row in
                      tablei.iterrows()]
    data_matrix = np.array([item[1] for item in extracted_data])

    # Perform KMeans clustering
    kmeans = KMeans(n_clusters=no_clusters, random_state=42)
    clusters = kmeans.fit_predict(data_matrix)
    # Calculate the mean line for each cluster
    cluster_means = {}
    for cluster in range(no_clusters):
        cluster_indices = np.where(clusters == cluster)[0]
        cluster_data = data_matrix[cluster_indices]
        cluster_mean = np.mean(cluster_data, axis=0)
        cluster_means[cluster] = cluster_mean

    # Score each line based on its deviation from the cluster mean
    deviation_scores = []
    for i, (pway, diff) in enumerate(extracted_data):
        cluster = clusters[i]
        cluster_mean = cluster_means[cluster]
        deviation_score = np.sum((diff - cluster_mean) ** 2)  # Sum of squared deviations
        deviation_scores.append((pway, diff, cluster, deviation_score))

    # Perform permutation testing to calculate p-values
    num_permutations = 500  # Adjust this for pvalue-calcualtion
    permutation_scores = {i: [] for i in range(no_clusters)}

    for _ in range(num_permutations):
        permuted_clusters = np.random.permutation(clusters)
        for cluster in range(no_clusters):
            cluster_indices = np.where(permuted_clusters == cluster)[0]
            cluster_data = data_matrix[cluster_indices]
            if len(cluster_data) > 0:
                cluster_mean = np.mean(cluster_data, axis=0)
                for idx in cluster_indices:
                    diff = data_matrix[idx]
                    deviation_score = np.sum((diff - cluster_mean) ** 2)
                    permutation_scores[cluster].append(deviation_score)

    p_values = []
    for i, (pway, diff, cluster, score) in enumerate(deviation_scores):
        permutation_score_distribution = permutation_scores[cluster]
        p_value = (np.sum(np.array(permutation_score_distribution) <= score) + 1) / (num_permutations + 1)
        p_values.append(p_value)

    # Sort the lines by their deviation score (ascending order)
    sorted_deviation_scores = sorted(zip(deviation_scores, p_values), key=lambda x: x[0][3])

    # Convert np.array to DataFrame
    first_col="KO" if cluster_level=="ko" else "Pathway"
    cluster_data = pd.DataFrame([{
        first_col: item[0][0],
        'Cluster': item[0][2],
        'Deviation_Score': item[0][3],
        'cluster_membership_pvalue': item[1]
    } for item in sorted_deviation_scores])
    ret_table=tablei.merge(cluster_data, on=first_col, how='left')
    ret_table.to_csv(clusterFile1, sep="\t",index=False)
    joblib.dump(kmeans, cluster_pickle)
    return ret_table, kmeans


if cluster_level =="ko":
    input_table =  input_from_file(ko_clustering_input) if os.path.exists(ko_clustering_input) else create_ko_input()
else:
    input_table = create_input()
print(f"input table after loading/creating: {len(input_table.index)}")
if filter_data:
    input_table2, sig_kos=filter_data_p_val(input_table)
    input_table3, sig_kos3 = filter_data_p_val(input_table, col_name="interaction_padj")
    #print(f"input table after filtering: {len(input_table2.index)}")



if os.path.exists(clusterFile1):
    merged_table = input_from_file(clusterFile1)
    kmeans = joblib.load(cluster_pickle)
else:
    merged_table, kmeans = cluster_and_pval(input_table, col_name_to_cluster, n_clusters)

# Create a DataFrame to store cluster results
results_df = merged_table
results_df[col_name_to_cluster] = results_df[col_name_to_cluster].apply(lambda x: np.array(x.split(';')).astype(float))
filtered_results_df = results_df#[results_df['interaction_padj'] < 0.05]
results_df["Sig_cluster"] = ""
results_df["TS_cluster"] = ""
results_df["IS_cluster"] = ""
#Adding Z-score lines to the input table
avg_c_df = pd.DataFrame(results_df['Avg_C'].tolist(), columns=time_pointsl)
avg_t_df = pd.DataFrame(results_df['Avg_T'].tolist(), columns=time_pointsl)

avg_c_means = avg_c_df.mean()
avg_c_stds = avg_c_df.std().fillna(0)
avg_t_means = avg_t_df.mean()
avg_t_stds = avg_t_df.std().fillna(0)
results_df['Avg_C_Zs'] = avg_c_df.apply(
    lambda row: calculate_z_scores_per_timepoint(row, avg_c_means, avg_c_stds), axis=1)
results_df['Avg_T_Zs'] = avg_t_df.apply(
    lambda row: calculate_z_scores_per_timepoint(row, avg_t_means, avg_t_stds), axis=1)

# Creating subplots
num_plots=4 if filter_data else 2
w_rat=[2,2,2,0.5] if filter_data else [1,0.5]
figusize=(16,12) if filter_data else (8,9)
fig, axs = plt.subplots(n_clusters, num_plots, figsize=figusize,width_ratios=w_rat)
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05, wspace=0.1, hspace=0.5)

for cluster in range(n_clusters):
    cluster_data = filtered_results_df[filtered_results_df['Cluster'] == cluster]
    threshold = cluster_data["Deviation_Score"].quantile(0.70)
    filtered_cluster_data = cluster_data[cluster_data["Deviation_Score"] <= threshold]
    filtered_cluster_data = filtered_cluster_data.sort_values(by='Deviation_Score', ascending=False)
    # Normalize Deviation_Score for color mapping
    norm = Normalize(vmin=filtered_cluster_data['Deviation_Score'].min(),
                     vmax=filtered_cluster_data['Deviation_Score'].max())
    sm = ScalarMappable(cmap='Blues_r', norm=norm)
    sm.set_array([])
    sm2 = ScalarMappable(cmap='Reds_r', norm=norm)
    sm.set_array([])
    sm3 = ScalarMappable(cmap='Greens_r', norm=norm)
    sm.set_array([])

    for _, row in filtered_cluster_data.iterrows():
        pway = _
        results_df.loc[pway, "Sig_cluster"] = cluster
        deviation_score = row['Deviation_Score']
        alpha = 0.8
        color = sm.to_rgba(deviation_score)
        color2 = sm2.to_rgba(deviation_score)
        color3 = sm3.to_rgba(deviation_score)
        #c_plots = row[col_name_to_cluster]#row['Avg_C']#row['Avg_C_Zs']#normalize_list(row['Avg_C'])#row['Avg_C_Zs']#normalize_list(row['Avg_C_Zs'])#
        #t_plots = row[col_name_to_cluster]#row['Avg_T']#row['Avg_T_Zs']#normalize_list(row['Avg_T'])#row['Avg_T_Zs']#normalize_list(row['Avg_T_Zs'])#
        axs[cluster, 0].plot(time_pointsl, row[col_name_to_cluster], color=color, alpha=alpha)
        axs[cluster, 0].set_xticks(time_pointsl)
        axs[cluster, 0].set_xticklabels(time_pointsl)
        plt.setp(axs[cluster, 0].get_xticklabels(), rotation=90)
        if filter_data:
            if pway in sig_kos:
                c_plots = row[col_name_to_cluster]
                axs[cluster, 1].plot(time_pointsl, c_plots, color=color2, alpha=alpha)
                results_df.loc[pway, "TS_cluster"]= cluster
            if pway in sig_kos3:
                t_plots = row[col_name_to_cluster]
                axs[cluster, 2].plot(time_pointsl, t_plots, color=color3, alpha=alpha)
                results_df.loc[pway, "IS_cluster"]= cluster
            axs[cluster, 1].set_xticks(time_pointsl)
            axs[cluster, 1].set_xticklabels(time_pointsl)
            plt.setp(axs[cluster, 1].get_xticklabels(), rotation=90)
            axs[cluster, 2].set_xticks(time_pointsl)
            axs[cluster, 2].set_xticklabels(time_pointsl)
            plt.setp(axs[cluster, 2].get_xticklabels(), rotation=90)

    axs[cluster, 0].plot(time_pointsl, kmeans.cluster_centers_[cluster], color='black', label='Cluster Mean',
                         linewidth=2)
    cbar = plt.colorbar(sm, ax=axs[cluster, 0], orientation='vertical')

    if filter_data:
        axs[cluster, 1].plot(time_pointsl, kmeans.cluster_centers_[cluster], color='black', label='Cluster Mean',
                             linewidth=2)
        axs[cluster, 2].plot(time_pointsl, kmeans.cluster_centers_[cluster], color='black', label='Cluster Mean',
                             linewidth=2)
    # Create a colorbar for the subplot

        cbar2 = plt.colorbar(sm2, ax=axs[cluster, 1], orientation='vertical')
        cbar3 = plt.colorbar(sm3, ax=axs[cluster, 2], orientation='vertical')

        #Just 3rd Coloumn:
    if cluster==3:
        if filter_data:
            cbar3.set_label('Average deviation from the mean')
            axs[cluster, 1].set_ylabel('Log2FC')
            axs[cluster, 2].set_ylabel('Log2FC')
        axs[cluster, 0].set_ylabel('Log2FC')

    if cluster==0:
        axs[cluster, 0].set_title(f'Log2 Fold Change Clustering')
        if filter_data:
            axs[cluster, 1].set_title(f'Treatment Significant')
            axs[cluster, 2].set_title(f'Interaction Significant')
    #else: axs[cluster, 1].set_title(f'Cluster {cluster}')

    if cluster == n_clusters-1:
        if filter_data:
            #axs[cluster,0].set_xlabel('Time Points (Hours)')
            axs[cluster, 1].set_xlabel('Time Points (Hours)')
            #axs[cluster, 2].set_xlabel('Time Points (Hours)')
        else:
            axs[cluster, 0].set_xlabel('Time Points (Hours)')
    else:
        axs[cluster, 0].set_xticklabels([])
        if filter_data:
            axs[cluster, 1].set_xticklabels([])
            axs[cluster, 2].set_xticklabels([])
    if not filter_data:
        title_line1 = f'Cluster {cluster}\n({len(filtered_cluster_data.index)} Members)'
    else:
        title_line1 = (f'Cluster {cluster}\n'
                       f'Members:\n'
                       f'Unfiltered List: {len(filtered_cluster_data.index)}\n'
                       f'Significant Treatmet: {len(list(set(filtered_cluster_data.index) & set(sig_kos)))}\n'
                       f'Significant Interaction: {len(list(set(filtered_cluster_data.index) & set(sig_kos3)))}')
    axs[cluster, 0].set_xlim(-5,250)
    axs[cluster, 1].set_xlim(-5, 250)
    axs[cluster, 2].set_xlim(-5, 250)
    axs[cluster, -1].text(0.5, 0.5, title_line1, ha='center', va='center', fontsize=12)
    axs[cluster, -1].axis('off')


plt.tight_layout()
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize='small')
plt.savefig(f'{output_file_pre}_{n_clusters}_clusters.png')

#plt.show()
plt.close()


results_df[col_name_to_cluster] = results_df[col_name_to_cluster].apply(
    lambda x: ";".join(str(x).replace('[', '').replace(']', '').split()))
results_df['Avg_C_Zs']=results_df['Avg_C_Zs'].apply(
    lambda x: ";".join(str(x).replace('[', '').replace(']', '').replace(',', ' ').split()))
results_df['Avg_T_Zs']=results_df['Avg_T_Zs'].apply(
    lambda x: ";".join(str(x).replace('[', '').replace(']', '').replace(',', ' ').split()))
results_df['Avg_C']=results_df['Avg_C'].apply(
    lambda x: ";".join(str(x).replace('[', '').replace(']', '').replace(',', ' ').split()))
results_df['Avg_T']=results_df['Avg_T'].apply(
    lambda x: ";".join(str(x).replace('[', '').replace(']', '').replace(',', ' ').split()))
results_df.to_csv(f"{output_file_pre}_{n_clusters}_final.tsv", sep="\t",index=True)
