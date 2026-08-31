# StringDB

import requests
import pandas as pd
import numpy as np
from scipy import stats
import numpy as np
from sklearn.neighbors import kneighbors_graph
import igraph as ig
import leidenalg
import matplotlib.pyplot as plt
from kneed import KneeLocator
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib import cm


__all__ = [
    "echo_to_format",
    "fdrsubset",
    "compute_modularity",
    "clustering_rawdata",
    "stringdb_analysis",
    "rescale",
    "graph_vis",
]

def echo_to_format(data):
    """
    Method to take data from ECHO and fit it for analysis, filtering for gene names and fitted timepoints only
    """
    data1 = data[[col for col in data if col.startswith('Fitted')]]

    IDs = data['Gene Name'].values

    output = pd.DataFrame(data1, index = IDs)
    return output

def fdrsubset(df1, df2):
  """
  Method to subset our datasets so we are only keeping the data that has FDR values below our acceptable threshold
  Note: make sure the index column is the name of the transcripts
  """
  arr1 = df1.index
  arr2 = df2.index
  set1 = set(arr1)
  set2 = set(arr2)
  overlap = set1.intersection(set2)
  overlap1 = list(overlap)
  indices_array1 = [i for i, element in enumerate(arr1) if element in overlap1]
  indices_array2 = [i for i, element in enumerate(arr2) if element in overlap1]
  sub1 = df1.iloc[indices_array1,:]
  sub2 = df2.iloc[indices_array2,:]

  return sub1, sub2

def compute_modularity(data, K, metric):
    """
    Compute the KNN graph adjacency matrix to
    """
    knn_graph = kneighbors_graph(data, n_neighbors=K, mode='connectivity', metric = metric, include_self=False)
    adjacency_matrix = knn_graph.toarray() #type: ignore 
    g = ig.Graph.Adjacency((adjacency_matrix > 0).tolist())
    partition = leidenalg.find_partition(g, leidenalg.ModularityVertexPartition,)
    modularity = g.modularity(partition)
    return modularity

def clustering_rawdata(df, ids, metric = 'correlation', knn = True, k_comp = 'knee', threshold = 0.05, show_plot=True):
    """
    Function that intakes raw data and clusters proteins based on their expression values using k means and leiden algorithm
    Proteins are assigned cluster labels as a result
    """
    data = np.array(df)

    if knn == True:
        np.random.seed(42)
        n_samples = data.shape[0]

        # Define range of K values
        K_values = [5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100]

        modularity_scores = []
        for K in K_values:
            modularity = compute_modularity(data, K, metric)
            modularity_scores.append(modularity)

        if k_comp == 'knee':
            kneedle = KneeLocator(K_values, modularity_scores, curve='convex', direction='decreasing')
            elbow = kneedle.elbow
            if show_plot:
                plt.figure(figsize=(10, 6))
                plt.plot(K_values, modularity_scores, marker='o')
                plt.axvline(x=elbow, color='r', linestyle='--', label=f'Elbow at K={elbow}') #type: ignore 
                plt.xlabel('K')
                plt.ylabel('Modularity')
                plt.title('Modularity vs K for KNN Graph')
                plt.grid(True)
                plt.legend()
                plt.show()

            K = elbow

        elif type(k_comp) == int:
            K = k_comp.copy()

        knn_graph = kneighbors_graph(data, n_neighbors=K, mode='connectivity', metric = metric, include_self=False) #type: ignore 
        adjacency_matrix = knn_graph.toarray() #type: ignore 
        g = ig.Graph.Adjacency((adjacency_matrix > 0).tolist())
        partition = leidenalg.find_partition(g, leidenalg.ModularityVertexPartition)

    elif knn == False:
        result = stats.spearmanr(data.T, axis=0, alternative='greater')
        corr = result.correlation #type: ignore
        p_mat = result.pvalue #type: ignore
        adjacency_matrix = np.where(p_mat < threshold, 1, 0)
        g = ig.Graph.Adjacency((adjacency_matrix > 0).tolist(), mode="undirected")
        partition = leidenalg.find_partition(g, leidenalg.ModularityVertexPartition)

    labels = np.array(partition.membership, dtype=int)
    label_df = pd.DataFrame([list(ids), labels]).T
    label_df.columns = ['Gene', 'Labels']
    return label_df

## StringDB Function
def stringdb_analysis(cl_df, protein, savefile = 'stringdb.csv', organism: int | str = 'mus musculus', cutoff = 0.05, database = 'all'):
    """
    Performs an analysis by connecting to STRING database and conducting a search for a protein or cluster
    All variable names are CASE SENSITIVE
    --------------
    Parameters:
    cl_df: DataFrame object with 2 columns: first column is gene names, second column is cluster labels corresponding to those gene names
    protein: The protein name of interest as it is in the dataset. Alternatively, enter an integer to automatically perform analysis of a specific cluster with that integer's label
    savefile: Desired name for the output analysis file.  Currently only supports .csv files, openable in excel
    cutoff: FDR threshold to output desired functions.  Set equal to 1 if all functional annotations are desired. 0.05 by default
    organism: Model organism to be studied, can either be scientific name of organism or NCBI taxon ID. Mus musculus by default.
    database: If you only want certain databases, e.g. KEGG or PMID (PubMed).  'all' by default.  Options are 'COMPARTMENTS', 'Component',
        'Function', 'InterPro', 'KEGG', 'Keyword', 'MPO', 'NetworkNeighborAL', 'PMID', 'Pfam', 'Process', 'RCTM', 'SMART', 'TISSUES', 'WikiPathways'
    --------------
    Return: df_export (DataFrame) String DB analysis .csv file for protein/cluster

    """
    ncbi_ids = {
        'mus musculus': 10090,
        'homo sapiens': 9606,
        'drosophila melanogaster': 7227,
        'neurospora crassa': 5141,
        'danio rerio': 7955,
        'arabidopsis thaliana': 3702,
    }

    if type(organism) is int:
      ncbi = organism
    elif organism.casefold() in [key.casefold() for key in ncbi_ids.keys()]: #type: ignore 
      ncbi = ncbi_ids[organism.casefold()] #type: ignore 
    else:
      print("Organism not detected. Please input NCBI Taxon ID or one of the following options: homo sapiens, drosophila melanogaster, neurospora crassa, danio rerio, or arabidopsis thaliana")
      return


    string_api_url = "https://version-12-0.string-db.org/api" #updated to version 12
    output_format = "tsv-no-header"
    method = "enrichment"

    if (type(protein) == int) == True:
        cluster = protein
    else:
        A = np.where(cl_df.iloc[:,0] == protein)[0][0]
        cluster = cl_df.iloc[A,1]

    idx = np.where(cl_df.iloc[:,1] == cluster)[0]
    proteins = list(cl_df.iloc[idx,0].values)

    request_url = "/".join([string_api_url, output_format, method])
    my_genes = proteins.copy()
    params = {

        "identifiers" : "%0d".join(my_genes), # your protein
        "species" : ncbi, # species NCBI identifier
        "caller_identity" : "chuah_lab_frappe" # REPLACE w GITHUB URL

    }

    response = requests.post(request_url, data=params)
    long_string = response.text.strip().split("\n")
    new_data = [s.split('\t') for s in long_string]
    data = pd.DataFrame(new_data, columns=['category', 'term', 'number_of_genes', 'number_of_genes_in_background',
       'ncbiTaxonId', 'inputGenes', 'preferredNames', 'p_value', 'fdr',
       'description'])
    dbs = np.unique(data.iloc[:,0].values) #type: ignore 
    data['fdr'] = data['fdr'].astype(float)
    data_fdr = data[data['fdr'] < cutoff]

    if database != 'all':
        if database not in dbs:
            print('Error, requested database not found in analysis! Available databases for this selection are:'+str(np.unique(data_fdr.iloc[:,0].values))) #type: ignore 
            return
        else:
            db_where = np.where(data_fdr.iloc[:,0] == database)[0]
            df_export = data_fdr.iloc[db_where,:]
    else:
        df_export = data_fdr.copy()

    df_export.to_csv(savefile)
    return df_export

def rescale(l,newmin,newmax):
    arr = list(l)
    return [(x-min(arr))/(max(arr)-min(arr))*(newmax-newmin)+newmin for x in arr]

def graph_vis(protein_list, center_protein, threshold, title, organism_id=10090):
    proteins = '%0d'.join(protein_list)
    url = 'https://string-db.org/api/tsv/network?identifiers=' + proteins + f'&species={organism_id}'
    r = requests.get(url)

    lines = r.text.split('\n')
    data = [l.split('\t') for l in lines]
    df = pd.DataFrame(data[1:-1], columns = data[0])

    interactions = df[['preferredName_A', 'preferredName_B', 'score']]

    # keep rows where the partner meets the threshold OR where Apoe is involved
    interactions = interactions[
        (interactions['score'].astype(float) >= threshold) |   # high scoring interactions
        (interactions['preferredName_A'] == center_protein) | 
        (interactions['preferredName_B'] == center_protein)
    ]

    # then filter down to only interactions involving Apoe
    interactions = interactions[
        (interactions['preferredName_A'] == center_protein) | 
        (interactions['preferredName_B'] == center_protein)
    ].sort_values('score', ascending=False)

    G=nx.Graph(name='Protein Interaction Graph')
    interactions = np.array(interactions)
    for i in range(len(interactions)):
        interaction = interactions[i]
        a = interaction[0]
        b = interaction[1]
        w = float(interaction[2])
        G.add_weighted_edges_from([(a,b,w)])

    def rescale(l,newmin,newmax):
        arr = list(l)
        return [(x-min(arr))/(max(arr)-min(arr))*(newmax-newmin)+newmin for x in arr]

    graph_colormap = cm.get_cmap('plasma', 12)
    c = rescale([G.degree(v) for v in G],0.0,0.9) # type: ignore
    c = [graph_colormap(i) for i in c]
    bc = nx.betweenness_centrality(G)
    s = rescale([v for v in bc.values()],1500,7000)
    ew = rescale([float(G[u][v]['weight']) for u,v in G.edges],0.1,4)
    ec = rescale([float(G[u][v]['weight']) for u,v in G.edges],0.1,1)
    ec = [graph_colormap(i) for i in ec]

    # override center node appearance
    node_list = list(G.nodes())
    center_idx = node_list.index(center_protein)
    s[center_idx] = 10000
    c[center_idx] = (0.6, 0, 1, 1)  # purple

    # pin center node at origin
    fixed_positions = {center_protein: (0, 0)}
    pos = nx.spring_layout(G, pos=fixed_positions, fixed=[center_protein], k=1.5, seed=42)

    plt.figure(figsize=(19,9),facecolor=[0.7,0.7,0.7,0.4]) #type: ignore
    nx.draw_networkx(G, pos=pos, nodelist=node_list, with_labels=True,
                 node_color=c, node_size=s, edge_color='black', width=ew,
                 font_color='white', font_weight='bold', font_size='9')
    plt.title(title, fontsize=16, fontweight='bold', pad=20) 
    plt.axis('off')
    plt.show()
