# Clustering Models 


import pandas as pd
import numpy as np
import leidenalg as la
import igraph as ig
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics.cluster import adjusted_rand_score


__all__ = ["corr_to_adj_graph", "clustering_models"]



def corr_to_adj(corr_matrix, threshold):
  ''' converts correlation matrix to adjacency matrix for given threshold
  '''
  adj_matrix = corr_matrix.applymap(lambda x: 1 if abs(x) > threshold else 0)
  np.fill_diagonal(adj_matrix.values, 0)  # Set diagonal to 0
  return adj_matrix

def clustering_models(): 
  file1 = pd.read_excel('FDSAD_ECHO_INTERSECT.csv') #FDSAD
  file2 = pd.read_excel('MDSLM Intersect 10-18-24 (1).xlsx') #MDSLM
  file3 = pd.read_excel('MDSAD Intersect 10-18-24 (1).xlsx') #MDSAD
  file4 = pd.read_excel('FDSLM Intersect 10-18-24 (1).xlsx') #FDSLM

  df1 = file1[[col for col in file1 if col.startswith('Fitted')]] #FDSAD
  df2 = file2[[col for col in file1 if col.startswith('Fitted')]] #MDSLM
  df3 = file3[[col for col in file1 if col.startswith('Fitted')]] #MDSAD
  df4 = file4[[col for col in file1 if col.startswith('Fitted')]] #FDSLM

  matrix1 = df1.T.corr(method = 'pearson')
  matrix2 = df2.T.corr(method = 'pearson')
  matrix3 = df3.T.corr(method = 'pearson')
  matrix4 = df4.T.corr(method = 'pearson')
  #makes adjacency matrices for all the data, adjustable threshold

  adjacency_matrix1 = corr_to_adj(matrix1, 0.8)
  adjacency_matrix2 = corr_to_adj(matrix2, 0.8)
  adjacency_matrix3 = corr_to_adj(matrix3, 0.8)
  adjacency_matrix4 = corr_to_adj(matrix4, 0.8)
  #starting to implement leidenalg


  adj_matrix_lists1 = np.array(adjacency_matrix1)
  adj_matrix_lists2 = np.array(adjacency_matrix2)
  adj_matrix_lists3 = np.array(adjacency_matrix3)
  adj_matrix_lists4 = np.array(adjacency_matrix4)

  graph_from_list1 = ig.Graph.Adjacency(adj_matrix_lists1)
  graph_from_list2 = ig.Graph.Adjacency(adj_matrix_lists2)
  graph_from_list3 = ig.Graph.Adjacency(adj_matrix_lists3)
  graph_from_list4 = ig.Graph.Adjacency(adj_matrix_lists4)


  #partition the data

  partition1 = la.find_partition(graph_from_list1, la.ModularityVertexPartition)
  partition2 = la.find_partition(graph_from_list2, la.ModularityVertexPartition)
  partition3 = la.find_partition(graph_from_list3, la.ModularityVertexPartition)
  partition4 = la.find_partition(graph_from_list4, la.ModularityVertexPartition)

  ig.plot(partition1)
  #clustering similar

  cluster_labels1 = partition1.membership #makes the clusters in the partition

  unique_elements1 = np.unique(cluster_labels1) #how many clusters are there?
  print(f"# of clusters in set 1:{unique_elements1}")

  cluster_labels2 = partition2.membership

  unique_elements2 = np.unique(cluster_labels2)
  print(f"# of clusters in set 2:{unique_elements2}")

  cluster_labels3 = partition3.membership

  unique_elements3 = np.unique(cluster_labels3)
  print(f"# of clusters in set 3:{unique_elements3}")

  cluster_labels4 = partition4.membership

  unique_elements4 = np.unique(cluster_labels4)
  print(f"# of clusters in set 4:{unique_elements4}")
  # Calculate modularity

  modularity_value1 = graph_from_list1.modularity(cluster_labels1)
  print(f"MSWAD Modularity: {modularity_value1}")

  modularity_value2 = graph_from_list2.modularity(cluster_labels2)
  print(f"MDSLM Modularity: {modularity_value2}")

  modularity_value3 = graph_from_list3.modularity(cluster_labels3)
  print(f"MDSAD Modularity: {modularity_value3}")

  modularity_value4 = graph_from_list4.modularity(cluster_labels4)
  print(f"MSWLM Modularity: {modularity_value4}")

  adjusted_rand_score1_2 = adjusted_rand_score(cluster_labels1, cluster_labels2)
  print(f"Random Index Score of FDSAD & MDSLM: {adjusted_rand_score1_2}")

  adjusted_rand_score1_3 = adjusted_rand_score(cluster_labels1, cluster_labels3)
  print(f"Random Index Score of FDSAD & MDSAD: {adjusted_rand_score1_3}")

  adjusted_rand_score1_4 = adjusted_rand_score(cluster_labels1, cluster_labels4)
  print(f"Random Index Score of FDSAD & FDSLM: {adjusted_rand_score1_4}")

  adjusted_rand_score2_3 = adjusted_rand_score(cluster_labels2, cluster_labels3)
  print(f"Random Index Score of MDSLM & MDSAD: {adjusted_rand_score2_3}")

  adjusted_rand_score2_4 = adjusted_rand_score(cluster_labels2, cluster_labels4)
  print(f"Random Index Score of MDSLM & FDSLM: {adjusted_rand_score2_4}")

  adjusted_rand_score3_4 = adjusted_rand_score(cluster_labels3, cluster_labels4)
  print(f"Random Index Score of MDSAD & FDSLM: {adjusted_rand_score3_4}")
  #elbow method to see how many clusters to use


  X = df1
  list1 = []

  for i in range(1,20):
      kmeans = KMeans(n_clusters=i)
      kmeans.fit(X)
      list1.append(kmeans.inertia_)

  plt.plot(range(1,20), list1, marker='o')
  plt.title('Elbow method')
  plt.xlabel('Number of clusters')
  plt.ylabel('Proteins')
  plt.show()
  #sample code for finding the period values of the clusters, with FDSAD partition data
  #can be replicated depending on how many clusters are in data

  cluster1_array = np.array(cluster_labels1)

  idx_0 = np.where(cluster1_array == 0)
  idx_0_list1 = [i.tolist() for i in idx_0]
  idx_0_list = idx_0_list1[0]
  period_0 = []
  for i in idx_0_list1:
    period_0.append(float(file1.iloc[i, 7]))
  print(f"Period values of cluster 0: {period_0}")
  print(len(period_0))

  #standard deviation of graphs
  std_dev1 = np.std(cluster_labels1)
  print(f"Standard Deviation of FDSAD Clusters: {std_dev1}")

  std_dev2 = np.std(cluster_labels2)
  print(f"Standard Deviation of MDSLM Clusters: {std_dev2}")

  std_dev3 = np.std(cluster_labels3)
  print(f"Standard Deviation of MDSAD Clusters: {std_dev3}")

  std_dev4 = np.std(cluster_labels4)
  print(f"Standard Deviation of FDSLM Clusters: {std_dev4}")
