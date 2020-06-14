#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
from sklearn.cluster import OPTICS
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA

#%%
df = pd.read_csv("scores_beta2.txt",sep=" , ", )
df.head(10)

#%%
nucleotide_ids = list(set(list(df["Source"])))
nucleotide_ids.sort()

# %%
# The entire idea...
# Allocate a matrix
# Find source and target index from the list
# Fill that position in the matrix

# %%
# Utility function for getting index from source, target pair
def get_index(source, target, data):
    row = data.index(source)
    col = data.index(target)
    return row,col

# %%
# Allocate the matrix
mat = np.zeros((len(nucleotide_ids), len(nucleotide_ids)), dtype='float32')

# %%
# Fill the matrix
for idx, row in df.iterrows():
    source = row["Source"]
    target = row["Target"]
    score = row["Score"]
    r,c = get_index(source, target, nucleotide_ids)
    mat[r][c] = float(score)

# %%
# Scale the matrix
scaler = MinMaxScaler()
scaler.fit(mat)
mat_sc = scaler.transform(mat)

# %%
# # Uncomment for K Means
# NUM_CLUSTERS = 3
# kmeans = KMeans(n_clusters=NUM_CLUSTERS)
# kmeans.fit(mat_sc)
# labels = kmeans.labels_

# %%
# OPTICS is similar to DBSCAN
optics = OPTICS(min_samples=2).fit(mat_sc)
labels = optics.labels_

# %%
# Reduce dimensions for visualization
pca = PCA(n_components=2)
pca.fit(mat_sc)
X=pca.transform(mat_sc)
print(np.sum(pca.explained_variance_ratio_))
df_pca = pd.DataFrame(data=X, columns=['x1', 'x2'])
df_pca["labels"] = labels
df_pca.head()

# %%
# Add the geolocation
df_map = pd.read_csv("map.csv")

def reduce_name(name):
    return str(name).split(":")[0].lower()

df_map["Geo_Location"] = df_map["Geo_Location"].apply(reduce_name)

geo = []
for seq in nucleotide_ids:
    id = seq.split(".")[0]          # remove versioning
    row = df_map.loc[df_map['Accession'] == str(id)]
    region = str(row["Geo_Location"].tolist()[0])
    geo.append(region)

# %%
df_pca["geo"] = geo
print(df_pca.groupby("geo")["geo"].count())
df_plt = df_pca
df_plt["labels"] = df_plt["labels"].apply(lambda x: x+1)

#%%
fig = px.scatter(df_plt, x="x1", y="x2", color="geo",
                 size='labels', hover_data=['labels'])
fig.show()