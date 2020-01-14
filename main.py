#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 14:01:58 2020


https://towardsdatascience.com/network-analysis-to-quickly-get-insight-into-an-academic-field-with-python-cd891717d547
@author: wenzhao


from Bio import Entrez
from Bio import Medline
from tqdm import tqdm

# Change this email to your email address
Entrez.email = "wenzhao.wei@student.unsw.edu.au"

keyword = "neuroscience"

result = Entrez.read(Entrez.esearch(db="pubmed", retmax=10, term=keyword))
print(
    "Total number of publications that contain the term {}: {}".format(
        keyword, result["Count"]
    )
)

# Fetch all ids but limit the number returned
MAX_COUNT = result["Count"]
result = Entrez.read(
    Entrez.esearch(db="pubmed", retmax=result["Count"], term=keyword)
)

ids = result["IdList"]

batch_size = 200
batches = [ids[x: x + 200] for x in range(0, len(ids), batch_size)]

record_list = []
for batch in tqdm(batches):
    h = Entrez.efetch(db="pubmed", id=batch, rettype="medline", retmode="text")
    records = Medline.parse(h)
    record_list.extend(list(records))
print("Complete.")
"""


'''
visualize the basic info
'''
from collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from ast import literal_eval
import re, random


sns.set_style("white")

#publication_data = pd.read_csv('hippocampus.csv')
publication_data = pd.read_csv("hippocampus.csv",
                               converters={"FAU": lambda x: x.strip("'[]").split("', '")},
                               #skiprows=lambda i: i>0 and random.random() > 0.7
                               )
'''
publication_data = pd.DataFrame(record_list)

publication_data.dropna(subset=['EDAT'], inplace=True)
publication_data["Year"] = (
    publication_data["EDAT"].astype(str).str[0:4].astype(int)
)
publication_data = publication_data.dropna(axis=1, how='all')
'''

plt.figure(figsize=(8, 8), dpi=100)

# Top 10 authors
plt.subplot(2, 2, 1)
authors_flat = [
    author
    for authors in list(publication_data["FAU"].dropna())
    for author in authors# literal_eval(authors)
]

top10authors = pd.DataFrame.from_records(
    Counter(authors_flat).most_common(10), columns=["Name", "Count"]
)
sns.barplot(x="Count", y="Name", data=top10authors, palette="RdBu_r")
plt.title("Top 10 Authors")


# Publications over Time
plt.subplot(2, 2, 2)
yearly = pd.DataFrame(publication_data["Year"].value_counts().reset_index())
yearly.columns = ["Year", "Count"]
sns.lineplot(x="Year", y="Count", data=yearly)
plt.title("Publications over Time")
plt.xlim([1990, 2020])


plt.subplot(2, 2, 3)

# TOP 10 Journals
top10journals = pd.DataFrame.from_records(
    Counter(publication_data["TA"]).most_common(10),
    columns=["Journal", "Count"],
)

sns.barplot(x="Count", y="Journal", data=top10journals, palette="RdBu_r")
plt.title("Top 10 Journals")

# Top associated keywords
plt.subplot(2, 2, 4)

'''
    regularize key words
'''
def tokenize(x):
    tags = re.compile('<.*?>') # remove the html tags
    pattern = r'[^a-zA-z\s]' # clean the non text first 
    x = re.sub(tags, '', x)
    x = re.sub(pattern, '', x)
    return x

flat_kw = [
    tokenize(kw.lower())
    for kws in list(publication_data["OT"].dropna())
    for kw in literal_eval(kws)
    # for _ in kw.split(" ")
]

top10kw = pd.DataFrame.from_records(
    Counter(flat_kw).most_common(10), columns=["Keyword", "Count"]
)

sns.barplot(x="Count", y="Keyword", data=top10kw, palette="RdBu_r")
plt.title("Top 10 Associated Keywords")
#plt.subplots_adjust(top=1, bottom=0, left=0, right=1, hspace=0.3, wspace=0.3)

plt.show()


'''
social network analysis
'''
from itertools import combinations
import networkx as nx
from nxviz.plots import CircosPlot

# Extract author connections
authors = publication_data["FAU"].dropna()
author_connections = list(
    map(lambda x: list(combinations(x[::-1], 2)), authors)
)
flat_connections = [item for sublist in author_connections for item in sublist]

# Create a dataframe with the connections
df_graph = pd.DataFrame(flat_connections, columns=["From", "To"])
df_graph = df_graph.groupby(["From", "To"]).size().reset_index()
df_graph.columns = ["From", "To", "Count"]


G = nx.from_pandas_edgelist(
    df_graph, source="From", target="To", edge_attr="Count"
)

# Limit to TOP 50 authors
top50authors = pd.DataFrame.from_records(
    Counter(authors_flat).most_common(50), columns=["Name", "Count"]
)

top50_nodes = (n for n in list(G.nodes()) if n in list(top50authors["Name"]))

G_50 = G.subgraph(top50_nodes)

for n in G_50.nodes():
    G_50.nodes[n]["publications"] = int(
        top50authors[top50authors["Name"] == n]["Count"]
    )

c = CircosPlot(
    G_50,
    dpi=600,
    node_grouping="publications",
    edge_width="Count",
    figsize=(20, 20),
    node_color="publications",
    node_labels=True,
)
c.draw()
plt.show()


# default plot
pos = nx.spring_layout(G_50)
betCent = nx.betweenness_centrality(G_50)
node_color = [10000.0 * G_50.degree(v) for v in G_50]
node_size =  [v * 10000 for v in betCent.values()]
plt.figure(figsize=(20,20))
nx.draw_networkx(G_50, 
                 pos=pos, 
                 with_labels=True,
                 node_color=node_color,
                 alpha=0.8,
                 font_size=8,
                 node_size=node_size )



# plotly
import igraph as ig
import plotly.graph_objects as go
from plotly.offline import plot

edge_trace = go.Scatter(
    x=[],
    y=[],
    line=dict(width=0.5, color='#888'),
    hoverinfo='none',
    mode='lines')

for edge in G_50.edges():
    x0, y0 = pos[edge[0]]
    x1, y1 = pos[edge[1]]
    edge_trace['x'] += tuple([x0, x1, None])
    edge_trace['y'] += tuple([y0, y1, None]) 


node_trace = go.Scatter(
    x=[],
    y=[],
    text=[],
    mode='markers',
    hoverinfo='text',
    marker=dict(
        showscale=True,
        # colorscale options
        #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
        #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
        #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
        colorscale='YlGnBu',
        reversescale=True,
        color=[],
        size=10,
        colorbar=dict(
            thickness=30,
            title='Node Connections',
            xanchor='left',
            titleside='right'
        ),
        line=dict(width=2)))  
        
for node in G_50.nodes():
    x, y = pos[node]
    node_trace['x'] += tuple([x])
    node_trace['y'] += tuple([y])    
    
# color node points
for node, degree in enumerate(G_50.degree()):
    node_trace['marker']['color']+=tuple([degree[1] * 2])
    node_trace['marker']['size'] = [ 5*c + 10 for c in node_trace['marker']['color'] ]
    node_info = degree[0] + ' has connections: '+str(degree[1])
    node_trace['text']+=tuple([node_info])
'''
node_trace.marker.color = node_adjacencies
node_trace.marker.size = 50
node_trace.text = node_text
'''

fig = go.Figure(data=[edge_trace, node_trace],
             layout=go.Layout(
                title='Neurosicne 70%',
                titlefont=dict(size=16),
                #hovermode='closest',
                margin=dict(b=20,l=5,r=5,t=40),
                annotations=[ dict(
                    text="WWZ",
                    showarrow=False,
                    xref="paper", yref="paper",
                    x=0.005, y=-0.002 ) ],
                xaxis=dict(showgrid=True, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=True, zeroline=False, showticklabels=False)))

#py.iplot(fig, filename='networkx')
    
plot(fig)
# fig.show()


'''
social influences
'''
deg = nx.degree_centrality(G_50)
bet = nx.betweenness_centrality(G_50)
# closeness_centrality = nx.closeness_centrality(G_50)

top_df = pd.DataFrame.from_dict(
    [deg, bet, dict(Counter(authors_flat).most_common(50))]
).T

top_df.columns = [
    "Degree Centrality",
    "Betweenness Centrality",
    "Publications",
]

for col in top_df.columns:
    top_df[col] = top_df[col] / max(top_df[col])

top_df = top_df.sort_values("Publications", ascending=False)[:10]
top_df = pd.DataFrame(top_df.stack())
top_df = top_df.reset_index()
top_df.columns = ["Name", "Parameter", "Value"]


fig, ax = plt.subplots(figsize=(10, 8), dpi=100)

sns.barplot(x="Value", y="Name", data=top_df, hue="Parameter", palette="Blues")

plt.title("Top 10 Authors, normalized parameters")
plt.show()
