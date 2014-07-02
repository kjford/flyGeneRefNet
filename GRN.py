'''
GRN.py

Construct GeneRefNet:
Network analysis of drosophila gene references in published research papers
uses local MySQL database (see buildGRNdb.py) to build

Copyright Kevin Ford (2014)

dependencies (many, you'll probably want an academic/full Canopy python install):
MySQLdb (http://mysql-python.sourceforge.net/MySQLdb.html)
networkx
sklearn
pandas
numpy
'''
import pandas as pd
import numpy as np
import MySQLdb as mdb
import networkx as nx
import nestedKmeans
from collections import Counter
import time


authent_local = {'host':'localhost',
           'user':'pyuser',
           'passwd':'testpass',
           'db':'flydb'}

def getGeneRefs(auth=authent_local,maxcite=500):
    # Restrict network to genes that have references
    # get all gene, reference pairs as a pandas dataframe
    # limit max number of genes in a reference to avoid references that
    # list every gene they screened/profiled which isn't useful
    sql='''
    SELECT flygenes.name AS symb, flygenes.FBid as FBid, generefs.pubid as pubid
    FROM generefs
    LEFT OUTER JOIN flygenes ON generefs.FBid=flygenes.FBid
    '''
    # use this with flatter table layout
    sql='''
    SELECT symb, FBid, pubid
    FROM generefs
    '''
    
    con=mdb.connect(**auth)
    print('Loading...')
    df=pd.io.sql.read_frame(sql,con)
    print('Done.')
    print('Cleaning out references with too many gene references...')
    todelete=[]
    for pub,gene in df.groupby('pubid'):
        # leave out publications that have more than maxcite genes
        # as these are usually not informative
        if len(gene)>maxcite:
            todelete.append(pub)
    for i in todelete:
        df=df.drop(df.index[df.pubid==i])
    return df


def getnedict(df):
    # create node-edge dictionary with weights
    # initialize dictionary of genes (nodes)
    genedict={}
    ugenes=df['symb'].unique()
    for g in ugenes:
        genedict[g]=[]
    for pub,gene in df.groupby('pubid'):
        # all genes in a publication update weights
        gl=gene.symb.tolist()
        for i in gl:
            toadd=gl[:]
            toadd.pop(gl.index(i))
            genedict[i].append(toadd)
                
    # loop back through genes and tally as weights
    for g in genedict.iterkeys():
        gl=genedict[g]
        minidict={}
        gl= [item for sublist in gl for item in sublist]
        for i in np.unique(gl):
            minidict[i]={'weight':gl.count(i)}
        # update genedict
        genedict[g]=minidict
    return genedict


def makeJaccardMat(df):
    # create similarity matrix from dataframe using Jaccard similarity index
    # which is A n B / A u B
    # A=[aij] [0,1]
    # note that the row and column ids can be retreived as df['symb'].unique()
    # with is sorted
    ugenes=df['symb'].unique()
    genedict={}
    for gene,pub in df.groupby('symb'):
        # add pubid's to gene
        genedict[gene]=set(pub.pubid.tolist())
    A=np.zeros((len(ugenes),len(ugenes)))
    for i in range(len(ugenes)):
        for j in range(i):
            jacind=len(set.intersection(genedict[ugenes[i]],genedict[ugenes[j]]))/ (1.0*len(\
            set.union(genedict[ugenes[i]],genedict[ugenes[j]])))
            A[i,j]=jacind
            A[j,i]=jacind
    return A


def makeNetwork_fromdict(genedict):
    # create a networkx network from the dictionary of gene references
    G=nx.from_dict_of_dicts(genedict)
    return G

def makeNetwork_fromAdj(df,Adj):
    # create a networkx network from adjacency array and labels in dataframe
    l=df.symb.unique()
    G=nx.from_numpy_matrix(Adj)
    for n in G.nodes():
        G.node[n]['label']=l[n]
    return G

def clusterNetwork_G(G,k=30):
    # cluster network using spectral clustering on precomputed affinity matrix
    # uses spectral clustering
    A=np.array(nx.to_numpy_matrix(G))
    sc=SpectralClustering(n_clusters=k,affinity='precomputed')
    y=sc.fit_predict(A)
    return y.tolist()


def clustNetwork_Adj(A,TOM,names,thr=15,k=30):
    # cluster network using spectral clustering
    # on Adjacency matrix, using only top thr% of nodes
    # and excluding top 0.1%
    # then assign rest as mode of closest labels of topological overlap matrix
    
    # get top nodes
    N=(A>0).sum(axis=0)
    cutoff=np.percentile(N,100-thr)
    topnodes=(N>cutoff) * (N<np.percentile(N,99.9))
    print('Clustering on %i of %i nodes'%(np.sum(topnodes),len(topnodes)))
    A2=A*1.0*topnodes.reshape(len(topnodes),1)*1.0*topnodes.reshape(1,len(topnodes))
    sc=SpectralClustering(n_clusters=k,affinity='precomputed')
    y=sc.fit_predict(A2)
    print('Done.')
    # labels as names of top genes
    yl=[]
    y=list(y)
    yarr=np.array(y)
    for i in range(k):
        topind=np.argmax(N*(yarr==i))
        if topnodes[topind]:
            yl.append(names[topind])
            for j in range(len(y)):
                if i==y[j]:
                    y[j]=names[topind]
        else:
            yl.append('none') # not in top inds
            for j in range(len(y)):
                if i==y[j]:
                    y[j]='none'
    print(yl)
    print('Assigning rest of labels')
    n0=y.count('none')
    del0=y.count('none')
    while del0:
        # assign labels
        for i in range(len(y)):
            if y[i]=='none':
                # take label of most similarly connected neighbor that isn't none
                sortlabs=np.argsort(-TOM[i,:])
                b=0
                while b>=0:
                    ltest=y[sortlabs[b]]
                    if ltest=='none':
                        if b>len(y):
                            yil='none'
                            b=-1
                        else:
                            b+=1
                    else:
                        yil=ltest
                        b=-1
                y[i]=yil
        del0=n0-y.count('none')
        n0=y.count('none')
    return y
    

def addGraphLabels(G,y,labeltitle):
    # add labels y to graph G
    c=0
    for n in G.nodes():
        G.node[n][labeltitle]=y[c]
        c+=1
    return G
    
def saveGraphML(G,fn):
    # save out network as a graphml file
    nx.write_graphml(G,fn)

if __name__=='__main__':
    timeint = np.int(time.time())
    # load genes and references into dataframe
    print('Loading gene-reference data')
    gr_df = getGeneRefs(auth=authent_local,maxcite=500)
    
    # create adjacency matrix
    print('Making adjacency matrix')
    A = makeJaccardMat(gr_df)
    np.save('savedfiles/adjmat_%i'%timeint,A)
    # cluster on adjacency matrix
    clustmodel=nestedKmeans.nkm(minclust=100,maxklevel=30,maxdepth=4)
    print('Clustering network')
    clustmodel.fit(A)
    L = clustmodel.labels_
    
    # create network X graph
    print('Creating network for visualization')
    G = makeNetwork_fromAdj(gr_df,A)
    
    # add labels
    G = addGraphLabels(G,L,'cluster')
    
    # save out
    print('Saving graph')
    saveGraphML(G,'savedfiles/graphml_%i'%timeint)
    print('Done')
    