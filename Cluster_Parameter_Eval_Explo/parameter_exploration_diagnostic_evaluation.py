#!/usr/bin/env python
# coding: utf-8

# In[2]:


import functools as ft
import numpy as np
import pandas as pd
import multiprocessing as mp
import sys
import re
from decimal import Decimal
import csv
import collections as co
import itertools as it
import umap
import hdbscan
import time 
import random
import progressbar as pb
import scipy.spatial.distance as ssd


# In[3]:


def cast_post_prior(num, step):
    if(type(num) == float):
        post_prior = [
            round(num - step, 1),
            round(num, 1),
            round(num + step, 1)
        ]
        return post_prior
    elif(type(num) == int):
        post_prior = [
            num - step,
            num,
            num + step
        ]
        return post_prior
    else:
        print('Wrong Input Format.')


# In[4]:


class vectorizer(object):
    
    #def __init__(self, procs = 8, k = 7, convert = 0, row = 0, index = [], exist = co.defaultdict(int)):
    def __init__(self, k = 7, convert = 0):
    
        self.k = k
        self.convert = convert
        self.index = [] 
        #self.index = index
        self.exist = co.defaultdict(int) 
        #self.exist = exist
        self.keys = list(self.exist.keys())
        self.col = len(self.keys)
        self.row = 0
        #self.row = row
        self.matrix = np.empty((self.row, self.col, ),dtype = "float32")
        self.amino = co.defaultdict(str, {
            'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAT':'N',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AGA':'R', 'AGC':'S', 'AGG':'R', 'AGT':'S',
            'ATA':'I', 'ATC':'I', 'ATG':'M', 'ATT':'I',
            'CAA':'Q', 'CAC':'H', 'CAG':'Q', 'CAT':'H',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'GAA':'E', 'GAC':'D', 'GAG':'E', 'GAT':'D',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',    
            'TAA':'Y', 'TAC':'*', 'TAG':'*', 'TAT':'Y',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TGA':'*', 'TGC':'C', 'TGG':'W', 'TGT':'C',
            'TTA':'L', 'TTC':'F', 'TTG':'L', 'TTT':'F'
        })
                
    def translate(self, read):
    
        chain = ''

        for i in range(len(read) - 2):
            trip = read[i:i+3]
            chain += self.amino[trip]

        return(chain)
    
    
    def adjust_to_data(self, infile, skip = 0):
    
        #widgets = ['Data adjustment:    ', pb.AnimatedMarker()] 
        #bar = pb.ProgressBar(widgets=widgets).start() 
        data = pd.read_csv(infile, skiprows=skip, chunksize = 10000, sep = ';', na_filter = False, header = None)
        self.row = 0
        #t = 0
        
        for split in data:

            for info, read in split.itertuples(index=False, name=None):
                
                name = info.split('|')[0][1:]

                self.index.append(name)
                del name
                
                if self.convert == 1:
                    #seq = self.translate(re.sub('[^ACGT]+', '', read))
                    seq = self.translate(read)
                    del read
                    
                    num = len(seq) - self.k + 1
                    
                    for i in range(num):
                        kmer = seq[i:i+self.k]
                        self.exist[kmer] = 0
                    
                else:
                    #seq = re.sub('[^ACGT]+', '', read)
                    seq = read
                    del read

                    num = len(seq) - self.k + 1

                    if re.match('^[ACGT]*$', seq): 
                        for i in range(num):
                            kmer = seq[i:i+self.k]
                            self.exist[kmer] = 0
                    else:
                        for i in range(num):
                            kmer = seq[i:i+self.k]
                            if re.match('^[ACGT]*$', kmer): 
                                self.exist[kmer] = 0
                
                #t = t + 1
                self.row = self.row + 1
                #bar.update(t)
            
        self.keys = list(self.exist.keys())
        self.col = len(self.keys)
        self.matrix = np.empty((self.row, self.col, ), dtype="float32")
        
        #bar.finish()
        del seq
    
    def calculate_frequence(self, infile, skip = 0):
        
        #widgets = [' [', pb.Timer(format= 'Vector calculation: %(elapsed)s'), '] ', pb.Bar('*'),' (', pb.ETA(), ') ', ] 
        #widgets = [pb.Timer(format= 'Vector calculation: '), pb.Bar('#'),' (', pb.ETA(), ')'] 
        #bar = pb.ProgressBar(max_value=self.row, widgets=widgets).start() 
        n = 0
        #t = 0
        data = pd.read_csv(infile, skiprows=skip, chunksize = 10000, sep = ';', na_filter = False, header = None)
        
        for split in data:

            for info, read in split.itertuples(index=False, name=None):

                if self.convert == 1:
                    #seq = self.translate(re.sub('[^ACGT]+', '', read))
                    seq = self.translate(read)
                    del read
                
                    counts = self.exist.copy()
                    num = len(seq) - self.k + 1

                    for i in range(num):
                        kmer = seq[i:i+self.k]
                        counts[kmer] += 1
                            
                else:
                    #seq = re.sub('[^ACGT]+', '', read)
                    seq = read
                    del read
                
                    counts = self.exist.copy()
                    num = len(seq) - self.k + 1

                    if re.match('^[ACGT]*$', seq): 
                        for i in range(num):
                            kmer = seq[i:i+self.k]
                            counts[kmer] += 1
                    else:
                        for i in range(num):
                            kmer = seq[i:i+self.k]
                            if re.match('^[ACGT]*$', kmer): 
                                counts[kmer] += 1

                vector = np.array(list(counts.values()), dtype = "float32")/num
                self.matrix[n] = vector
                
                n = n + 1
                #t = t + 1
                #bar.update(t)
                
                counts.clear()
                del vector
                del seq
                del counts

        #bar.finish()
            
    def get_index(self):
        
        return(self.index)
    
    
    def get_keys(self):
        
        return(self.keys)
    
    
    def get_matrix(self):
        
        return(self.matrix)


# In[5]:


class parameter(object):
    
    def __init__(self, infile, outfile, V, W, X, Y, A, B, C, skip = 0, text = 'Parameter exploration: '):
        
        self.infile = infile
        self.outfile = outfile
        self.skip = skip
        self.text = text
        self.V = V
        self.W = W
        self.X = X
        self.Y = Y
        self.A = A
        self.B = B
        self.C = C
        self.length_U = len(self.V)*len(self.W)*len(self.X)*len(self.Y)
        self.length_H = len(self.A)*len(self.B)*len(self.C)
        self.res = np.empty((self.length_U * self.length_H, 3, ), dtype="float32")
        self.index = []
        
    def exploration(self):
        
        widgets = [pb.Timer(format= self.text), pb.Bar('#'),' (', pb.ETA(), ')'] 
        bar = pb.ProgressBar(max_value=self.length_U, widgets=widgets).start() 
        z=0

        for v in self.V:
            for w in self.W:
                for x in self.X:
                    for y in self.Y:

                        if v <= 2:
                            v2 = 2
                        else:
                            v2 = v
                            
                        if w <= 1:
                            w2 = 1
                        else:
                            w2 = w
                        
                        if x <= 2:
                            x2 = 2
                        else: 
                            x2 = x
                        
                        if y <= 1:
                            y2 = 1
                        else:
                            y2 = y
                        
                        try:
                            freq_nt = vectorizer(k = 7, convert = 0)
                            freq_nt.adjust_to_data(self.infile, self.skip)
                            freq_nt.calculate_frequence(self.infile, self.skip)

                            matrix_nt = freq_nt.get_matrix()
                            index_nt = freq_nt.get_index()   
                            keys_nt = freq_nt.get_keys()

                            del freq_nt

                            matrix_nt_red = umap.UMAP(
                                n_neighbors = v2,
                                min_dist = w2,
                                random_state = 42,
                                n_components = 20,
                                metric = 'cosine',
                            ).fit_transform(matrix_nt)

                            del matrix_nt

                            freq_aa = vectorizer(k = 5, convert = 1)
                            freq_aa.adjust_to_data(self.infile, self.skip)
                            freq_aa.calculate_frequence(self.infile, self.skip)

                            matrix_aa = freq_aa.get_matrix()
                            index_aa = freq_aa.get_index()
                            keys_aa = freq_aa.get_keys()

                            del freq_aa

                            matrix_aa_red = umap.UMAP(
                                n_neighbors = x2,
                                min_dist = y2,
                                random_state = 42,
                                n_components = 20, #teste mit etwas mehr dimensions bei AS
                                metric = 'cosine',
                            ).fit_transform(matrix_aa)

                            del matrix_aa

                            matrix_aa_ind = pd.DataFrame(matrix_aa_red, index = index_aa)
                            matrix_nt_ind = pd.DataFrame(matrix_nt_red, index = index_nt)

                            matrix = pd.concat([matrix_nt_ind, matrix_aa_ind], axis=1, copy = False, ignore_index = True) #falsches Ergebnis? checken ob ignore_index = Fehler

                            d = 0
                            for a in self.A:
                                for b in self.B:
                                    for c in self.C:

                                        pos = ( z * self.length_H ) + d

                                        if a <= 1:
                                            a2 = 1
                                        else:
                                            a2 = a
                                            
                                        if b <= 2:
                                            b2 = 2
                                        else:
                                            b2 = b
                                        
                                        if c <= 0:
                                            c2 = 0
                                        else:
                                            c2 = c

                                        matrix_clust = hdbscan.HDBSCAN(
                                            min_samples = a2, #larger the value the more conservative the clustering (more points will be declared as noise)
                                            min_cluster_size = b2, #minimum size that can become a cluster
                                            cluster_selection_epsilon = c2, #don't seperate clusters with a distance less than value
                                            alpha = 1.0, #don't mess with this
                                        ).fit(matrix)

                                        clusterlabel = matrix_clust.labels_

                                        clusters = pd.DataFrame(zip(clusterlabel, ['false'] * len(clusterlabel)), index = index_nt, columns = ['cluster', 'centroid'])

                                        num = clusters['cluster'].max()+1
                                        values = ['true']*num
                                        accessions = []

                                        inner=0
                                        for i in range(num):

                                            query = clusters[clusters.cluster == i]
                                            match = query.index.values.tolist()
                                            sub = matrix.filter(items = match, axis=0)
                                            dist = ssd.cdist(sub, sub, metric = 'cosine')
                                            accessions.append(pd.DataFrame(dist, columns = match, index = match, dtype = 'float32').mean().idxmin())
                                            inner = inner + pd.DataFrame(dist, columns = match, index = match, dtype = 'float32').mean().mean()

                                        centroids = pd.DataFrame(values, columns=['centroid'], index = accessions)

                                        clusters.update(centroids)

                                        #outfile = 'Results/' + str(v) + '_' + str(w)  + '_' + str(x)  + '_' + str(y)  + '_' + str(a)  + '_' + str(b)  + '_' + str(c) + '.csv'

                                        #stop = time.perf_counter()

                                        diagnostic = co.Counter(clusterlabel)

                                        pd.DataFrame([[str(diagnostic[-1]), str(len(set(diagnostic))), str(inner/num)]], index = [str(v) + '_' + str(w)  + '_' + str(x)  + '_' + str(y)  + '_' + str(a)  + '_' + str(b)  + '_' + str(c)], columns = ['n_unclustered', 'n_cluster', 'mean_cluster_mean']).to_csv(self.outfile, mode='a', header = False)

                                        self.index.append(str(v) + '_' + str(w)  + '_' + str(x)  + '_' + str(y)  + '_' + str(a)  + '_' + str(b)  + '_' + str(c))

                                        vector = np.array([[str(diagnostic[-1]), str(len(set(diagnostic))), str(inner/num)]], dtype = "float32")

                                        self.res[pos] = vector

                                        d = d + 1
                        
                        except:
                            self.index.extend(['zero']*self.length_H)
                        
                        z = z + 1
                        
                        bar.update(z)

        bar.finish()
        
    
    def evaluation(self):
        
        n = sum(1 for l in open(self.infile))
    
        data = pd.DataFrame(self.res, index=self.index, columns = ['n_unclustered', 'n_cluster', 'mean_cluster_mean'], dtype = 'float32')
        
        try:
            data.drop('zero', inplace=True)
        except:
            pass
        
        data['rating_unclust'] = data['mean_cluster_mean'] * ((data['n_unclustered']+1)/n)
        
        return(data)


# In[6]:


def main():

    ####### Building a solid random Subset #######
    
    infile = sys.argv[1]
    outfile = sys.argv[2]
    
    lines = sum(1 for l in open(infile))
    size = 10000
    skip = random.sample(range(1, lines), lines - size)
    
    ####### Initial Parameter Exploration #######
    
    ex = parameter(infile = infile, outfile = outfile, V = [40, 50, 60], W = [0.0, 0.1, 0.2], X = [10, 20, 30], Y = [0.0, 0.1, 0.2], A = [1, 5, 10], B = [5, 10, 15], C = [0.25, 0.5, 0.75], skip = skip, text = 'Initial Exploration: ')
    ex.exploration()

    res = ex.evaluation()
    best = res.sort_values(by=['rating_unclust']).head(1)
    para = best.index.values.tolist()[0].split("_")
    
    V_new = cast_post_prior(int(para[0]), 10)
    W_new = cast_post_prior(float(para[1]), 0.1)
    X_new = cast_post_prior(int(para[2]), 10)
    Y_new = cast_post_prior(float(para[3]), 0.1)
    A_new = cast_post_prior(int(para[4]), 5)
    B_new = cast_post_prior(int(para[5]), 5)
    C_new = cast_post_prior(float(para[6]), 0.1)
    
    loop = 0
    matching = 0
    
    ####### Iterative Parameter Refinement #######
    
    while matching < 7 and loop < 10:
        
        ###or if statement einfügen
        
        del ex
        del res
        del best
        del para
        
        matching = 0
        
        loop_check = " " * (3 - len(str(loop))) + str(loop) 
        ex = parameter(infile = infile, outfile = outfile, V = V_new, W = W_new, X = X_new, Y = Y_new, A = A_new, B = B_new, C = C_new, skip = skip, text = 'Refinement Loop ' + loop_check + ': ')
        ex.exploration()

        res = ex.evaluation()
        best = res.sort_values(by=['rating_unclust']).head(1)
        para = best.index.values.tolist()[0].split("_")
        
        if V_new[1] == int(para[0]):
            matching += 1
        else:
            V_new = cast_post_prior(int(para[0]), 10)
        
        if W_new[1] == float(para[1]):
            matching += 1
        else:
            W_new = cast_post_prior(float(para[1]), 0.1)

        if X_new[1] == int(para[2]):
            matching += 1
        else:
            X_new = cast_post_prior(int(para[2]), 10)
           
        if Y_new[1] == float(para[3]):
            matching += 1
        else:
            Y_new = cast_post_prior(float(para[3]), 0.1)

        if A_new[1] == int(para[4]):
            matching += 1
        else:
            A_new = cast_post_prior(int(para[4]), 5)
       
        if B_new[1] == int(para[5]):
            matching += 1
        else:
            B_new = cast_post_prior(int(para[5]), 5)
                
        if C_new[1] == float(para[6]):
            matching += 1
        else:
            C_new = cast_post_prior(float(para[6]), 0.1)
        
        loop += 1
      
    ####### Rerun Pipeline with best Set of Parameters #######
    
    if loop == 10:
        print('Given Limit of Iterations reached.')
    elif matching == 7:
        print('Success.')
    
    print('Best Parameter: ' + para[0] + ' ' + para[1] + ' ' + para[2] + ' ' + para[3] + ' ' + para[4] + ' ' + para[5] + ' ' + para[6])


# In[7]:


if __name__ == "__main__":

    main()


# In[ ]:




