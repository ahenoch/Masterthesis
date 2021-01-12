#!/usr/bin/env python


import functools as ft
import numpy as np
import pandas as pd
import multiprocessing as mp
import sys
import re
import csv
import collections as co
import itertools as it
import umap
import hdbscan
import time 
import progressbar


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
    
    
    def adjust_to_data(self, infile):
    
        widgets = ['Data adjustment:    ', progressbar.AnimatedMarker()] 
        bar = progressbar.ProgressBar(widgets=widgets).start() 
        data = pd.read_csv(infile, chunksize = 10000, sep = ';', na_filter = False, header = None)
        self.row = 0
        t = 0
        
        for split in data:

            for info, read in split.itertuples(index=False, name=None):
                
                name = info.split('|')[0][1:]

                self.index.append(name)
                del name
                
                if self.convert == 1:
                    seq = self.translate(re.sub('[^ACGT]+', '', read))
                    del read
                else:
                    seq = re.sub('[^ACGT]+', '', read) 
                    del read

                num = len(seq) - self.k + 1
                    
                for i in range(num):
                    kmer = seq[i:i+self.k]
                    self.exist[kmer] = 0
                
                t = t + 1
                self.row = self.row + 1
                bar.update(t)
            
        self.keys = list(self.exist.keys())
        self.col = len(self.keys)
        self.matrix = np.empty((self.row, self.col, ), dtype="float32")
        
        bar.finish()
        del seq
    
    def calculate_frequence(self, infile):
        
        #widgets = [' [', progressbar.Timer(format= 'Vector calculation: %(elapsed)s'), '] ', progressbar.Bar('*'),' (', progressbar.ETA(), ') ', ] 
        widgets = [progressbar.Timer(format= 'Vector calculation: '), progressbar.Bar('#'),' (', progressbar.ETA(), ')'] 
        bar = progressbar.ProgressBar(max_value=self.row, widgets=widgets).start() 
        n = 0
        t = 0
        data = pd.read_csv(infile, chunksize = 10000, sep = ';', na_filter = False, header = None)
        
        for split in data:

            for info, read in split.itertuples(index=False, name=None):

                if self.convert == 1:
                    seq = self.translate(re.sub('[^ACGT]+', '', read))
                    del read
                else:
                    seq = re.sub('[^ACGT]+', '', read) 
                    del read
                
                counts = self.exist.copy()
                num = len(seq) - self.k + 1
                
                for i in range(num):
                    kmer = seq[i:i+self.k]
                    counts[kmer] += 1

                vector = np.array(list(counts.values()), dtype = "float32")/num
                self.matrix[n] = vector
                
                n = n + 1
                t = t + 1
                bar.update(t)
                
                counts.clear()
                del vector
                del seq
                del counts

        bar.finish()
            
    def get_index(self):
        
        return(self.index)
    
    
    def get_keys(self):
        
        return(self.keys)
    
    
    def get_matrix(self):
        
        return(self.matrix)


def main():
    
    infile = sys.argv[1]
    outfile = sys.argv[2]

    freq_nt = vectorizer(k = 7, convert = 0)
    freq_nt.adjust_to_data(infile)
    freq_nt.calculate_frequence(infile)

    matrix_nt = freq_nt.get_matrix()
    index_nt = freq_nt.get_index()   
    keys_nt = freq_nt.get_keys()

    del freq_nt

    print('Running UMAP.')
    
    matrix_nt_red = umap.UMAP(
        n_neighbors = 50,
        min_dist = 0.25,
        n_components = 20,
        random_state = 42,
        metric = 'cosine',
    ).fit_transform(matrix_nt)

    del matrix_nt

    freq_aa = vectorizer(k = 5, convert = 1)
    freq_aa.adjust_to_data(infile)
    freq_aa.calculate_frequence(infile)

    matrix_aa = freq_aa.get_matrix()
    index_aa = freq_aa.get_index()
    keys_aa = freq_aa.get_keys()

    del freq_aa

    print('Running UMAP.')
    
    matrix_aa_red = umap.UMAP(
        n_neighbors = 50,
        min_dist = 0.25,
        n_components = 20,
        random_state = 42,
        metric = 'cosine',
    ).fit_transform(matrix_aa)

    del matrix_aa

    matrix_aa_ind = pd.DataFrame(matrix_aa_red, index = index_aa)
    matrix_nt_ind = pd.DataFrame(matrix_nt_red, index = index_nt)

    matrix = pd.concat([matrix_nt_ind, matrix_aa_ind], axis=1)

    print('Running HDBscan.')
    
    matrix_clust = hdbscan.HDBSCAN().fit(matrix)

    clusterlabel = matrix_clust.labels_

    cluster = pd.DataFrame(clusterlabel, index = index_nt, columns = ['cluster'])
    pd.DataFrame(cluster).to_csv(outfile, index_label='accession', index=True, header=True, sep='\t')
    
    print('Finished.')


if __name__ == "__main__":

	main()
