#!/usr/bin/env python
# coding: utf-8

import re
import numpy as np
import pandas as pd
import multiprocessing as mp
import glob
from Bio import SeqIO
import sys
import csv
import collections
import itertools
import difflib


def search(accession, gb_vrl):
    
    vers = 1
    while True:
        try:
            entry = gb_vrl[accession + '.' + str(vers)].description
            return entry
        except KeyError:
            if vers >= 10:
                entry = 'Null'
                return entry
            else:
                vers = vers + 1


def quest (question, mode):
    if mode == 1:    
        while True:
            rem = input(question)
            if rem == "y":
                return rem
            elif rem == "n":
                return rem
            else:
                print("Wrong input. [y/n]")
    else:
        rem = "y"
        return rem


def genome(df, A_dict, B_dict, C_dict, gb_vrl, mode):
    
    for row in df.itertuples(index=False, name=None):
        
        row = list(map(str, row))
        
        if row[8] == 'Pass':
            entry = search(row[0][1:], gb_vrl).upper()
            
            if entry != 'Null':   
                if row[4] == 'A':
                    seg_dict = A_dict
                elif row[4] == 'B':
                    seg_dict = B_dict
                elif row[4] == 'C':
                    seg_dict = C_dict
                else:
                    print("Error!") 
                
                try:
                    for num_key in ["1", "2", "3", "4", "5", "6", "7", "8"]:
                        if re.search(rf"SE(T|G)MENT[: ]+{num_key}", entry):
                            num = num_key
                            if num:
                                break
                    
                    if 'num' not in locals():
                        for seg_key in seg_dict:
                            if re.search(rf"{seg_key}", entry):
                                num = seg_dict.get(seg_key)
                                if num:
                                    break
                        
                        if 'num' not in locals():   
                            rem = quest(row[0][1:] 
                                        + ": No protein information found. Remove? [y/n].\nIRD entry: segment number: " 
                                        + row[2] + ".\nGenBank entry: " + entry + ".", mode)
                            if rem == "y":
                                print(row[0][1:] + ": Removed.")  
                            elif rem == "n":
                                yield row[0] + "|" + row[1] + "|" + row[2] + "|" + row[3] + "|" + row[4] + "|" + row[5] + "|" + row[6]  + "|" + row[7] + "|" + row[8] + ";" + row[9]
                        
                        else:
                            if num == row[2]:
                                yield row[0] + "|" + row[1] + "|" + row[2] + "|" + row[3] + "|" + row[4] + "|" + row[5] + "|" + row[6]  + "|" + row[7] + "|" + row[8] + ";" + row[9]
                            else:
                                rem = quest(row[0][1:] 
                                            + ": Wrong segment number. Change? [y/n].\nIRD entry: segment number: " 
                                            + row[2] + ".\nGenBank entry segment number: " + num + ".", mode)
                                if rem == "y":
                                    row[2] = num
                                    print(row[0][1:] + ": Segment number changed.")
                                    yield row[0] + "|" + row[1] + "|" + row[2] + "|" + row[3] + "|" + row[4] + "|" + row[5] + "|" + row[6]  + "|" + row[7] + "|" + row[8] + ";" + row[9]
                                elif rem == "n":
                                    yield row[0] + "|" + row[1] + "|" + row[2] + "|" + row[3] + "|" + row[4] + "|" + row[5] + "|" + row[6]  + "|" + row[7] + "|" + row[8] + ";" + row[9]
                    
                    else:
                        if num == row[2]:
                            yield row[0] + "|" + row[1] + "|" + row[2] + "|" + row[3] + "|" + row[4] + "|" + row[5] + "|" + row[6]  + "|" + row[7] + "|" + row[8] + ";" + row[9]
                        else:
                            rem = quest(row[0][1:] 
                                        + ": Wrong segment number. Change? [y/n].\nIRD entry: segment number: " 
                                        + row[2] + ".\nGenBank entry segment number: " + num + ".", mode)
                            if rem == "y":
                                row[2] = num
                                print(row[0][1:] + ": Segment number changed.")  
                                yield row[0] + "|" + row[1] + "|" + row[2] + "|" + row[3] + "|" + row[4] + "|" + row[5] + "|" + row[6]  + "|" + row[7] + "|" + row[8] + ";" + row[9]
                            elif rem == "n":
                                yield row[0] + "|" + row[1] + "|" + row[2] + "|" + row[3] + "|" + row[4] + "|" + row[5] + "|" + row[6]  + "|" + row[7] + "|" + row[8] + ";" + row[9]
                
                except:
                    rem = quest(row[0][1:] 
                                + ": No protein information found. Remove? [y/n].\nIRD entry: segment number: " 
                                + row[2] + ".\nGenBank entry: " + entry + ".", mode)
                    if rem == "y":
                        print(row[0][1:] + ": Removed.")  
                    elif rem == "n":
                        yield row[0] + "|" + row[1] + "|" + row[2] + "|" + row[3] + "|" + row[4] + "|" + row[5] + "|" + row[6]  + "|" + row[7] + "|" + row[8] + ";" + row[9]
            
            else:
                rem = quest(row[0][1:] 
                            + ": No GenBank entry found. Remove? [y/n]. \nIRD entry: segment number: " 
                            + row[2] + ".", mode)
                if rem == "y":
                    print(row[0][1:] + ": Removed.")  
                elif rem == "n":
                    yield row[0] + "|" + row[1] + "|" + row[2] + "|" + row[3] + "|" + row[4] + "|" + row[5] + "|" + row[6]  + "|" + row[7] + "|" + row[8] + ";" + row[9]
        
        else:
            #rem = quest(row[0][1:] 
            #            + ": Not passed curation. Remove? [y/n].", mode)
            #if rem == "y":
            #    print(row[0][1:] + ": Removed.")  
            #elif rem == "n":
            #    yield row[0] + "|" + row[1] + "|" + row[2] + "|" + row[3]
            pass
            
            #to many not curation pass so omit this dialog
            
        if 'num' in locals():
            del num
        if 'entry' in locals():
            del entry
        if 'rem' in locals():
            del rem


def main():
         
    infile = sys.argv[1]
    outfile = sys.argv[2]
    gene = sys.argv[3]
    switch = int(sys.argv[4])
    chunk = 10000
        
    A_dict = {"[( ]P[B]?2[ .,)]":"1", "[ ]?POLYMERASE( PROTEIN| GENE)?( BASIC| B)( PROTEIN| GENE)? 2[,. ]":"1", 
              "[( ]P[B]?1[ .,)]":"2", "[ ]?POLYMERASE( PROTEIN| GENE)?( BASIC| B)( PROTEIN| GENE)? 1[,. ]":"2", 
              "[( ]PA[ .,)]":"3", "[ ]?POLYMERASE( PROTEIN| GENE)?( ACIDIC| A| ACID)( PROTEIN| GENE)?[,. ]":"3", 
              "[( ]HA[ .,)]":"4", "[ ]?H[A]?EMAGGLUTININ[,. ]":"4", 
              "[( ]NP[ .,)]":"5", "[ ]?NUCLEOPROTEIN[,. ]":"5", "[ ]?NUCLEOCAPSID[,. ]":"5",
              "[( ]NA[ .,)]":"6", "[ ]?NEURAMINIDASE[,. ]":"6",
              "[( ]M[12]?[ .,)]":"7", "[ ]?MATRIX[., ]":"7", "[ ]MEMBRANE?[., ]":"7",
              "[( ]NS[1]?[ .,)]":"8", "[( ]NEP[ .,)]":"8", "[ ]?NON-STRUCTURAL[., ]":"8", "[ ]?NONSTRUCTURAL[., ]":"8", "[ ]?NUCLEAR EXPORT[., ]":"8"}
    
    B_dict = {"[( ]P[B]?1[ .,)]":"1", "[ ]?POLYMERASE( PROTEIN| GENE)?( BASIC| B)( PROTEIN| GENE)? 1[,. ]":"1", 
              "[( ]P[B]?2[ .,)]":"2", "[ ]?POLYMERASE( PROTEIN| GENE)?( BASIC| B)( PROTEIN| GENE)? 2[,. ]":"2", 
              "[( ]PA[ .,)]":"3", "[ ]?POLYMERASE( PROTEIN| GENE)?( ACIDIC| A| ACID)( PROTEIN| GENE)?[,. ]":"3", 
              "[( ]HA[ .,)]":"4", "[ ]?H[A]?EMAGGLUTININ[,. ]":"4", 
              "[( ]NP[ .,)]":"5", "[ ]?NUCLEOPROTEIN[,. ]":"5", "[ ]?NUCLEOCAPSID[,. ]":"5",
              "[( ]NA[ .,)]":"6", "[ ]?NEURAMINIDASE[,. ]":"6",
              "[( ][B]?M[12]?[ .,)]":"7", "[ ]?MATRIX[., ]":"7", "[ ]MEMBRANE?[., ]":"7",
              "[( ]NS[12]?[ .,)]":"8", "[( ]NEP[ .,)]":"8", "[ ]?NON-STRUCTURAL[., ]":"8", "[ ]?NONSTRUCTURAL[., ]":"8", "[ ]?NUCLEAR EXPORT[., ]":"8"}
    
    C_dict = {"[( ]P[B]?2[ .,)]":"1", "[ ]?POLYMERASE( PROTEIN| GENE)?( BASIC| B)( PROTEIN| GENE)? 2[,. ]":"1", 
              "[( ]P[B]?1[ .,)]":"2", "[ ]?POLYMERASE( PROTEIN| GENE)?( BASIC| B)( PROTEIN| GENE)? 1[,. ]":"2", 
              "[( ]P[A3][ .,)]":"3", "[ ]?POLYMERASE( PROTEIN| GENE)?( ACIDIC| A| ACID)( PROTEIN| GENE)?[,. ]":"3", 
              "[( ]H[AE][F]?[ .,)]":"4", "[ ]?H[A]?EMAGGLUTININ[,. ]":"4", 
              "[( ]NP[ .,)]":"5", "[ ]?NUCLEOPROTEIN[,. ]":"5", "[ ]?NUCLEOCAPSID[,. ]":"5",
              "[( ][C]?M[12]?[ .,)]":"6", "[ ]?MATRIX[., ]":"6", "[ ]MEMBRANE?[., ]":"6",
              "[( ]NS[1]?[ .,)]":"7", "[( ]NEP[ .,)]":"7", "[ ]?NON-STRUCTURAL[., ]":"7", "[ ]?NONSTRUCTURAL[., ]":"7", "[ ]?NUCLEAR EXPORT[., ]":"7"}
    
    files = glob.glob(gene + "gbvrl*.seq")
    gb_vrl = SeqIO.index_db(gene + "gbvrl.idx", files, "genbank")
    reader = pd.read_csv(infile, chunksize = chunk, sep = '|', na_filter = False, header=None)
    
    for df in reader:
        
        #correction mode einfügen als input
        result = genome(df, A_dict, B_dict, C_dict, gb_vrl, switch)

        with open(outfile, 'a', newline='') as file:

            for out in result:

                out_write = csv.writer(file)
                out_write.writerow([out])


if __name__ == '__main__':
    
    main()
