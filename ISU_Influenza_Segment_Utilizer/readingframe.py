#!/usr/bin/env python


import functools as ft
import numpy as np
import pandas as pd
import multiprocessing as mp
import re
import sys
import csv
import collections
import itertools
import difflib


def search(entry, cities, subcountries, countries, cities_names, subcountries_names, countries_names):

    match_subcountry = difflib.get_close_matches(entry, subcountries_names, 1, 0.9)
    if not match_subcountry:
        match_city = difflib.get_close_matches(entry, cities_names, 1, 0.9)
        if not match_city:
            match_country = difflib.get_close_matches(entry, countries_names, 1, 0.9)
            if not match_country:
                result = ['null', 'null', 'null']
            else:
                match = match_country[0]
                result = countries[match]
        else:
            match = match_city[0]
            result = cities[match]
    else:
        match = match_subcountry[0]
        result = subcountries[match]

    if any(isinstance(i, list) for i in result):
        output = result[0]
    else:
        output = result

    return(output)


#def orf(entry):
    
#    trip_end = ['tga', 'taa' ,'tag', 'TGA', 'TAA' ,'TAG', 'end']
#    trip_start = ['atg', 'gtg, 'ttg']
#    trip_start = ['atg', 'ATG']
#    pos = 1
#    orfstart = 0
#    orfend = 0
#    orflength = 0
#    neworf = 0
#    inframe = 0
#    newgenome = ''
#    frame=''

#    trips = re.findall('...', entry)

#    trips.append('end')

#    for trip in trips:
#        if inframe == 0:
#            if trip in trip_start:
#                orfstart = pos
#                inframe = 1
#                newgenome = trip
#        else:
#            if not trip == 'end':
#                newgenome += trip
#            if trip in trip_end:
#                orfend = pos
#                neworf = orfend - orfstart
#                if neworf > orflength:
#                    orflength = neworf
#                    frame = newgenome
#                inframe = 0
#            pos += 1
#    return(frame)


def worker(cities, subcountries, countries, cities_names, subcountries_names, countries_names, col):
    
    col = list(map(str, col))
    
    head = col[0].split('|')
    
    acc = head[0][1:]
    name = head[1]
    seg = head[2]
    org = head[4]
    sub = head[5]
    if sub == 'NA' or sub == 'nan':
        sub = 'null'
    host = head[7]
    entry = col[1]

    info = name.split('/')

    spec = info[0]
    del info[0]
    year = info[-1]
    del info[-1]

    if year.isdecimal():
        if len(year) == 2:
            year = '19'+year 
    else:
        year = 'null'

    if not info:
        pos = ['null', 'null', 'null']
    else:
        for i in info:
            pos = search(i, cities, subcountries, countries, cities_names, subcountries_names, countries_names)
            if not all([item == 'null' for item in pos]):
                break

    #entry_rev = entry[::-1].translate(str.maketrans('acgtACGT', 'tgcaTGCA'))

    #frame_list = [orf(entry), orf(entry[1:]), orf(entry[2:])]
    #frame_rev_list = [orf(entry_rev), orf(entry_rev[1:]), orf(entry_rev[2:])]

    #frame = max(frame_list, key=len)
    #frame_rev = max(frame_rev_list, key=len)

    #if len(frame) >= len(frame_rev):
    #    if frame == '' or len(frame) <= 30:
    #        frame = entry
    #        utr = ['null', 'null']
    #        genome = entry_rev
    #    else:
    #        utr = entry.split(frame)
    #        if utr[0] == '':
    #            utr[0] = 'null'
    #        if utr[1] == '':
    #            utr[1] = 'null'             
    #        genome = entry_rev
    #else:
    #    if frame_rev == '' or len(frame_rev) <= 30:
    #        frame_rev = entry_rev
    #        utr = ['null', 'null']
    #        genome = entry
    #    else:
    #        utr = entry_rev.split(frame_rev)
    #        if utr[0] == '':
    #            utr[0] = 'null'
    #        if utr[1] == '':
    #            utr[1] = 'null'
    #        genome = entry
    
    #out = [acc, name, seg, utr[0], frame, utr[1], genome + ';' + name, spec, pos[0], pos[1], pos[2], year, sub, host]
    out = [[acc, name, seg, entry], [name, spec, pos[0], pos[1], pos[2], year, sub, host]]
    
    return out


def main():

    worldfile = sys.argv[1]
    reader = pd.read_csv(sys.argv[2], chunksize = 10000, sep = ';', na_filter = False, header = None)
    procs = int(sys.argv[3])
    
    cities = collections.defaultdict(list)
    subcountries = collections.defaultdict(list)
    countries = collections.defaultdict(list)

    with open(worldfile, newline='') as csvfile:
        world = csv.reader(csvfile, delimiter=',')
        header = next(world)
        if header != None:
            for row in world:
                city = row[0]
                country = row[1]
                subcountry = row[2]
                cities[city].append([city, subcountry, country])
                subcountries[subcountry]. append(['null', subcountry, country])
                countries[country].append(['null', 'null', country])

    cities_names = list(cities.keys())
    subcountries_names = list(subcountries.keys())
    countries_names = list(countries.keys())
                
    for df in reader:

        with mp.Pool(procs) as pool:

            result = pool.imap(ft.partial(worker, cities, subcountries, countries, cities_names, subcountries_names, countries_names), df.itertuples(index=False, name=None), chunksize=10)

            with open(sys.argv[4], 'a', newline='') as genomes, open(sys.argv[5], 'a', newline='') as strains:

                for out in result:

                    out_gnms = csv.writer(genomes, delimiter="\t")
                    out_gnms.writerow(out[0])
                    
                    out_strns = csv.writer(strains, delimiter="\t")
                    out_strns.writerow(out[1])


if __name__ == "__main__":

    main()
