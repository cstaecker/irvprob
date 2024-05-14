# ---------------------------------------------------------------------- #
# Functions for parsing data from the AK 2022 general election using RCV 
# Used in the paper:
# Ahead of the Count: An Algorithm for Probabilistic Prediction of Instant Runoff (IRV) Elections
# by Kapoor & Staecker
#
# Creates PreferenceSchedule objects, imported from irvprob
#
# This code is written by Chris Staecker.
# Please contact us with any questions! cstaecker@fairfield.edu
# ---------------------------------------------------------------------- #

import json
import random

from sys import path
path.append("../imports")

from irvprob import *
from irvsimulator import *

# ---------------------------------------------------------------------- #
# AlaskaElections object
# 
# properties: 
#  contestnames: dict with keys contestIds (numbers) and string values. 
#                 These are the individual races
#  candidateneames: dict with keys candidateIds (numbers) and string values.
#  contests: dict with keys contest ids, values a dict with keys precinctportions, values PreferenceSchedule objects
# 
# methods:
# fullprefschedule(cid) - gives a full pref schedule for a certain contest id
# prefschedule(cid,pp) - gives a pref schedule for a certain contest id and precinct portion
# 
# For example to get a preference schedule for AK Senate District E (Holland/Giessel/Cacy), 
# this is contest ID 10, so we do:
# 
#  ae = alaskaElection_from_simple_file("alaska2022-simplified.json")
#  ps = ae.fullprefschedule(10)
#  print(ps.dict)
# 
# The candidates are given just as numbers 129, 130, 131. To see the real names, do like:
#  ae.candidatenames[129]
#
# The method above uses our simplified json summary file. If you want to 
# use the AK data directly from them, see function build_simple_file below
#
# Our simplified json file makes choices about weird edge cases (overvotes, writins, etc)
# which result in slightly different vote tallies than the "official" numbers. 
# If you need super-accurate counts, look very closely at the details of 
# build_simple_file. (or just do it yourself!)
# ---------------------------------------------------------------------- #

class AlaskaElections:
    def __init__(self,d):
        self.contests = d["contests"]
        
        cids = list(d["contests"].keys())
        cids.sort()
        self.contestIds = cids
        
        self.contestnames = d["contestnames"]
        self.candidatenames = d["candidatenames"]
        
    def prefschedule(self,cid,pp):
        return self.contests[cid][pp]
        
    def fullprefschedule(self,cid):
        pps = list(self.contests[cid].keys())
        p = self.prefschedule(cid,pps[0])
        others = [self.prefschedule(cid,qs) for qs in pps[1:]]
        for q in others:
            p = prefScheduleSum(p,q)
        return p
            
def has_duplicates(l):
    return len(l) != len(set(l))

#-----------------------------------------------------------------------------
# Builds the AlaskaElections object from the "alaska2022-simplified.json" file
#-----------------------------------------------------------------------------
def alaskaElection_from_simple_file(filename):
    print("Loading election data from " + filename)
    with open(filename) as json_file:
        data = json.load(json_file)
    
    d = {}
    d["contests"] = {int(c):{int(p):PreferenceSchedule(data["contests"][c][p]) for p in data["contests"][c].keys()} for c in data["contests"].keys()}
    d["contestnames"] = {int(id):data["contestnames"][id] for id in data["contestnames"].keys()}
    d["candidatenames"] = {int(id):data["candidatenames"][id] for id in data["candidatenames"].keys()}
    
    return(AlaskaElections(d))


#---------------------------------------------------------------------
# Original code used to create the "alaska2022-simplified.json" file
# uses the CVR_Export file (1.5GB) obtained from Alaska Board of Elections
# (click link for "Cast Vote Record")
#---------------------------------------------------------------------

# Reads a list of raw CVR JSON files into an AlaskaElection object
def alaskaElection_from_rawfiles(jsons):
    votes = []
    
    for jsonfilename in jsons:
        with open(jsonfilename) as json_file:
            data = json.load(json_file)
        print("Tallying " + jsonfilename)
        
        for s in data["Sessions"]:
            pp = s["Original"]["PrecinctPortionId"]
        
            for card in s["Original"]["Cards"]:
                for contest in card["Contests"]:
                    v = {}
                    
                    v["contestId"] = contest["Id"]
                    v["precinctPortion"] = pp
                    #cid = contest["Id"]
                    
                    ranking = []
                    for m in contest["Marks"]:
                        # ignore writein candidates
                        if "WriteinIndex" not in m.keys():
                            ranking.append(str(m['CandidateId']))
                    
                    # some people have repeats in their rankings- throw these out
                    if not has_duplicates(ranking):
                        v["ranking"] = ranking
                        votes.append(v)
    
    # make these votes into a dict of dicts of dicts
    d = {}
    for v in votes:
        cid = v["contestId"]
        pp = v["precinctPortion"]
        sr = list_to_stringlist(v["ranking"])
        if cid not in d.keys():
            d[cid] = {}
        if pp not in d[cid].keys():
            d[cid][pp] = {}
        if sr not in d[cid][pp].keys():
            d[cid][pp][sr] = 0
            
        d[cid][pp][sr] += 1
        
    # dict to be fed into AlaskaElections constructor
    dd = {}
    dd["contests"] = {}
    for cid in d.keys():
        dd["contests"][cid] = {}
        for pp in d[cid].keys():
            dd["contests"][cid][pp] = PreferenceSchedule(d[cid][pp])
    
    # retrieve contest IDs
    dd["contestnames"] = {}
    with open("CVR_Export/ContestManifest.json") as json_file:
        data = json.load(json_file)
        for d in data["List"]:
            dd["contestnames"][d["Id"]] = d["Description"]
    
    dd["candidatenames"] = {}
    with open("CVR_Export/CandidateManifest.json") as json_file:
        data = json.load(json_file)
        for d in data["List"]:
            dd["candidatenames"][d["Id"]] = d["Description"]
    ae = AlaskaElections(dd)
    return ae
    
# read all raw CVR files and build the alaska2022-simplified.json file.
def build_simple_file(filename):
    print("Building simple file " + filename)
    jsons = ["CVR_Export/CvrExport_" + str(n) + ".json" for n in range(1,1887)]
    #jsons = ["CVR_Export/CvrExport_" + str(n) + ".json" for n in range(1,18)]
    
    
    ae = alaskaElection_from_rawfiles(jsons)
    
    d = {}
    d["contestnames"] = ae.contestnames
    d["candidatenames"] = ae.candidatenames
    d["contests"] = { c: {pp: ae.contests[c][pp].dict for pp in ae.contests[c].keys()} for c in ae.contestIds}

    with open(filename, 'w') as fp:
        json.dump(d,fp)


#---------------------------------------------------------------------
# For finding interesting examples to be used in the paper
# runs irvsimulator on every contest in the dataset with 3 or 4 candidates, 
# does each one with 5 different set randomization seeds, 
# builds plots for all of them
#
# Plots go in directory called "plots", which should exist
#---------------------------------------------------------------------

def big_alaska_simulation():
    simpledatafile = "alaska2022-simplified.json"
    ae = alaskaElection_from_simple_file(simpledatafile)
    
    # each simulation uses roughly this many buckets
    numbuckets = 100
    
    numseeds = 5
    
    contests = ae.contestIds

    for cid in contests:
        ps = ae.fullprefschedule(cid)
        
        if len(ps.candidates) > 2 and ps.votecount > 100:
            for seed in range(numseeds):
                print("Simulating CID:" + str(cid) + " Seed " + str(seed) + " (" + str(len(ps.candidates)) + " candidates)")
                random.seed(seed)
                
                # arrange votelist sorted by precinct portion
                pps = list(ae.contests[cid].keys())
                random.shuffle(pps)
                
                ppsizes = []
                votelist = []
                for p in pps:
                    
                    votedump = []
                    for r in ps.srankings:
                        if r in ae.prefschedule(cid,p).dict.keys():
                            votedump += [r for k in range(ae.prefschedule(cid,p).dict[r])]
                    random.shuffle(votedump)
                    ppsizes += [len(votedump)]
                    
                    votelist += votedump
                                    
                wps = irvsimulator(ps,numbuckets, votelist=votelist)
                
                id = "cid" + str(cid) + "-seed" + str(seed)
                filename = "plots/" + id + ".png"
                
                title = ae.contestnames[cid] + "\n" + str(ps.votecount) + " votes"
                plotwps(wps,filename,ae.candidatenames,title)

                if len(ps.candidates) == 3:
                    tfilename = "plots/" + id + "-triangle.png"
                    triangleplotwps(wps,tfilename,ae.candidatenames,title)
