# ---------------------------------------------------------------------- #
# Examples from:
# Ahead of the Count: An Algorithm for Probabilistic Prediction of Instant Runoff (IRV) Elections
# by Kapoor & Staecker
#
# This code is written by Chris Staecker.
# Please contact us with any questions! cstaecker@fairfield.edu
#
# Many of these will dump output into files in the current directory
#
# simpleexample()    - Simple example in the introduction
# house18()          - Example 1 (AK House District 18)
# ussenate()         - Example 2 (AK Senate Murkowski race)
# recountexample()   - Example 3 (simple recount example)
# house15_recount()  - Example 4 (AK House District 15 recount)
# senate_e_recount() - Example 5 (AK Senate District E recount)
# timingtest()       - Timed calculations on random input, 
#                      results quoted in paper conclusion
#
# Also interesting:
# randompt(n) gives a random probability table for n candidates, 
# based on normal distributions (this is used in timingtest)
# ---------------------------------------------------------------------- #

import random
import time

from irvprob import *
from irvsimulator import *
from alaska2022 import *

def simpleexample():
    # SIMPLE EXAMPLE
    d = {"[A,B]":[.01,.17,.34,.25,.13,.1], "[A,C]":[.5,.4,.07,.03], 
         "[B,A]":[0,0,.45,.31,.2,.04], "[B,C]":[0,.1,.3,.27,.19,.14],
         "[C,A]":[.09,.17,.37,.21,.1,.06], "[C,B]":[.17,.75,.08],
         "[A]":[.5,.5], "[B]":[.1,.3,.3,.2,.1],
         "[C]":[.02,.33,.21,.2,.15,.09]}
     
    simplept = ProbabilityTable(d, 100)
    
    print(simplept.dict)
    
    print(simplept.winning_probabilities())

    t = export_elimination_tree(simplept, "simpleexample.tex")
    print(t)

# ALASKA EXAMPLES

aksimpledatafile = "alaska2022-simplified.json"

# cid is the "Contest ID" from the AK data
# fileid is used for output plot graphic filenames
def akplots(cid, fileid,seed=0):
    random.seed(seed)  # to get repeatable RNG

    ae = alaskaElection_from_simple_file(aksimpledatafile)

    # each simulation uses roughly this many buckets in first step
    numbuckets = 100
    
    ps = ae.fullprefschedule(cid)
    
    # full names
    names = {n:ae.candidatenames[int(n)] for n in ps.candidates}
    
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

    wps = irvsimulator(ps,numbuckets, votelist=votelist,stepsize=.005,bucketcap=100,silent=True)
    
    filename = "plot" + fileid + "wps-seed" + str(seed) + ".png"
    plotwps(wps,filename,ae.candidatenames,filename)
    
    printwps(wps,fileid + "wps-seed" + str(seed), names)
    

def ussenate_plots():
    akplots(3,"ussenate-seed4",4)
    akplots(3,"ussenate-seed2",2)
    akplots(3,"ussenate-seed3",3)
    akplots(3,"ussenate-seed0",0)


def house18():
    seed = 4
    #random.seed(45)
    random.seed(seed)
    
    ae = alaskaElection_from_simple_file(aksimpledatafile)
    
    bucketsize = 75
    
    cid = 42 # AK House district 18
    
    fullps = ae.fullprefschedule(cid)
    
    fullps.irv_results()
    
    # arrange votelist sorted by precinct portion
    pps = list(ae.contests[cid].keys())
    random.shuffle(pps)
    
    print(str(len(pps)) + " precinct portions")
    
    ppsizes = []
    votelist = []
    for p in pps:
        
        votedump = []
        for r in fullps.srankings:
            if r in ae.prefschedule(cid,p).dict.keys():
                votedump += [r for k in range(ae.prefschedule(cid,p).dict[r])]
        random.shuffle(votedump)
        ppsizes += [len(votedump)]
        
        votelist += votedump
    
    print(len(votelist),"total votes cast")
    
    # take only half of them
    halfvotelist = votelist[:(int(len(votelist)/2))]
    
    d = {r:0 for r in fullps.srankings}
    for r in halfvotelist:
        d[r] += 1
    
    halfps = PreferenceSchedule(d)
    
    chalfps = halfps.consolidated()

    # full names
    names = {n:ae.candidatenames[int(n)] for n in fullps.candidates}
    
    #{'50': 'Franks, Lyn D.', '51': 'Nelson, David', '52': 'Groh, Cliff'}

    # but we want abbreviated
    names = {'50': "F", '51': "N", '52':"G"}
    
    # the various rankings (consolidated)
    rs = [r for r in chalfps.rankings if chalfps.dict[list_to_stringlist(r)] != 0]
    rs = rs[1:] # remove [] at front of rs
    
    print("\n")
    print(r"\begin{tabular}{r|ccccccccccc}")
    print("ranking: & ",end="")
    for r in rs:
        print(" " + "".join([names[s] for s in r]) + " & ",end="")
    print(r"\\")
    print(r"\hline")
    
    print("votes: & ",end="")
    for r in rs:
        print(" " + str(chalfps.dict[list_to_stringlist(r)]) + " & ",end="")
    print(r"\\")
    print(r"\end{tabular}")
    print("\n")
    
    
    halfpt = ptfrombound(halfps.consolidated(),len(votelist)-len(halfvotelist),bucketsize)
    print(halfpt.winning_probabilities())
    
    print("\n")
    
    print(halfpt.latex(names))
    
    wp50 = halfpt.winning_probabilities()
    
    print("\n")
    print(r"\[")
    for c in names.keys():
        print(r" \qquad " + names[c] + ": " + str(round(100*wp50[c],1)) + r"\% ", end="")
    print("\n" + r"\]")
    print("\n")
    
    export_elimination_tree(halfpt, "akhouse18seed" + str(seed) + ".tex",names)

    random.seed(45)
    wps = irvsimulator(fullps,20,votelist=votelist,stepsize=.005,silent=True)

    # print data points for copy-paste into pgf plots
    printwps(wps,"akhouse18wpsseed" + str(seed), names)
    
    printwpstriangle(wps,"akhouse18wpstriangleseed" + str(seed),names)
    
    title = "AK House District 18"
    
    plotwps(wps,"akhouse18seed" + str(seed) + ".png",ae.candidatenames,"title")
    triangleplotwps(wps,"akhouse18triangleseed" + str(seed) + ".png",ae.candidatenames,"title")

def recountexample():
    d = {"[A,B]":501, "[A,C]":300, 
         "[B,A]":400, "[B,C]":400,
         "[C,A]":200, "[C,B]":600,
         "[A]":500, "[B]":400,
         "[C]":500}
         
    ps = PreferenceSchedule(d)
    print(ps.dict)
    pt = recount_probtable(ps)

    print(pt.latex({'A':'A','B':'B','C':'C'}))
    
    print(pt.winning_probabilities())
    
def house15_recount():
    # Results of this election are very sensitive to slight differences, and the 
    # vote counts in the published AK file are not exactly what were used in the real election
    #  (they apparently handled writins, overvotes, etc slightly differently)
    # so we cannot reproduce exactly the votes which appear in the official results
    # if we use the AK file. 
    
    # Official results are only reported as totals round-by-round, published here:
    #   https://www.elections.alaska.gov/results/22GENR/15.pdf
    # Round 1:
    #  W: 3384
    #  E: 1039
    #  M: 2839
    
    # Round 2:
    #  W: 3476
    #  M: 3483
    
    # We will have to guess some values in the official preference schedule. 
    # We know for sure: 
    #  [E]: 303, [E,W]: 92, [E,M]: 644
    # We know the total of all votes with W on top in round 1 is 3384. Using ratios
    # from my AK file, W should scale by 1.000355, M by 0.99225. 
    # So we guess a preference schedule of:
    # [W]: 2361, [W,E]: 667, [W,M]: 356
    # [M]: 1004, [M,W]: 97, [M,E]: 1738

    
    ps = PreferenceSchedule(
        {'[E]': 303, '[E,W]': 92, '[E,M]': 644,
         '[W]': 2361, '[W,E]': 667, '[W,M]': 356,
         '[M]': 1004, '[M,W]': 97, '[M,E]': 1738})
    pt = recount_probtable(ps)

    print(str(pt.winning_probabilities()))
    
def senate_e_recount():
    # As in House 15, we need to do this by hand.
    
    # Official results are only reported as totals round-by-round, published here:
    #   https://www.elections.alaska.gov/results/22GENR/E.pdf
    
    # Round 1:
    #  G (130): 5652
    #  H (129): 5532
    #  C (131) : 5518
    
    # Round 2:
    #  G: 7881
    #  H: 5949
    
    # We cook up the preference table as above
    # We know: 
    # [C]: 2871, [C,H]: 417, [C,G]: 2229
    # We guess the others
    # and we make [C]: 2872 instead to make the totals work out. 
    # There is a weird "overvote" in the data here
    
    ps = PreferenceSchedule(
        {'[C]': 2872, '[C,H]': 417, '[C,G]': 2229,
         '[H]': 2679, '[H,G]': 2614, '[H,C]': 239,
         '[G]': 1995, '[G,C]': 1444, '[G,H]': 2213 })
    
    pt = recount_probtable(ps)         
    print(str(pt.winning_probabilities()))
    
# generates a random ProbabilityTable for n candidates
# works by generating a random Preference schedule with vote totals random ints in [0,10000],
# then calls ptfrombound on this, imagining that this preference schedule represents
# half of all eventual votes cast.
def randompt(n):
    # ['A','B','C',...] of length n
    candidates = [chr(x) for x in range(ord('A'), ord('Z')+1)][0:n]

    rankings = subsetperms(candidates)
    
    d = {}
    for r in rankings:
        d[list_to_stringlist(r)] = int(random.random() * 10000)
    
    ps = PreferenceSchedule(d)
    votes = sum(d.values())
    
    pt = ptfrombound(ps, votes * 2, 100)
    
    return pt

# tests used for execution times quoted in the paper
# does 5 runs each, prints times in seconds
def timingtest():
    for n in range(3,6):
        pt = randompt(n)
        pt.bucketcap = 1500
        print("Testing " + str(n) + " candidates")
        for i in range(5):
            t0 = time.time()
            pt.winning_probabilities()
            t1 = time.time()
            print(t1-t0)

# ---------------------------------------------------------------------- #
# VARIOUS HELPERS
# ---------------------------------------------------------------------- #

# print data points for copy-paste into pgf plots
def printwps(wps,fileid,names):
    s = ""
    for c in wps.keys():
        s += "Probs for " + names[c] + ":\n"
        l = len(wps[c])
        for i in range(l):
            s += "(" + f'{(i+1)/l:.3f}' + "," + f'{wps[c][i]:.3f}' + ") "
        s += "\n\n"
    
    with open(fileid + ".txt", 'w') as fp:
        fp.write(s)
        fp.close()
        
def printwpstriangle(wps,fileid,names):
    s = ""
    [a,b,c] = list(wps.keys())

    for i in range(len(wps[a])):
        s += f'({wps[a][i]:.3f},{wps[b][i]:.3f},{wps[c][i]:.3f}) '
    s += "\n\n"
    
    with open(fileid + ".txt", 'w') as fp:
        fp.write(str(names)+"\n")
        fp.write(s)
        fp.close()