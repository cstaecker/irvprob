# irvprob
Probabilistic predictions of IRV elections

Functions for probabilistic analysis of IRV elections, implementing the algorithm from: *Ahead of the Count: An Algorithm for Probabilistic Prediction of Instant Runoff (IRV) Elections* by Nick Kapoor & Chris Staecker

Read the paper at [arxiv:2405.09009](https://arxiv.org/abs/2405.09009)

This code is written by Chris Staecker. Please contact us with any questions!

Includes files:
 * irvprob.py (main functions here)
 * irvsimulator.py (for simulating realtime election-night predictions)
 * examples.py (examples used in the paper)
 * alaska2022.py (required to analyze examples from Alaska 2022 elections)
 * alaska2022-simplified.json (JSON containing preference schedules for all AK 2022 election races)
 * alaska2022-candidates.txt (Big text file of all candidate names and their internal candidate ID#)
 * alaska2022-contests.txt   (Big text file of all individual contests and their internal contest ID#)

For a basic example, try:
```
from irvprob import *
from examples import randompt

pt = randompt(4)                  # generate a random probability table with 4 candidates
print(pt.winning_probabilities()) # calculate probabilities of wins for each one
```

If all you care about is looking at Alaska vote data, try:
```
from irvprob import *
from alaska2022 import *

ae = alaskaElection_from_simple_file("alaska2022-simplified.json")

# gets preference schedule for AK contest #42, which is AK House District 18
ps = ae.fullprefschedule(42) 

print(ps.dict)   # show vote totals for various rankings

ps.irv_results() # find IRV result for this election

