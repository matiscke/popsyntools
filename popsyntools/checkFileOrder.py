#!/usr/bin/env python
# coding: utf-8

# ### Check tracks files for consistency
# #### Usage:
# - run "checkPopulation()" from the population folder
# - optional: save output if you are interested in exact positions of errors
# - collect the file 'brokenSims.txt' which contains all erroneous simulations

# In[131]:


import glob
import shutil


# In[84]:


def checkFileOrder(path):
    """ Check if the beginnings of lines in a file are a continuous
    sequence of integers without gaps, starting from 1.

    Returns
    =======
    order : bool
        'True' if everything correct, 'False' otherwise
    """
    with open(path, 'r', errors='replace') as f:
        step = 0
        for line in f.readlines():
            step += 1
            if line.strip().startswith(str(step)):
                continue
            else:
                print('issue at step {}'.format(step))
                return False
    return True


# Run the following from the population folder (the one with the 'SIMXXXXX' folders in it)

# In[90]:


def checkPopulation():
    brokenSimList = set([])
    for path in sorted(glob.glob('**/tracks*.outputdat')):
        if checkFileOrder(path) == False:
            print(path)
            brokenSimList.add(path.split('/')[0])

    with open('brokenSims.txt', 'w') as f:
        for sim in brokenSimList:
            f.write(sim + '\n')
    return


# ### Check broken sims in NGM8:

# In[126]:


def get_brokenSims(brokenSimsFile='brokenSims.txt'):
    with open(brokenSimsFile, 'r') as f:
        brokenSims = sorted([sim.strip() for sim in f.readlines()])
    return brokenSims


# ### Find directories _not_ contained in brokenSims

# In[127]:


def get_goodSims(brokenSims):
    goodSims = sorted([d for d in glob.glob('SIM0*') if not d in brokenSims])
    return goodSims


# ### Delete folders containing bad simulations

# In[136]:


def delete_badSims(brokenSims):
    for sim in brokenSims:
        try:
            shutil.rmtree(sim)
        except FileNotFoundError:
            print('No such file or directory: {}'.format(sim))
            continue
