"""Simulation list manipulations.

The functions in this script are meant to create or manipulate a list of
initial conditions for the Planet Population Synthesis
Code 'Planete' by the Bern planet formation group.

Written by: Martin Schlecker
schlecker@mpia.de
"""
import pandas as pd
import numpy as np

def read_simlist(filename, varlen=17):
    """Read a simulation list from a file.

    Parameters
    ----------
    filename : string
        path to the simulation list
    varlen : int
        Number of characters for each variable

    Returns
    -------
    simlist : Pandas dataframe
        simulation list as a dataframe
    """

    simlist = []
    with open(filename) as f:
        # get line length without escape characters and reset iterator
        lineLen = len(f.readline().rstrip('\n'))
        f.seek(0)
        for line in f:
            """read line by line with new parameter every varlen characters.
            In each parameter, the first 3 characters are omitted (they
            contain parameter designations such as "EX_").
            """
            simParams = [line[i+3:i+varlen] for i in range(0, lineLen, varlen)]
            simlist.append(simParams)

    # drop last row (contains only "END")
    simlist = simlist[:-1]

    # turn list into a pandas DataFrame
    columns = [
            "CDnumber",
            "fgp",
            "diskM",
            "a_in",
            "a_out",
            "expo",
            "windM",
            "simName",
            "a_start",
            "t_start"]
    simlist = pd.DataFrame(simlist, columns=columns)

    # set the correct data types
    simlist[['fgp','diskM','a_in','a_out','expo','windM',
             'a_start','t_start']] = simlist[['fgp','diskM','a_in','a_out',
             'expo','windM','a_start','t_start']].apply(pd.to_numeric)
    return simlist


def changeListCol(simlist, colname, func, *funcArgs, **funcKwargs):
    """Fill a single column in a simulation list with new values.

    Parameters
    ----------
    simlist : Pandas dataframe
        simulation list
    colname : string
        name of the column to change
    func : function handle
        function to call that returns values for the column
    *funcArgs : arbitrary
        positional arguments for func
    **funcKwargs : arbitrary
        keyword arguments for func

    Returns
    -------
    simlist : pandas dataframe
        edited simulation list
    """
    try:
        # vectorial method for functions accepting a "size" argument
        simlist[colname] = func(*funcArgs, size=len(simlist), **funcKwargs)
    except:
        # list comprehension for functions returning single values
        import numpy as np
        simlist[colname] = np.array([func(*funcArgs, **funcKwargs)
                                    for i in range(len(simlist))])
    return simlist


def write_singleSim2File(fileHandle, singleSim):
    """Writes a single simulation (a "line" of a DataFrame) to a file.

    Parameters
    ----------
    fileHandle : file object
        handle to an opened output file
    singleSim : pandas DataFrame row object
        one line of a simulation list
    """
    def var2str(variable):
        """helper function to distinguish between floats and other data types.
        Returns a string of length 11.
        """
        try:
            # if variable is a float, return in scientific format
            return '{:1.8E}'.format(variable)
        except ValueError:
            # if variable is a string/object, return string of length 14
            return '{:11s}'.format(variable)

    line = 'CD_' + var2str(singleSim['CDnumber']) + 'FP_' + var2str(singleSim['fgp']) + 'SI_' + var2str(singleSim['diskM'])+ 'AI_' + var2str(singleSim['a_in'])+ 'AO_' + var2str(singleSim['a_out'])+ 'EX_' + var2str(singleSim['expo'])+ 'MW_' + var2str(singleSim['windM'])+ 'SIM' + var2str(singleSim['simName'])+ 'AS_' + var2str(singleSim['a_start'])+ 'ST_' + var2str(singleSim['t_start'])
    line += '\n'
    fileHandle.write(line)


def write_simlist(filename, simlist):
    """Write a simulation list to a file.

    Parameters
    ----------
    filename : string
        path to the output file
    simlist : Pandas dataframe
        simulation list as a dataframe
    """
    simlist.to_csv(filename)
