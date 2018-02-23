"""Simulation list manipulations.

The functions in this script are meant to create or manipulate a list of
initial conditions for the Planet Population Synthesis
Code 'Planete' by the Bern planet formation group.

The functions are meant to replace the ones in the "make_liste.f" FORTRAN Code
that was being used to create the same kind of lists.

Written by: Martin Schlecker
schlecker@mpia.de
"""
import pandas as pd



def read_simlist(filename):
    """Read a simulation list from a file.

    Parameters
    ----------
    filename : string
        path to the simulation list

    Returns
    -------
    simlist : Pandas dataframe
        simulation list as a dataframe
    """
    from io import StringIO

    columns = {
            "CDname": "str",
            "CDnumber": "float64",
            "fgp": "float64",
            "diskM": "float64",
            "a_in": "float64",
            "a_out": "float64",
            "expo": "float64",
            "windM": "float64",
            "simName": "str",
            "a_start": "float64",
            "t_start": "float64"}

    # dirty hack because of inconsistent usage of separators in file
    simlist = pd.read_csv(StringIO(''.join(l.replace('SIM', '_SIM')
                          for l in open(filename))), sep='_', dtype='str',
                          names=columns.keys())

    # drop last row
    simlist.drop(simlist.tail(1).index, inplace=True)

    # decontaminate numeric values from non-numeric characters
    simlist["CDnumber"] = simlist["CDnumber"].str.extract('(\d+)', expand=False)
    simlist["fgp"] = simlist["fgp"].str[:-2]
    simlist["diskM"] = simlist["diskM"].str[:-2]
    simlist["a_in"] = simlist["a_in"].str[:-2]
    simlist["a_out"] = simlist["a_out"].str[:-2]
    simlist["expo"] = simlist["expo"].str[:-2]
    simlist["a_start"] = simlist["a_start"].str[:-2]

    # finally we can set the correct data types
    simlist = simlist.astype(columns)
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
    singleSim : pandas DataFrame line
        one line of a simulation list
    """
    line = 'CD_' + singleSim['CDnumber'] + 'FP_' + singleSim['fgp'] +
           'SI_' + singleSim['diskM'] + 'AI_' + singleSim['a_in'] +
           'AO_' + singleSim['a_out'] + 'EX_' + singleSim['expo'] +
           'MW_' + singleSim['windM'] + 'SIM' + singleSim['simName'] +
           'AS_' + singleSim['a_start'] + 'ST_' + singleSim['t_start']
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
