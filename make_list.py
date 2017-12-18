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

    simlist = pd.DataFrame()
    with open(filename, 'r') as listfile:
        for i,line in enumerate(listfile):
            # format consistent with driver.py by Alexandre Emsenhuber

            print("line {}".format(i))

            data = {
                    "CDname": line[0:17],
                    "CDnumber": line[3:17],
                    "fgp": line[20:34],
                    "diskM": line[37:51],
                    "a_in": line[54:68],
                    "a_out": line[71:85],
                    "expo": line[88:102],
                    "windM": line[105:119],
                    "simName": line[119:136],
                    "simNumber": line[122:136],
                    "a_start": line[139:153],
                    "t_start": line[156:170]}
            simlist = simlist.append(data, ignore_index=True)
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
        function to call that returns a single value for the column
    *funcArgs : arbitrary
        positional arguments for func
    **funcKwargs : arbitrary
        keyword arguments for func

    Returns
    -------
    simlist : pandas dataframe
        edited simulation list
    """

    for val in simlist[colname]:
        val = func(*funcArgs, **funcKwargs)

    return simlist
