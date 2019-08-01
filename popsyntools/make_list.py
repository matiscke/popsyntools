"""Simulation list manipulations.

The functions in this script are meant to create or manipulate a list of
initial conditions for the Planet Population Synthesis
Code 'Planete' by the Bern planet formation group.

Written by: Martin Schlecker
schlecker@mpia.de
"""
import pandas as pd
import numpy as np
import warnings

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

    try:
        # old simlist format (roughly until summer 2019)
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
                "sigma0",
                "a_in",
                "a_out",
                "expo",
                "mWind",
                "simName",
                "a_start",
                "t_start"]
        simlist = pd.DataFrame(simlist, columns=columns)

        # set the correct data types
        simlist[['fgp','sigma0','a_in','a_out','expo','mWind',
                 'a_start','t_start']] = simlist[['fgp','sigma0','a_in','a_out',
                 'expo','mWind','a_start','t_start']].apply(pd.to_numeric)

    except AssertionError:
        # new format (from early "NG" populations on (summer 2019))
        simlist = []
        with open(filename) as f:
            for line in f:
                line_split = line.split()
                line = [item.split('=')[1] for item in line_split[1:]]
                simlist.append(line)

        # drop last row (contains only "END")
        simlist = simlist[:-1]

        # turn list into a pandas DataFrame
        columns = [
                "Mstar",
                "sigma0",
                "expo",
                "a_in",
                "a_out",
                "fgp",
                "mWind",]
        simlist = pd.DataFrame(simlist, columns=columns, dtype=np.float64)

    return simlist


def write_header(filename, colname, func, *funcArgs, **funcKwargs):
    import os
    filename += '.hdr'
    # if filename in os.listdir():
    headertext = 'Column "{}" filled with values from function {}'.format(colname,
        str(func))


    with open(filename) as f:
        f.write(headertext)


def grid(start, stop, size, nsamples=50, log=False):
    """ Return a one-dimensional grid of repeating, evenly distributed values.

    The values in the grid are repeated until the specified length of the array
    "size" is reached.

    Parameters
    ----------
    start : scalar
        The starting value of the sequence.
    stop : scalar
        The end value of the sequence.
    size : int
        length of resulting array
    nsamples : int, optional
        Number of samples to generate. Default is 50. Must be non-negative.
    log : boolean, optional
        if True, samples are spaced evenly on a log scale.

    Returns
    -------
    grid : numpy array
        one-dimensional array of length `size` consisting of repeating sequences
        of `nsamples` equally spaced samples.

    Examples
    --------
    >>> grid(1, 10, 10, nsamples=4)
    array([ 1.,  4.,  7., 10.,  1.,  4.,  7., 10.,  1.,  4.])
    """
    if nsamples > size:
        warnings.warn("number of samples > length of data column.")

    if log:
        grid = np.tile(np.logspace(np.log10(start), np.log10(stop), nsamples),
            int(np.ceil(size/nsamples)))
    else:
        grid = np.tile(np.linspace(start, stop, nsamples),
            int(np.ceil(size/nsamples)))
    return grid[:size]


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
    if not colname in simlist.columns:
        warnings.warn('Table has no column "{}". A new column is created.')

    try:
        # vectorial method for functions accepting a "size" argument
        simlist[colname] = func(*funcArgs, size=len(simlist), **funcKwargs)
    except:
        # or a "num" argument
        try:
            simlist[colname] = func(*funcArgs, num=len(simlist), **funcKwargs)
        except:
            # list comprehension for functions returning single values
            simlist[colname] = np.array([func(*funcArgs, **funcKwargs)
                                        for i in range(len(simlist))])

    # # write reference to header
    # write_header(filename, colname, func, *funcArgs, **funcKwargs)
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

    line = 'CD_' + var2str(singleSim['CDnumber'])\
    + 'FP_' + var2str(singleSim['fgp'])\
    + 'SI_' + var2str(singleSim['sigma0'])\
    + 'AI_' + var2str(singleSim['a_in'])\
    + 'AO_' + var2str(singleSim['a_out'])\
    + 'EX_' + var2str(singleSim['expo'])\
    + 'MW_' + var2str(singleSim['mWind'])\
    + 'SIM' + var2str(singleSim['simName'])\
    + 'AS_' + var2str(singleSim['a_start'])\
    + 'ST_' + var2str(singleSim['t_start'])
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
    with open(filename, 'w') as fileHandle:
        for i in range(len(simlist)):
            row = simlist.iloc[i]
            write_singleSim2File(fileHandle, row)
        fileHandle.write('END')
    print('{} simulations written to "{}".'.format(i, filename))
