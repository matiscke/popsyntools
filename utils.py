""" This module includes helper functions.

Written by: Martin Schlecker
schlecker@mpia.de
"""
import numpy as np
import pandas as pd

def get_M0(rc, Sigma0, expo, r0=5.2):
    """Compute the total disk mass from initial condition parameters.

    Parameters
    ----------
    rc : float
        characteristic radius [au]
    Sigma0 : float
        gas surface density at 5.2 au [g/cm^2]
    expo : float
        Power law slope
    r0 : float
        reference radius [au], in general 5.2 au

    Returns
    -------
    M0 : float
        total disk mass in solar masses
    """
    M0 = (2*np.pi)/(2-expo)*Sigma0*(r0*au)**expo*(rc*au)**(2-expo)
    return M0/Msol


def get_Sigma0(rc, M0, expo, r0=5.2):
    """Compute Sigma0, the gas surface density at the reference radius r0
    necessary for an initial total disk mass of M0.

    Parameters
    ----------
    rc : float
        characteristic radius [au]
    M0 : float
        total disk mass in solar masses
    expo : float
        Power law slope
    r0 : float
        reference radius [au], in general 5.2 au

    Returns
    -------
    Sigma0 : float
        gas surface density at 5.2 au [g/cm^2]
    """
    return (r0*au)**(-expo)*(rc*au)**(-2+expo)*(2-expo)*M0*Msol/(2*np.pi)


def get_orbitalPeriod(population, MstarRel=0.1):
    """ Compute the orbital period P from the semi-major axis a and Mstar.

    get_orbitalPeriod uses Kepler's Third Law to calculate the orbital period of
    planets from their semi-major axis and a given stellar mass Mstar. It adds a
    new column 'period' to the population table.

    Entries with negative semi-major axes are removed.

    The semi-major axis must be given in [au], period will be in [d].

    Parameters
    ----------
    population : Pandas DataFrame
        Table with the population. Has to contain a column 'a' (semi-major axis)
    MstarRel : float
        Mass of the stellar host in solar Masses

    Returns
    -------
    pop_posSma : Pandas DataFrame
        Table with additional column `period` and entries with negative semi-
        major axis removed

    Example
    -------
    >>> MstarRel = 1.0
    >>> a_Earth = 1.0
    >>> a_Mars = 1.523662
    >>> test = pd.DataFrame({'a' : [a_Earth, a_Mars]})
    >>> get_orbitalPeriod(test, MstarRel)
            a      period
    0  1.000000  365.257762
    1  1.523662  686.961516
    """
    # convert a from au to cm
    sma_cm = lambda sma_au : sma_au*au

    # Remove entries with negative semi-major axis
    pop_posSma = population[population['a'] > 0.].copy()

    Mstar = MstarRel*Msol
    KeplerConst = 4*np.pi**2/(G*Mstar)
    pop_posSma['period'] = np.sqrt(KeplerConst*sma_cm(pop_posSma['a'])**3)

    # convert period from seconds to days
    pop_posSma['period'] = pop_posSma['period']/86400
    return pop_posSma


def replace_line(filename, pattern, replacement, newFile=None, backup=True):
    """ Replace a single line in a file.

    Parameters
    ----------
    filename : string
        path to file
    pattern : string
        pattern to replace in the file
    replacement : string
        content of the new line
    newFile : string, optional
        filename of changed file
    backup : bool, optional
        create a backup file before overwriting. The backup file has the name
        of the original file with '.bak' added. Not active if newFile is given.
    """
    from tempfile import mkstemp
    from shutil import move, copy2
    from os import fdopen, remove

    # Create temp file
    fh, abs_path = mkstemp()

    with fdopen(fh,'w') as new_file:
        with open(filename) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, replacement))

    if not newFile:
        newFile = filename
        if backup:
            copy2(filename, filename + '.bak')
        remove(filename)

    move(abs_path, newFile)


def linearScale(x1, x2, y1, y2, x):
    """ Evaluate y(x) of a linear function going through (x1,y1) and (x2,y2).

    The function returns y = ((x-x2)/(x1-x2))*y1+((x-x1)/(x2-x1))*y2

    Parameters
    ----------
    x1 : float
        first x value
    x2 : float
        second x value
    y1 : float
        first y value
    y2 : float
        second y value
    x : float
        x value at which to evaluate the function

    Returns
    -------
    y : float
        y(x) of above function
    """
    return ((x-x2)/(x1-x2))*y1+((x-x1)/(x2-x1))*y2


def get_diskFractions(population, timeColumn='t', Ntimes=100):
    """ Compute fractions of remaining disks from disk dispersal times.

    For a list of disk dispersal times, the function evaluates for a grid of
    times between zero and the maximum dispersal time how many disks are left.

    Parameters
    ----------
    population : Pandas DataFrame
        a data frame with a column containing the disk dispersal times
    timeColumn : string
        name of the column containing the dispersal times
    Ntimes : integer
        number of times at which to evaluate the disk fraction

    Returns
    -------
    times : numpy array
        times at which the disk fractions were evaluated
    fractions : numpy array
        fractions of remaining disks for each time

    Example
    -------
    >>> disks = pd.DataFrame([3799270., 3495080., 7174510., 8329200., 4125000.],columns=['t'])
    >>> get_diskFractions(disks, Ntimes=5)
    (array([      0., 2082300., 4164600., 6246900., 8329200.]),
     array([1. , 1. , 0.4, 0.4, 0. ]))
    """
    nDisks = len(population)
    times = np.linspace(0., max(population[timeColumn]), Ntimes)
    return times, np.array([len(population[population[timeColumn] > t])/nDisks
        for t in times])
