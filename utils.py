""" This module includes helper functions.

Written by: Martin Schlecker
schlecker@mpia.de
"""

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


def replace_line(filename, pattern, replacement):
    """ Replace a single line in a file.

    Parameters
    ----------
    filename : string
        path to file
    pattern : string
        pattern to search for in the line to replace
    replacement : string
        content of the new line

    Returns
    -------
    """
    from tempfile import mkstemp
    from shutil import move
    from os import fdopen, remove

    #Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(filename) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, replacement))
    #Remove original file
    remove(filename)
    #Move new file
    move(abs_path, filename)