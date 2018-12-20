""" This module includes helper functions.

Written by: Martin Schlecker
schlecker@mpia.de
"""
import numpy as np
import pandas as pd

# Define au, Msol, Gravitational Constant in cm, g, s
au = 1.496e13
Msol = 1.98855e33
G = 6.6740831e-8

def get_M0(Sigma0, rc, expo, r0=5.2):
    """Compute the total gas disk mass from initial condition parameters.

    Parameters
    ----------
    Sigma0 : float
        gas surface density at 5.2 au [g/cm^2]
    rc : float
        characteristic radius [au]
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


def get_Sigma0(M0, rc, expo, r0=5.2):
    """Compute Sigma0, the gas surface density at the reference radius r0
    necessary for an initial total disk mass of M0.

    Parameters
    ----------
    M0 : float
        total disk mass in solar masses
    rc : float
        characteristic radius [au]
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

def get_Msolid(Sigma0, rc, expo, dgr, r0=5.2, condensationFactor=0.65):
    """ Compute solid disk mass in Earth masses using disk initial conditions.

    Parameters
    ----------
    Sigma0 : float
        gas surface density at 5.2 au [g/cm^2]
    rc : float
        characteristic radius [au]
    expo : float
        Power law slope
    dgr : float
        dust-to-gas ratio.
    r0 : float, optional
        reference radius [au], in general 5.2 au
    condensationFactor : float, optional
        compensates for solid material inside the ice line which is not
        completely condensed. The default value was derived by comparison with
        the planetesimal disk mass 'mdiskp' at early simulation times.
    Returns
    -------
    Msolid : float
        total solid disk mass in Earth masses
    """
    M0 = get_M0(Sigma0, rc, expo, r0=5.2)

    # correct for not fully condensed material within ice line
    M0 *= condensationFactor

    # multiply by gas-to-dust ratio and transfer from Solar masses to Earth masses
    return M0*332948.6*dgr


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


def linear(x, a, b):
    """Evaluate a two-parameter linear function."""
    y = a*np.array(x) + b
    return y


def exponential(t, A, tau, C):
    """Evaluate a three-parameter exponential function."""
    return A*np.exp(-t/tau) + C


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
    >>> disks = pd.DataFrame([3.8e6, 3.5e6, 7e6, 9e5, 1.1e6],columns=['t'])
    >>> get_diskFractions(disks, Ntimes=5)
    (array([      0., 1750000., 3500000., 5250000., 7000000.]),
    array([1. , 0.6, 0.4, 0.2, 0. ]))
    """
    nDisks = len(population)
    times = np.linspace(0., max(population[timeColumn]), Ntimes)
    return times, np.array([len(population[population[timeColumn] > t])/nDisks
        for t in times])


def fit_diskFractions(times, fractions, func=exponential,
                      paramsBounds=([0.,1e3,-np.inf], [10.,1e9, np.inf])):
    """ fit a three-parameter function to disk fractions.

    Uses the non-linear least squares fit scipy.optimize.curve_fit.

    Parameters
    ----------
    times : numpy array
        times at which the disk fractions were evaluated
    fractions : numpy array
        fractions of remaining disks for each time
    func : function handle
        function to use for the fit
    paramsBounds : tuple
        list of lower (first element) and upper (second element) limits for
        fit parameters

    Returns
    -------
    params : array
        list with the best-fit parameters
    covariance : array
        covariance matrix of the fit
    """
    from scipy.optimize import curve_fit

    params, covariance = curve_fit(exponential, times, fractions,
                                   bounds=paramsBounds)
    params_std = np.sqrt(np.diag(covariance))
    print(params)
    print('std of parameters: {}'.format(params_std))

    return params, covariance


def get_Trappist1data():
    """
    return a table with the planet parameters of the TRAPPIST-1 system.

    Masses, semimajor axes, and radii are from Grimm et al. 2018.
    Ice mass fractions are computed by the internal structure model of the
    COMPLETO code of Bern (2018).

    Returns
    -------
    df : pandas DataFrame
        Data frame containing the planet parameters, indexed by planet names
    """
    T1data = {'name': ['b', 'c', 'd', 'e', 'f', 'g', 'h'],
              'a': [0.01155, 0.01582, 0.02228, 0.02928, 0.03853, 0.04688, 0.06193],
             'r': [1.12, 1.10, .766, .913, 1.05, 1.15, .775],
             'm': [1.017,1.156,0.297, 0.772,0.934,1.148,0.331],
             'ice': [0.26,0.09,0.29,0.0,0.11,0.22,0.19]}
    return pd.DataFrame(T1data, index=T1data['name'])

def convert_dgr2metallicity(population):
    """ convert dust/gas ratio into metallicity, add it to the population

    Assumes a metallicity distribution with a mean of 0.02 (see Mordasini et al.
    (2009). Extrasolar planet population synthesis II. Statistical comparison
    with observations, 1184, 1161–1184.)

    Parameters
    ----------
    population : Pandas DataFrame
        Table with the population. Has to contain a column 'dust2gas'.

    Returns
    -------
    population: Pandas DataFrame
        Table with additional column `metallicity'.
    """
    population['metallicity'] = np.log10(population['dust2gas']/0.02)
    return population


def get_RVsemiamplitude(Mp, P, e=0., Mstar=1.0):
    """ Compute the radial velocity semi-amplitude K.

    Uses the formula in Cumming et al. 1999 (Eqn. 1).

    Parameters
    ----------
    Mp : scalar or array
        planet mass in Mearth
    P : scalar or array
        orbital period in d
    e : scalar or array
        planet eccentricity
    Mstar : scalar
        mass of the host star in Msol

    Returns
    -------
    K : scalar or array
        RV semi-amplitude in m/s
    """
    # some conversions (into m/kg/s SI units)
    Mearth = 5.9722e24
    Msol = 1.9891e30
    d = 86400
    G = 6.67e-11

    return (2*np.pi*G/(P*d))**(1/3)*Mp*Mearth/((Mstar*Msol + Mp*Mearth)**(2/3))*(1-e**2)**(-1/2)


def r_Jup2r_Earth(r):
    """ Transform a radius given in Jupiter radii into Earth radii."""
    return r*10.973


def get_msini(m):
    """ transform a given true m into m*sin(i) with random isotropic inclination i.

    Assuming isotropic orientation of orbits, the distribution of inclination
    angles i follows sin(i) (where this is already normalized for 0 < i < pi/2).
    In order to easily draw x from a uniform distribution in [0,1), we transform
    the distribution to i = arccos(1 - x)

    It is not very reasonable to draw the inclination randomly for individual
    planets in the same system.

    Parameters
    ----------
    m : scalar or array
        planet mass

    Returns
    -------
    msini : scalar or array
        mass times a random (isotropic) sin(i)
    """
    sini = np.sin(np.arccos(1 - np.random.random(len(m))))
    return m*sini


def get_sini(dummy=None):
    """ same as 'get_msini', but returns a single sin(i) ignoring the mass
    """
    return np.sin(np.arccos(1 - np.random.random()))


def add_msini(pop):
    """ add a new column 'msini' from planet masses and random isotropic i.

    This function goes through a population system by system and draws the same
    random inclination for all planets in a given system. The sine of it is then
    multiplied by each planet's mass and written to a column 'msini'.

    Parameters
    ----------
    pop : pandas DataFrame
        planet population

    Returns
    -------
    pop : pandas DataFrame
        planet population including a new column 'msini'
    """
    pop['msini'] = np.nan
    for i in pop.isystem:
        pop.loc[pop.isystem == i, 'msini'] = get_sini()*pop[pop.isystem == i]['m']
    return pop


def get_label(pType):
    """ Translate planet types to nicely formatted labels (e.g. for figure legends)"""
    dictionary = {
    'SuperEarth' : 'Super Earth',
    'SuperEarth_ejected' : 'Ejected Super Earth',
    'ColdJupiter' : 'Cold Jupiter',
    'ColdJupiter_ejected' : 'Ejected Cold Jupiter',
    'HotJupiter' : 'Hot Jupiter',
    'HotJupiter_ejected' : 'Ejected Hot Jupiter',
    'WarmJupiter' : 'Warm Jupiter',
    'WarmJupiter_ejected' : 'Ejected Warm Jupiter',
    'all' : 'All',
    'allTypes' : 'All'}
    return dictionary[pType]


def get_semiMinorAxis(a, e):
    """ Compute semi-minor axis from semi-major axis a and eccentricity e."""
    return a*np.sqrt(1 - e**2)


def get_ApoPeri(a, e):
    """ Compute orbital radius at apoapsis and periapsis from semi-major axis a
    and eccentricity e.
    """
    r_per = (1-e)*a
    r_apo = (1+e)*a
    return r_per, r_apo
