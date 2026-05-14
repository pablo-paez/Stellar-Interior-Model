from .Optimization_functions import refinement
from .simple_model import model

def optimized_model(Mtot: float,Ltot: float,Tc: float,Rtot: float,X: float,Y: float,n_capas: int= 100,Rini=None,tol=0.0001,tol2=1e-5,g=5/3, r_inf = 11.5,r_sup = 12.5,L_inf = 35,L_sup = 50,itL=11,itr=11,it: int=3, debug: bool =True, debug2: bool =True, debug_opt_param: bool =True):
    """
    Function that calculates the complete model, including optimise the radius, brightness and temperature settings.

    Parameters
    ----------
    Mtot: float
        Total mass [10^33 g].
    Ltot: float
        Total brightness [10^33 erg/s].
    Tc: float
        Central temperature [10^7 K].
    Rtot: float
        Total radio [10^10 cm].
    X : float
        Mass fraction of hydrogen.
    Y : float
        Mass fraction of helium.
    n_capas: int
         Numero de capas del Number of layers in the model.
    Rini : float
        Initial radius from which the model calculation is made [10^10 cm].
    tol : float
        the tolerance at which two values are considered equal.
    tol2 : float
        tolerance of the temperature.
    g : float
        gamma constant.
    r_inf :
        Lower radius of the mesh [10^10 cm]
    r_sup :
        Upper radius of the mesh [10^10 cm]
    L_inf :
        Lower brightness of the mesh [10^33 erg/s]
    L_sup :
        Upper brightness of the mesh [10^33 erg/s]
    itL :
        Number of brightness values in the mesh
    itr :
        Number of radius values in the mesh
    it : int
        The number of meshes created, resulting in a greater zoom level with each iteration.
    debug : bool
        If this is true, the function will print both the model results and the colour graphs and maps.
    debug2 : bool
        If this is true, the function will generate, for each grid cell created, a plot showing the error and the central temperature in relation to the radius-luminosity pair, if not, just for the last one
    debug_opt_param : bool
        If this is true, the function will print the optimal parameters of the model
    debug_spectral : bool
        If this is true, the function will print the spectral star's type information.
    debug_HR : bool
        If this is true, the function will print the HR diagram position of the star.
    all_debug : bool
        If it is false, all the others debug are false and so the function does not print anything


    Returns
    -------
    E_fin :
        The generation process of each layer; None means --, 1 means PP and 0 means CN
    fase_fin :
        List with the fase of each layer; None==ˆˆˆˆˆˆ, 0==RADIAT, 1==CONVEC, 2==CENTRO and 3==INICIO
    r_fin :
        List of the radios.
    P_fin :
        List of the pressures.
    T_fin :
        List of the temperatures.
    L_fin :
        List of the brightness.
    M_fin :
        List of the masses.
    par_rad_fin :
        List showing the radiative parameter for each layer.
    error :
        Relative model error.
    Teff :
        Effective temperature of the star.
    mallado_min [R,L,T,Error]:
        List showing the radius, luminosity, core temperature and model error.
    """
    k = 1.380649 * 10 ** -16
    H = (6.02214076 * 10 ** 23) ** -1  # reciprocal of Avogadro's number
    Z = 1 - X - Y
    mu = 1 / (2 * X + 3 * Y / 4 + Z / 2)

    if Rini is None:
        Rini = 0.9 * Rtot

    h = -Rini / n_capas

    estrella = {
        'k': k,
        'H': H,
        'g': g,
        'tol': tol,
        'tol2': tol2,
        'X': X,
        'Y': Y,
        'Z': Z,
        'Mtot': Mtot,
        'Ltot': Ltot,
        'Tc': Tc,
        'Rtot': Rtot,
        'n_capas': n_capas,
        'Rini': Rini,
        'mu': mu,
        'h': h
    }
    mallado_min=refinement(r_inf,r_sup,L_inf,L_sup,itL,itr,it,estrella,debug=debug, debug2=debug2, debug_opt_param=debug_opt_param)

    E_fin, fase_fin, r_fin, P_fin, T_fin, L_fin, M_fin, par_rad_fin, error = model(Mtot, mallado_min[1], mallado_min[2], mallado_min[0],X ,Y,n_capas=n_capas, tol=tol)
    return E_fin, fase_fin, r_fin, P_fin, T_fin, L_fin, M_fin, par_rad_fin, error, mallado_min