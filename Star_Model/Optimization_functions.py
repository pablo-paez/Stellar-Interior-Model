import numpy as np
import matplotlib.pyplot as plt
from simple_model import model

def golden_section(Mtot: float,Ltot: float,Rtot: float,X: float,Y: float,Ti: float=1.9,Tf: float=2.0,Rini: float=None,n_capas: int=100,tol: float=0.0001,tol2: float=1e-5):
    """
    This function uses the golden ratio method to calculate the minimum of the function.

    The method takes two values chosen from the boundary values, evaluating the function at those two points and taking the interval that includes the value that gives the smallest function value between c and d and that includes both c and d. This is done until the interval is less than 10^-5.

    Parameters
    ----------
    Mtot: float
        Total mass.
    Ltot: float
        Total brightness.
    Rtot: float
        Total radio.
    X : float
        Mass fraction of hydrogen.
    Y : float
        Mass fraction of helium.
    Ti: float
        Lower temperature limit.
    Tf: float
        Upper temperature limit.
    Rini : float
        Initial radio.
    n_capas : int
        Number of layers that the model will have, used to calculate the integration step h.
    tol : float
        Relative tolerance of the function evaluacion.
    tol2 : float
        Relative tolerance of the optimization method.
    down : list
        List of limit values for integration from the surface.

    Returns
    -------
    (a+b)/2
        The optimum temperature.
    modelo(Mtot,Ltot,(a+b)/2,Rtot,n_capas,X,Y,Rini,tol)
        The parameters of the star and the minimum relative error value at that temperature.
    """
    if Rini is None:
        Rini = 0.9 * Rtot
    phi = (1 + 5**0.5) / 2
    resphi = 2 - phi
    a=Ti
    b=Tf
    c = b - (b - a) / phi
    d = a + (b - a) / phi
    f_c = model(Mtot,Ltot,c,Rtot,X,Y,Rini=Rini,n_capas=n_capas,tol=tol)[-1]
    f_d = model(Mtot,Ltot,d,Rtot,X,Y,Rini=Rini,n_capas=n_capas,tol=tol)[-1]
    while abs(b - a) > tol2:
        if f_c < f_d:
            b = d
            d = c
            f_d = f_c
            c = a + resphi * (b - a)
            f_c = model(Mtot, Ltot, c, Rtot,X,Y,Rini=Rini,n_capas=n_capas,tol=tol)[-1]
        else:
            a = c
            c = d
            f_c = f_d
            d = b - resphi * (b - a)
            f_d = model(Mtot, Ltot, d, Rtot,X,Y,Rini=Rini,n_capas=n_capas,tol=tol)[-1]
    return (a + b) / 2, model(Mtot,Ltot,(a+b)/2,Rtot,X,Y,Rini=Rini,n_capas=n_capas,tol=tol)

def colours_maps(R, L, T, Error):
    """
    It generates two colour maps: one showing the central temperature and the other the relative error for the radiometric-luminance value pairs.

    Parameters
    ----------
    R : list[float]
        List of radio.
    L : list[float]
        List of brightness.
    T : list[float]
        List of temperatures.
    Error : list[float]
        List of error values.
    """
    R = np.array(R)
    L = np.array(L)
    T = np.array(T)
    Error = np.array(Error)

    r_vals = np.unique(R)
    L_vals = np.unique(L)

    Tz = np.full((len(L_vals), len(r_vals)), np.nan)
    Z  = np.full((len(L_vals), len(r_vals)), np.nan)

    for k in range(len(R)):
        i = np.searchsorted(L_vals, L[k])
        j = np.searchsorted(r_vals, R[k])
        Tz[i, j] = T[k]
        Z[i, j]  = Error[k]

    # Mínimo del error (en índices)
    idx = np.argmin(Z)
    i, j = np.unravel_index(idx, Z.shape)
    r_min_idx = j
    L_min_idx = i


    fig, axs = plt.subplots(1, 2, figsize=(13, 8))

    for ax, data, title in zip(
        axs,
        [Z, Tz],
        ["Minimum total error(%)", "Best central temperature"]
    ):
        im = ax.imshow(data, origin='lower', aspect='auto',
                       cmap='viridis', interpolation='none')

        ax.set_xticks(np.arange(len(r_vals)))
        ax.set_xticklabels([f"{v:.3f}" for v in r_vals])
        ax.set_yticks(np.arange(len(L_vals)))
        ax.set_yticklabels([f"{v:.2f}" for v in L_vals])

        ax.scatter(r_min_idx, L_min_idx, color='red', s=50)

        ax.set_title(title)
        ax.set_xlabel(r'Radio ($10^{10}$ cm)')
        ax.set_ylabel(r'Luminosidad ($10^{33}$ erg $s^{-1}$)')

        cbar = fig.colorbar(im, ax=ax, orientation='horizontal',
                            location='top', pad=0.08)
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_position('top')

    plt.tight_layout()
    plt.show()

def refinement(r_inf,r_sup,L_inf,L_sup,itL,itr,it: int,estrella,debug: bool =True, debug2: bool =True, debug_opt_param: bool=True):
    """
    A function that generates a radiance-based grid to find the value that minimises the relative error of the model.

    Parameters
    ----------
    r_inf :
        Lower radius of the mesh
    r_sup :
        Upper radius of the mesh
    L_inf :
        Lower brightness of the mesh
    L_sup :
        Upper brightness of the mesh
    itL :
        Number of brightness values in the mesh
    itr :
        Number of radius values in the mesh
    it : int
        The number of meshes created, resulting in a greater zoom level with each iteration.
    estrella
        Dictionary containing all the star’s parameters
    debug : bool
        If this is true, the function will generate, for the last grid cell created, a plot showing the error and the central temperature in relation to the radius-luminosity pair
    debug2 : bool
        If this is true, the function will generate, for each grid cell created, a plot showing the error and the central temperature in relation to the radius-luminosity pair, if not, just for the last one
    debug_opt_param : bool
        If this is true, the function will print the optimal parameters of the model

    Returns
    -------
    mallado_min [R,L,T,Error]:
        List showing the optimum radius, brightness, core temperature and relative error.
    """
    tol = estrella['tol']
    tol2 = estrella['tol2']
    X = estrella['X']
    Y = estrella['Y']
    Mtot = estrella['Mtot']
    n_capas = estrella['n_capas']
    for k in range(it):
        mallado=[]
        rmallado = []
        Lmallado = []
        Tmallado = []
        Errormallado = []

        if k==0:
            vr=np.linspace(r_inf,r_sup,itr)
            vL=np.linspace(L_inf,L_sup,itL)
            difr=(vr[1]-vr[0])*1
            difL=(vL[1]-vL[0])*1

        else:
            difr=(vr[1]-vr[0])*1
            difL=(vL[1]-vL[0])*1

            if mallado_min[0]<r_inf+difr:
                vr=np.linspace(r_inf,r_inf+4*difr,itr)
            elif mallado_min[0]>r_sup-difr:
                vr=np.linspace(r_sup-4*difr,r_sup,itr)
            else:
                vr=np.linspace(mallado_min[0]-difr,mallado_min[0]+difr,itr)

            if mallado_min[1]<L_inf+difL:
                vL=np.linspace(L_inf,L_inf+4*difL,itL)
            elif mallado_min[1]>L_sup-difL:
                vL=np.linspace(L_sup-4*difL,L_sup,itL)
            else:
                vL=np.linspace(mallado_min[1]-difL,mallado_min[1]+difL,itL)

        for i in range(len(vr)):
            Rtot_=vr[i]
            for j in range(len(vL)):
                Ltot_=vL[j]
                Rini_=Rtot_*0.9
                T,ress=golden_section(Mtot,Ltot_,Rtot_,X,Y,Ti=1.5,Tf=2,Rini=Rini_,n_capas=n_capas,tol=tol,tol2=tol2)
                Err=ress[-1]
                mallado.append([Rtot_,Ltot_,T,Err])
                rmallado.append(Rtot_)
                Lmallado.append(Ltot_)
                Tmallado.append(T)
                Errormallado.append(Err)
        if debug2:
            colours_maps(rmallado,Lmallado,Tmallado,Errormallado)
        mallado_min = min(mallado, key=lambda x: x[3])
    if debug and not debug2:
        colours_maps(rmallado,Lmallado,Tmallado,Errormallado)
    if debug_opt_param:
        print(f"Los parámetros óptimos de la estrella corresponden a una temperatura central {mallado_min[2]:.5f}·10⁷ K con una tolerancia de {tol2}, un radio total de {mallado_min[0]:.3f}·10¹⁰ cm con una tolerancia de {difr:.3f} y una luminosidad total de {mallado_min[1]:.2f}·10³³ erg/s con una tolerancia de {difL:.2f} todo ello con un error relativo del {mallado_min[3]:.5f}%.")
    return mallado_min