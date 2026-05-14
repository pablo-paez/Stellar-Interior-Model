from .surface_integration import outer_layers, surf_int
from .centre_integration import first_layers, centre_int
import math

def border(parametro_rad, i_2: int,r,P,T,L,M):
    """
    Calculate the model values at the boundary between the radioactive and convective zones.

    Parameters
    ----------
    parametro_rad: list
        The different values of the radioactive parameter that indicate when the energy transport method becomes convective, when it drops below two and a half.
    i_2 : int
        The index of the convective zone.

    Returns
    -------
    i_limite :
        The index in which the radioactive parameter wuold be 2.5
    P_limite :
        The pressure at that index.
    T_limite :
        The temperature at that index.
    L_limite :
        The brightness at that index.
    M_limite :
        The mass at that index.
    r_limite :
        The radio at that index.
    """
    m=parametro_rad[i_2]-parametro_rad[i_2-1]
    n=parametro_rad[i_2]-i_2*m
    i_limite=(2.5-n)/m
    m_P=P[i_2]-P[i_2-1]
    m_T=T[i_2]-T[i_2-1]
    m_L=L[i_2]-L[i_2-1]
    m_M=M[i_2]-M[i_2-1]
    m_r=r[i_2]-r[i_2-1]
    n_P=P[i_2-1]-(i_2-1)*m_P
    n_T=T[i_2-1]-(i_2-1)*m_T
    n_L=L[i_2-1]-(i_2-1)*m_L
    n_M=M[i_2-1]-(i_2-1)*m_M
    n_r=r[i_2-1]-(i_2-1)*m_r
    P_limite=i_limite*m_P+n_P
    T_limite=i_limite*m_T+n_T
    L_limite=i_limite*m_L+n_L
    M_limite=i_limite*m_M+n_M
    r_limite=i_limite*m_r+n_r
    return i_limite,P_limite, T_limite, L_limite, M_limite, r_limite


def model(Mtot: float,Ltot: float,Tc: float,Rtot: float,X: float,Y: float,Rini=None,n_capas: int=100,tol=0.0001):
    """
    Function that calculates the model for the parameters given, returning all the values of the model.

    Parameters
    ----------
    Mtot: float
        Total mass.
    Ltot: float
        Total brightness.
    Tc: float
        Central temperature.
    Rtot: float
        Total radio.
    X : float
        Mass fraction of hydrogen.
    Y : float
        Mass fraction of helium.
    Rini : float
        Initial radius from which the model calculation is made.
    n_capas: int
         Numero de capas del Number of layers in the model.
    tol : float
        tolerance.
    g : float
        gamma constant.

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
    error:
        Relative model error.
    """
    if Rini is None:
        Rini = 0.9 * Rtot
    Z=1-X-Y
    mu=1/(2*X+3*Y/4+Z/2)
    h=-Rini/n_capas
    n_sobra=int(0.1*Rtot/abs(h))#number of layers I need to calculate to complete the star with the outer layers.
    par_rad_fin=[0]*n_sobra
    E_fin=[None]*n_sobra#Save the energy generation process (None means --, 1 means PP and 0 means CN)
    fase_fin=[None]*n_sobra+[3]*3+[0]*(n_capas-5)+[2]*3# None==ˆˆˆˆˆˆ, 0==RADIAT, 1==CONVEC, 2==CENTRO and 3==INICIO
    r,P,T,L,M=outer_layers(Rtot,Rini,h,mu,X,Z,Mtot,Ltot)
    r_fin=r[::-1]
    P_fin=P[::-1]
    T_fin=T[::-1]
    L_fin=L[::-1]
    M_fin=M[::-1]
    r,T,P,L,M,parametro_rad,E,i_frontera,K=surf_int(n_capas, Mtot, Ltot, Rini, Rtot, mu, X, Z, tol)
    r_fin.extend(r)
    P_fin.extend(P)
    T_fin.extend(T)
    L_fin.extend(L)
    M_fin.extend(M)
    par_rad_fin.extend(parametro_rad)
    E_fin.extend(E)
    i_limite,P_limite, T_limite, L_limite, M_limite, r_limite=border(parametro_rad, i_frontera,r,P,T,L,M)
    down=[r_limite,P_limite,T_limite,L_limite,M_limite]
    i_max=math.ceil(abs(r_limite/h))#It is the maximum value that i can take from the centre so that the radius is not greater than the limit radius.
    r_fin=r_fin[:i_frontera+n_sobra]
    T_fin=T_fin[:i_frontera+n_sobra]
    P_fin=P_fin[:i_frontera+n_sobra]
    L_fin=L_fin[:i_frontera+n_sobra]
    M_fin=M_fin[:i_frontera+n_sobra]
    par_rad_fin=par_rad_fin[:i_frontera+n_sobra]
    E_fin=E_fin[:i_frontera+n_sobra]
    fase_fin=fase_fin[:i_frontera+n_sobra]
    fase_con=[2]*3+[1]*(i_max-3)
    par_rad_con=[0]*i_max
    rcon,Mcon,Tcon,Lcon,Pcon,Econ, dMcon,dTcon,dLcon,dPcon=first_layers(i_max+1,Tc, -h, K, X, Z, mu)
    rcon, Tcon, Pcon, Lcon, Mcon, Econ, rext,Mext,Text,Lext,Pext=centre_int(i_max+1,rcon,Tcon,Pcon,Mcon,Lcon,Econ,mu,-h,dTcon, dMcon, dLcon, K, tol, X, Z)
    m=rext-rcon[-1]
    n=rext-i_max*m
    i_limite_con=(r_limite-n)/m
    up=[r_limite,i_limite_con*(Pext-Pcon[-1])+Pext-i_max*(Pext-Pcon[-1]),i_limite_con*(Text-Tcon[-1])+Text-i_max*(Text-Tcon[-1]),i_limite_con*(Lext-Lcon[-1])+Lext-i_max*(Lext-Lcon[-1]),i_limite_con*(Mext-Mcon[-1])+Mext-i_max*(Mext-Mcon[-1])]
    errores_relativos = []
    for i in range(len(down)-1):
        errores_relativos.append(abs(up[i+1]/down[i+1]-1)*100)
    error=math.sqrt(sum(x**2 for x in errores_relativos))
    r_fin.extend(reversed(rcon))
    T_fin.extend(reversed(Tcon))
    P_fin.extend(reversed(Pcon))
    L_fin.extend(reversed(Lcon))
    M_fin.extend(reversed(Mcon))
    par_rad_fin.extend(par_rad_con)
    E_fin.extend(reversed(Econ))
    fase_fin.extend(reversed(fase_con))
    return E_fin,fase_fin,r_fin,P_fin,T_fin,L_fin,M_fin,par_rad_fin,error