from Auxiliary_functions import Erel, gen_rate
from radiative import radiative3
from surface_integration_aux import presandtemp_layer_noP

def conv_center_eq(mu: float,T: float,r: float,K: float,X: float,Z: float):
    """
    Calculation of parameters in the central layers of a star assuming convective transport.

    Parameters
    ----------
    mu : float
        Average molecular weight.
    T : float
        Central temperature of the star [10^7 K].
    r : float
        Radio for which the parameters are calculated [10^10 cm].
    K : float
        Politropic constant [10^(15-7·gamma/(gamma-1)) din cm^-2 K^(gamma/(gamma-1))].
    X : float
        Mass fraction of hydrogen.
    Z : float
        Mass fraction of heavy metals.

    Returns
    -------
    tuple
        (M, L, Tcon, Pcon)
        - M: mass interior [10^33 g]
        - L: luminosity interior [10^33 erg/s]
        - Tcon: convective temperature [10^7 K]
        - Pcon: central pressure [10^15 din cm^-2]
    """
    Tcon=T-0.008207*mu**2*K*T**1.5*r**2
    eps, eps1,nu,X1,X2,_=gen_rate(Tcon,X,Z)
    M=0.005077*mu*K*T**1.5*r**3
    L=0.00615*eps1*X1*X2*10**nu*mu**2*K**2*T**(3+nu)*r**3
    Pcon=K*Tcon**2.5
    return M, L, Tcon, Pcon

def first_layers(i_max: int,Tc: float, h: float, K:float, X:float, Z:float, mu: float):
    """This function calculates the three innermost layers of the star.

    Parameters
    ----------
    i_max : int
        It is the maximum value that i can take from the centre so that the radius is not greater than the limit radius where the layer goes to radioactive.
    Tc : float
        Central temperature.
    h : float
        Step between layers in positive.
    K : float
        Polytropic constant.
    X : float
        Mass fraction of hydrogen.
    Z : float
        Mass fraction of heavy metals.
    mu : float
        Average molecular weight.

    Returns
    -------
    rcon :
        Vector showing the positions of the star’s layers.
    Mcon :
        Vector showing the mass of the star’s layers.
    Tcon :
        Vector showing the temperature of the star’s layers.
    Lcon :
        Vector showing the brightness of the star’s layers.
    Pcon :
        Vector showing the pressure of the star’s layers.
    Econ :
        Vector showing the energy generation process of the star’s layers.
    dMcon :
        Vector showing the mass gradient in the star’s layers.
    dTcon :
        Vector showing the temperature gradient in the star’s layers.
    dLcon :
        Vector showing the brightness gradient in the star’s layers.
    dPcon :
        Vector showing the pressure gradient in the star’s layers.
    """
    Mcon=[None]*(i_max)
    Tcon=[None]*(i_max)
    Lcon=[None]*(i_max)
    Pcon=[None]*(i_max)
    rcon=[None]*(i_max)
    dMcon=[None]*(i_max)
    dTcon=[None]*(i_max)
    dLcon=[None]*(i_max)
    dPcon=[None]*(i_max)
    Econ=[None]*3
    for i in range(3):
        rcon[i]=i*h
        M1,L1,T1,P1=conv_center_eq(mu, Tc, rcon[i], K,X,Z)
        Mcon[i]=M1
        Tcon[i]=T1
        Lcon[i]=L1
        Pcon[i]=P1
        eps,eps1,nu,X1,X2,_=gen_rate(T1,X,Z)
        if i==0:
            dMcon[i]=0
            dTcon[i]=0
            dLcon[i]=0
            dPcon[i]=0
        else:
            dMcon[i]=0.01523*mu*K*T1**1.5*rcon[i]**2
            dTcon[i]=-3.234*mu*M1/rcon[i]**2
            dPcon[i]=-8.084*mu*K*T1**1.5*M1/rcon[i]**2
            dLcon[i]=0.01845*eps1*X1*X2*10**nu*mu**2*K**2*T1**(3+nu)*rcon[i]**2
    return rcon,Mcon,Tcon,Lcon,Pcon,Econ, dMcon,dTcon,dLcon,dPcon

def centre_int(i_max:int,rcon,Tcon,Pcon,Mcon,Lcon,Econ,mu: float,h: float,dTcon, dMcon, dLcon, K: float, tol: float, X: float, Z: float):
    """This function calculates the layers integrating from the centre.

    Parameters
    ----------
    i_max : int
        It is the maximum value that i can take from the centre so that the radius is not greater than the limit radius where the layer goes to radioactive.
    rcon :
        Radio list.
    Tcon :
        Temperature list.
    Pcon :
        Pressure list.
    Mcon :
        Mass list.
    Lcon :
        Brightness list.
    Econ :
        Save the energy generation process (None means --, 1 means PP and 0 means CN)
    mu : float
        Average molecular weight.
    h : float
        Step between layers in positive.
    dTcon :
        List of temperature errors.
    dMcon :
        List of mass errors.
    dLcon :
        List of brightness errors.
    K : float
        Polytropic constant.
    tol : float
        Tolerance.
    X : float
        Mass fraction of hydrogen.
    Z : float
        Mass fraction of heavy metals.

    Returns
    -------
    rcon :
        Vector showing the positions of the star’s layers.
    Tcon :
        Vector showing the temperature of the star’s layers.
    Pcon :
        Vector showing the pressure of the star’s layers.
    Lcon :
        Vector showing the brightness of the star’s layers.
    Mcon :
        Vector showing the mass of the star’s layers.
    rcon[-1] :
        The final radius value is not used, as it lies within the area calculated from the surface; it will be useful for calculating the values at the boundary from the interior.
    Mcon[-1] :
        The final mass value is not used, as it lies within the area calculated from the surface; it will be useful for calculating the values at the boundary from the interior.
    Tcon[-1] :
        The final temperature value is not used, as it lies within the area calculated from the surface; it will be useful for calculating the values at the boundary from the interior.
    Lcon[-1] :
        The final brightness value is not used, as it lies within the area calculated from the surface; it will be useful for calculating the values at the boundary from the interior.
    Pcon[-1] :
        The final pressure value is not used, as it lies within the area calculated from the surface; it will be useful for calculating the values at the boundary from the interior.
    """
    i=1
    while i<i_max-2:
        i+=1
        rcon[i]=h*i
        Test=presandtemp_layer_noP(Tcon, dTcon, h, i)
        Pest=K*Test**2.5
        dMcon,Mcal=radiative3(Mcon, dMcon, Pest, Test, rcon, mu, h, i)
        dTcon[i+1]=-3.234*mu*Mcal/((i+1)*h)**2
        Tcal=Tcon[i]+h*dTcon[i+1]-h/2*(dTcon[i+1]-dTcon[i])#r=0 only for i=0, which will not happen here
        while not Erel(Tcal,Test,tol):
            Test=Tcal
            Pest=K*Test**2.5
            dMcon,Mcal=radiative3(Mcon, dMcon, Pest, Test, rcon, mu, h, i)
            dTcon[i+1]=-3.234*mu*Mcal/((i+1)*h)**2
            Tcal=Tcon[i]+h*dTcon[i+1]-h/2*(dTcon[i+1]-dTcon[i])#r=0 only for i=0, which will not happen here
        Pcal=K*Tcal**2.5
        _,eps1,nu,X1,X2,_=gen_rate(Tcal,X,Z)
        dLcon[i+1]=0.01845*eps1*X1*X2*10**nu*mu**2*Pcal**2*Tcal**(nu-2)*(h*(i+1))**2
        Lcal=Lcon[i]+h*dLcon[i+1]-h/2*(dLcon[i+1]-dLcon[i])-h/12*(dLcon[i+1]-2*dLcon[i]+dLcon[i-1])
        Pcon[i+1]=Pcal
        Tcon[i+1]=Tcal
        Mcon[i+1]=Mcal
        Lcon[i+1]=Lcal
        rcon[i+1]=h*(i+1)
        Econ.append(gen_rate(Tcal,X,Z)[-1])
    return rcon[:-1], Tcon[:-1], Pcon[:-1], Lcon[:-1], Mcon[:-1], Econ[:-1], rcon[-1],Mcon[-1],Tcon[-1],Lcon[-1],Pcon[-1]