from Auxiliary_functions import gen_rate

def radiative3(M,dM: list[float],Pcal: float,Tcal: float,r: list[float],mu: float,h: float,i: int):
    """
    Calculation of calibration mass.

    This value will be compared with the estimated mass, becoming the estimated value if they do not match and recalculating a new calibration value.

    Parameters
    ----------
    M :
        Mass c.
    dM : list[float]
        Mass gradient.
    Pcal : float
        Pressure calibration [10^15 din cm^-2].
    Tcal : float
        Temperature calibration [10^7 K].
    r : list[float]
        Radio [10^10 cm] array.
    mu : float
        Average molecular weight.
    h : float
        Step size [10^10 cm].
    i : int
        Current layer index (uses i, i-1, i-2).

    Returns
    -------
    [dM,Mcal]
        dM :
            List of errors in the mass.
        Mcal :
            List of calibration mass.

    Notes
    -----
    M colud be either a float or a list, in this second case, we just take the float that corresponds to our layer.
    """
    if isinstance(M,list):
        M=M[i]
    else:
        pass
    dM[i+1]=0.01523*mu*Pcal/Tcal*(r[i]+h)**2
    Mcal=M+h*dM[i+1]-h/2*(dM[i+1]-dM[i])
    return dM,Mcal


def radiative4(P: list[float],Pest: float,dP: list[float],Test: float,M: float,r: list[float],h: float,mu: float,i: int):
    """
    Calculation of calibration pressure.

    This value will be compared with the estimated pressure, becoming the estimated value if they do not match and recalculating a new calibration value.

    Parameters
    ----------
    P : list[float]
        Pressure [10^15 din cm^-2] array.
    Pest : float
        Pressure estimation [10^15 din cm^-2].
    dP : list[float]
        Gradient pressure array.
    Test : float
        Temperature estimation [10^7 K].
    M : float
        Mass [10^33 g].
    r : list[float]
        Radio [10^10 cm] array.
    h : float
        Step size [10^10 cm].
    mu : float
        Average molecular weight.
    i : int
        Current layer index (uses i, i-1, i-2).

    Returns
    -------
    [dP,Pcal]
        dP :
            List of errors in the pressure.
        Pcal:
            List of calibration pressure.
    """
    dP[i+1]=-8.084*mu*Pest/Test*M/(r[i]+h)**2
    Pcal=P[i]+h*dP[i+1]-h*(dP[i+1]-dP[i])/2
    return dP, Pcal

def radiative6(L: list[float],dL: list[float],Test: float,X: float,Z: float,Pcal:float,mu: float,r: list[float],i: int,h: float):
    """
    Calculation of calibration brightness.

    This value will be compared with the estimated brightness, becoming the estimated value if they do not match and recalculating a new calibration value.

    Parameters
    ----------
    L : list[float]
        Brightness [10^33 erg/s] array.
    dL : list[float]
        Brightness derivative [10^33 erg/s] array.
    Test : float
        Temperature estimation [10^7 K].
    X : float
        Mass fraction of hydrogen.
    Z : float
        Mass fraction of heavy metals.
    Pcal : float
        Pressure calibration [10^15 din cm^-2].
    mu : float
        Average molecular weight.
    r : list[float]
        Radio [10^10 cm] array.
    i : int
        Current layer index (uses i, i-1, i-2).
    h : float
        Step size [10^10 cm].

    Returns
    -------
    [dL,Lcal]
        dL:
            List of errors in the brightness.
        Lcal:
            List of calibration brightness.
    """
    eps,eps1,nu,X1,X2,E=gen_rate(Test,X,Z)
    dL[i+1]=0.01845*eps1*X1*X2*10**nu*mu**2*Pcal**2*Test**(nu-2)*(r[i]+h)**2
    Lcal=L[i]+h*dL[i+1]-h/2*(dL[i+1]-dL[i])-h/12*(dL[i+1]-2*dL[i]+dL[i-1])
    return dL,Lcal

def radiative7(T,dT,Pcal,Test,L,r,X,Z,mu,h,i):
    """
    Calculation of calibration temperature.

    This value will be compared with the estimated temperature, becoming the estimated value if they do not match and recalculating a new calibration value.

    Parameters
    ----------
    T : list[float]
        Temperature [10^7 K] array.
    dT : list[float]
        Gradient temperature array.
    Pcal : float
        Pressure calibration [10^15 din cm^-2].
    Test : float
        Temperature estimation [10^7 K].
    L : float
        Brightness [10^33 erg s^-1].
    r : list[float]
        Radio [10^10 cm] array.
    X : float
        Mass fraction of hydrogen.
    Z : float
        Mass fraction of heavy metals.
    mu : float
        Average molecular weight.
    h : float
        Step size [10^10 cm].
    i : int
        Current layer index (uses i, i-1, i-2).

    Returns
    -------
    [dT,Tcal]
        dT:
            List of errors in the temperature.
        Tcal:
            List of calibration temperature.
    """
    dT[i+1]=-0.01679*Z*(1+X)*mu**2*Pcal**2/Test**8.5*L/(r[i]+h)**2
    Tcal=T[i]+h*dT[i+1]-h/2*(dT[i+1]-dT[i])
    return dT,Tcal

def radiative9(Tcal: float,Pcal: float,dT: list[float],dP: list[float],i: int):
    """
    Calculation of the parameter that indicates whether transport is radioactive or convective (it is convective if it is less than or equal to 2.5).

    Parameters
    ----------
    Tcal : float
        Temperature calibration [10^7 K].
    Pcal : float
        Pressure calibration [10^15 din cm^-2].
    dT : list[float]
        Gradient temperature array.
    dP : list[float]
        Gradient pressure array.
    i : int
        Current layer index (uses i, i-1, i-2).
    """
    return Tcal/Pcal*dP[i+1]/dT[i+1]