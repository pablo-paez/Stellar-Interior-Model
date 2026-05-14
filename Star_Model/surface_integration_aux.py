import math

def der_pres_rad(r: float,M: float,P: float,T: float,mu: float) -> float:
    """
    Calculation of the derivative of pressure with respect to position in the radioactive case.

    Parameters
    ----------
    r : float
        Star radius at which the derivative is calculated [10^10 cm].
    M : float
        Mass contained within that radius [10^33 g].
    P : float
        Pressure contained within that radius [10^15 din cm^-2].
    T : float
        Temperature contained within that radius [10^7 K].
    mu : float
        Average molecular weight.

    Returns
    -------
    float
        Derivative of the pressure [10^5 din cm^-3].
    """
    return -8.084*mu*P/T*M/r**2

def der_mass_rad(P: float,T: float,r: float,mu: float) -> float:
    """
    Calculation of the mass gradient in radioactive cases.

    Parameters
    ----------
    P: float
        Pressure [10^15 din cm^-2].
    T : float
        Temperature [10^7 K].
    r : float
        Radio [10^10 cm].
    mu : float
        Average molecular weight.

    Returns
    -------
    float
        Gradient of the mass.
    """
    return 0.01523*mu*P/T*r**2

def der_temp_rad(r: float,L: float,P: float,T: float,X: float,Z: float,mu: float) -> float:
    """
    Calculation of the derivative of temperature with respect to position in the radioactive case.

    Parameters
    ----------
    r : float
        Star radius at which the derivative is calculated [10^10 cm].
    L : float
        Brightness contained within that radius [10^33 erg/s].
    P : float
        Pressure contained within that radius [10^15 din cm^-2].
    T : float
        Temperature contained within that radius [10^7 K].
    X : float
        Mass fraction of hydrogen.
    Z : float
        Mass fraction of heavy metals.
    mu : float
        Average molecular weight.

    Returns
    -------
    float
        Derivative of the temperature [10^-3 k/cm].
    """
    return -0.01679*Z*(1+X)*mu**2*P**2/T**8.5*L/r**2

def temp_rad(r: float,Rtot: float,mu: float,Mtot: float) -> float:
    """
    Calculation of the temperature of the outer layers of a star, assuming radioactive energy transport.

    Parameters
    ----------
    r : float
        Star radius at which the temperature is calculated [10^10 cm].
    Rtot : float
        Star radius [10^10 cm].
    mu : float
        Average molecular weight.
    Mtot : float
        Total mass of the star [10^33 g].

    Returns
    -------
    float
        Temperature [10^7 K].
    """
    A1=1.9022*mu*Mtot
    return A1*(1/r-1/Rtot)
def pres_rad(T: float,mu: float,X: float,Z: float,Mtot: float,Ltot: float) -> float:
    """
    Calculation of the pressure of the outer layers of a star, assuming radioactive energy transport.

    Parameters
    ----------
    T : float
        Star temperature at which the pressure is calculated [10^7 K].
    mu : float
        Average molecular weight.
    X : float
        Mass fraction of hydrogen.
    Z : float
        Mass fraction of heavy metals.
    Mtot : float
        Total mass of the star [10^33 g].
    Ltot : float
        Total brightness of the star [10^33 erg/s].

    Returns
    -------
    float
        Pressure [10^15 din cm^-2].
    """
    A2=10.645*math.sqrt(Mtot/(Ltot*mu*Z*(1+X)))
    return A2*T**4.25

def presandtemp_layer(P: list[float], T: list[float], dP: list[float],dT: list[float], h: float, i: int) -> tuple[float, float]:
    """
    Calculation of estimated pressure and temperature based on the previous layer.

    We use the difference method, at most of second order, to calculate the evolution of temperature and pressure from the values already calculated.

    Parameters
    ----------
    P, T : list[float]
        Pressure [10^15 dyne/cm²], Temperature [10^7 K] arrays.
    dP, dT : list[float]
        Gradients dP/dr , dT/dr arrays.
    h : float
        Step size Δr [10^10 cm].
    i : int
        Current layer index (uses i, i-1, i-2).

    Returns
    -------
    tuple[float, float]
        (Pest, Test): predicted values at i+1.
    """
    Pest=P[i]+h*dP[i]+h/2*(dP[i]-dP[i-1])+5/12*h*(dP[i]-2*dP[i-1]+dP[i-2])
    Test=T[i]+h*dT[i]+h/2*(dT[i]-dT[i-1])
    return Pest, Test

def presandtemp_layer_noP(T: list[float],dT: list[float],h: float,i: int) -> float:
    """
    Calculation of estimated temperature based on the previous layer.

    We use the difference method, at most of second order, to calculate the evolution of temperature and pressure from the values already calculated.

    Parameters
    ----------
    T : list[float]
        Temperature [10^7 K] arrays.
    dT : list[float]
        Gradient dT/dr arrays.
    h : float
        Step size Δr [10^10 cm].
    i : int
        Current layer index (uses i, i-1, i-2).

    Returns
    -------
    float
        Test: predicted values at i+1.
    """
    Test=T[i]+h*dT[i]+h/2*(dT[i]-dT[i-1])
    return  Test