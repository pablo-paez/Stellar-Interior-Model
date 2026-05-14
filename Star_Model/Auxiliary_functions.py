def Erel(cal: float,est: float,tol: float) -> bool:
    """
    Checks whether two values have a relative error less than a given tolerance.

    Parameters
    ----------
    cal : float
        The calibration value (reference).
    est : float
        The estimated value to compare.
    tol : float
        Maximum allowed relative error.

    Returns
    -------
    bool
        Returns true if both values have a relative error less than the tolerance, false if not.

    Raises
    ------
    ValueError
        If tolerance is not positive.
    """
    if tol<=0:
        raise ValueError('Tolerance must be positive.')
    return abs(cal-est)/abs(cal)<tol

def opacity(rho: float,T: float, X: float, Z:float) -> float:
    """
    Calculates the opacity of the radiance.

    Parameters
    ----------
    rho : float
        The radiance.
    T : float
        The temperature.
    X : float
        Mass fraction of hydrogen.
    Z : float
        Mass fraction of heavy metals.

    Returns
    -------
    Opacity
    """
    return 4.34 / 3.162 * 10 ** 25 * Z * (1 + X) * rho / (T * 1e7) ** (7 / 2)

def rho_fun(mu: float,P: float,T: float) -> float:
    """
    Calculation of density using the equation of state for non-degenerate ideal gases.

    Parameters
    ----------
    mu : float
        Average molecular weight.
    P : float
        Star pressure [10^15 din cm^-2].
    T : float
        Star temperature [10^7 K].

    Returns
    -------
    float
        Density value for the given parameters [g cm^-3].
    """
    k=1.380649*10**-16 #erg/K
    H=(6.02214076*10**23)**-1 #g
    return mu*H*P*1e8/(k*T)

def gen_rate(T: float,X: float,Z: float):
    """
    Calculation of the star's energy generation rate.

    Calculation of the rate of energy generation by the two existing processes, proton-proton chain and CN cycle, given that the star only burns hydrogen by one of the two processes, both are compared and the one that gives the highest value is chosen.

    Parameters
    ----------
    T : float
        Star temperature [10^7 K].
    X : float
        Mass fraction of hydrogen.
    Z : float
        Mass fraction of heavy metals.

    Returns
    -------
    tuple
        (ε, eps1, ν, X_factor1, X_factor2, cycle_id)
        - ε: energy generation rate [erg/(g·s)].
        - eps1: normalization constant.
        - ν: temperature exponent.
        - X_factor1, X_factor2: fuel dependencies (X*X for pp, X*Z/3 for CNO).
        - cycle_id: 0 for CNO, 1 for pp.
    """
    T=T*10 #Our temperatures will be given in units of 10^7, and we want them in units of 10^6.
    if 12<=T<16:
        eps1CN=10**-22.2
        vCN=20
    elif 16<=T<22.5:
        eps1CN=10**-19.8
        vCN=18
    elif 22.5<=T<27.5:
        eps1CN=10**-17.1
        vCN=16
    elif 27.5<=T<36:
        eps1CN=10**-15.6
        vCN=15
    elif 36<=T<=50:
        eps1CN=10**-12.5
        vCN=13
    else:
        eps1CN=0
        vCN=0
    CN=eps1CN*X*Z/3*(T)**vCN
    if 4<=T<6:
        eps1pp=10**-6.84
        vpp=6
    elif 6<=T<9.5:
        eps1pp=10**-6.04
        vpp=5
    elif 9.5<=T<12:
        eps1pp=10**-5.56
        vpp=4.5
    elif 12<=T<16.5:
        eps1pp=10**-5.02
        vpp=4
    elif 16.5<=T<=24:
        eps1pp=10**-4.4
        vpp=3.5
    else:
        eps1pp=0
        vpp=0
    pp=eps1pp*X*X*(T)**vpp
    if CN>pp:
        return CN, eps1CN,vCN, X, Z/3, 0
    else:
        return pp, eps1pp, vpp, X, X, 1 #The number 0 or 1 is to identify the responsible cycle.

def polytropic(T: float, P: float, gamma: float = 5 / 3) -> float:
        """
        Calculation of the polytropic constant.

        Parameters
        ----------
        T : float
            Temperature [10^7 K].
        P : float
            Pressure [10^15 din cm^-2].
        gamma : float
            adiabatic index (if not specified, it is assumed to be 5/3, corresponding to a perfect monatomic gas).

        Returns
        -------
        float
            Politropic constant [10^(15-7·gamma/(gamma-1)) din cm^-2 K^(gamma/(gamma-1))].
        """
        return P / T ** (gamma / (gamma - 1))