from .surface_integration_aux import der_pres_rad, der_temp_rad, temp_rad, pres_rad, presandtemp_layer, presandtemp_layer_noP
from .radiative import radiative3, radiative4, radiative6, radiative7, radiative9
from .Auxiliary_functions import Erel, gen_rate, polytropic

def outer_layers(Rtot: float,Rini: float,h: float,mu: float,X: float,Z: float,Mtot: float,Ltot: float):
    """This function calculates the outer layers of the model which cannot be calculated earlier.

    Parameters
    ----------
    Rtot : float
        Total radius of the star.
    Rini: float
        Radius from which the layers have been calculated.
    h: float
        Step between layers in negative.
    mu : float
        Average molecular weight.
    X : float
        Mass fraction of hydrogen.
    Z : float
        Mass fraction of heavy metals.
    Mtot : float
        Total mass of the star.
    Ltot : float
        Total brightness of the star.

    Returns
    -------
    [r,P,T,L,M]
    """
    n_sobra=int(0.1*Rtot/abs(h))
    i=0
    r=[None]*n_sobra
    P=[None]*n_sobra
    T=[None]*n_sobra
    M=[Mtot]*n_sobra
    L=[Ltot]*n_sobra
    while i<n_sobra:
        r[i]=Rini-h*(i+1)
        T[i]=temp_rad(r[i],Rtot,mu,Mtot)
        P[i]=pres_rad(T[i],mu,X,Z,Mtot,Ltot)
        i+=1
    return r,P,T,L,M


def surf_int(n_capas: int, Mtot: float, Ltot: float, Rini: float, Rtot: float, mu: float, X: float, Z: float, tol: float):
    """
    This function calculates the integration from the surface of the stellar model, stopping the radioactive code until the so-called radioactive parameter is less than 2.5 and then activating the convective one. The first three layers are calculated different to the others.

    Parameters
    ----------
    n_capas: int
        Number of layers in our model.
    Mtot : float
        Total mass of the star.
    Ltot : float
        Total brightness of the star.
    Rini: float
        Radius from which the layers have been calculated.
    Rtot : float
        Total radius of the star.
    mu : float
        Average molecular weight.
    X : float
        Mass fraction of hydrogen.
    Z : float
        Mass fraction of heavy metals.
    tol : float
        Relative tolerance for two values to be considered equal.

    Returns
    -------
    r :
        Vector showing the positions of the star’s layers.
    T :
        Vector showing the temperatures in the star's layers.
    P :
        Vector showing the pressures in the star's layers.
    L :
        Vector showing the brightness in the star's layers.
    M :
        Vector showing the masss in the star's layers.
    parametro_rad :
        Vector showing the radioactive parameter in the star's layers.
    E :
        Vector showing the energy generation process in the star's layers.
    i_2 :
        Layer in which radiation is converted into convection.
    K :
        Polytropic constant.
    """
    h=-Rini/n_capas
    r=[None]*(n_capas+1)
    P=[None]*(n_capas+1)
    T=[None]*(n_capas+1)
    M=[Mtot,Mtot,Mtot]+[None]*(n_capas-2)
    L=[Ltot,Ltot,Ltot]+[None]*(n_capas-2)
    dP=[None]*(n_capas+1)
    dT=[None]*(n_capas+1)
    dM=[0,0,0]+[None]*(n_capas-2)
    dL=[0,0,0]+[None]*(n_capas-2)
    parametro_rad=[0]*3 #n+1 (Is 0 in the first 3 layers)
    E=[None]*3 #Save the energy generation process (None has been set in the first three layers.)
    for i in range(0,3):
        r[i]=Rini+h*i
        T[i]=temp_rad(r[i],Rtot,mu,Mtot)
        P[i]=pres_rad(T[i],mu,X,Z,Mtot,Ltot)
        dP[i]=der_pres_rad(r[i],Mtot,P[i],T[i],mu)
        dT[i]=der_temp_rad(r[i],Ltot,P[i],T[i],X,Z,mu)
    i-=1 #Given that the last calculated layer is not valid
    while i<n_capas:
        i+=1
        r[i]=Rini+h*i
        Pest,Test=presandtemp_layer(P, T, dP, dT, h, i)
        dM,Mcal=radiative3(M,dM,Pest,Test,r,mu,h,i)
        dP,Pcal=radiative4(P, Pest, dP, Test, Mcal, r, h, mu, i)
        while not Erel(Pcal,Pest,tol):
            Pest=Pcal
            dM,Mcal=radiative3(M,dM,Pest,Test,r,mu,h,i)
            dP,Pcal=radiative4(P, Pest, dP, Test, Mcal, r, h, mu, i)
        dL,Lcal=radiative6(L,dL, Test, X, Z, Pcal, mu, r, i, h)
        dT,Tcal=radiative7(T, dT, Pcal, Test, Lcal, r, X, Z, mu, h, i)
        while not Erel(Tcal,Test,tol):
            Test=Tcal
            Pest=Pcal
            dM,Mcal=radiative3(M,dM,Pest,Test,r,mu,h,i)
            dP,Pcal=radiative4(P, Pest, dP, Test, Mcal, r, h, mu, i)
            while not Erel(Pcal,Pest,tol):
                Pest=Pcal
                dM,Mcal=radiative3(M,dM,Pest,Test,r,mu,h,i)
                dP,Pcal=radiative4(P, Pest, dP, Test, Mcal, r, h, mu, i)
            dL,Lcal=radiative6(L,dL, Test, X, Z, Pcal, mu, r, i, h)
            dT,Tcal=radiative7(T, dT, Pcal, Test, Lcal, r, X, Z, mu, h, i)
        n=radiative9(Tcal, Pcal, dT, dP, i)
        if n>2.5:
            P[i+1]=Pcal
            T[i+1]=Tcal
            M[i+1]=Mcal
            L[i+1]=Lcal
            if Tcal<0.4:
                E.append(None)
            else:
                E.append(gen_rate(Tcal,X,Z)[-1])
            parametro_rad.append(n)
        else:
            dP[i+1]=None
            dT[i+1]=None
            dL[i+1]=None
            dM[i+1]=None
            K=polytropic(Tcal,Pcal) #It will be a constant in the next phase.
            break
    i-=1
    i_2=i+2
    while i<n_capas:
        i+=1
        r[i]=Rini+h*i
        Test=presandtemp_layer_noP(T, dT, h, i)
        if Test<0:
            Pest=0
        else:
            Pest=K*Test**2.5
        dM,Mcal=radiative3(M, dM, Pest, Test, r, mu, h, i)
        if i<n_capas-1:#equivalent to r[i+1]=0
            dT[i+1]=-3.234*mu*Mcal/(r[i]+h)**2
            Tcal=T[i]+h*dT[i+1]-h/2*(dT[i+1]-dT[i])
        else:
            Tcal=Test
        while not Erel(Tcal,Test,tol):
            Test=Tcal
            if Test<0:
                Pest=0
            else:
                Pest=K*Test**2.5
            dM,Mcal=radiative3(M, dM, Pest, Test, r, mu, h, i)
            if r[i]+h>0:
                dT[i+1]=-3.234*mu*Mcal/(r[i]+h)**2
                Tcal=T[i]+h*dT[i+1]-h/2*(dT[i+1]-dT[i])
            else:
                Tcal=Test
        if Tcal<0:
            Pcal=0
        else:
            Pcal=K*Tcal**2.5
        dL,Lcal=radiative6(L,dL, Tcal, X, Z, Pcal, mu, r, i, h)
        if i<n_capas-1:
            P[i+1]=Pcal
            T[i+1]=Tcal
            M[i+1]=Mcal
            L[i+1]=Lcal
            E.append(gen_rate(Tcal,X,Z)[-1])
            parametro_rad.append(n)
        else:
            P[i+1]=Pcal
            T[i+1]=Tcal
            M[i+1]=Mcal
            L[i+1]=Lcal
            r[i+1]=Rini+h*(i+1)
            E.append(gen_rate(Tcal,X,Z)[-1])
            parametro_rad.append(n)
            break
    return r,T,P,L,M,parametro_rad,E,i_2,K