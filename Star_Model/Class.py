import numpy as np
import matplotlib.pyplot as plt
from simple_model import border
from Optimized_model import optimized_model
from Optimization_functions import golden_section
import math
import pandas as pd
from Auxiliary_functions import gen_rate, rho_fun, opacity
from Charts import parameter_charts, HRDiagram

class StarModel:
    def __init__(self, Mtot: float, X: float, Y: float, Ltot: float, Rtot: float, Tc: float, unidades_solares:bool = False):
        """
        Initialise the star object with its basic physical parameters.

        Parameters
        ----------
        Mtot: float
            Total mass [10^33 g].
        X : float
            Mass fraction of hydrogen.
        Y : float
            Mass fraction of helium.
        Ltot: float
            Total brightness [10^33 erg/s].
        Rtot: float
            Total radio [10^10 cm].
        Tc: float
            Central temperature [10^7 K].
        unidades_solares: bool
            If true, the units of the mass, brightness and radius introduced will be transformed from solar units.
        """
        if X + Y > 1:
            raise ValueError("X + Y cannot exceed 1")
        if unidades_solares:
            Mtot = Mtot * 1.989 #Solar mass=1.989 *10**33 g
            Ltot = Ltot * 3.828
            Rtot = Rtot * 6.957
        self.Mtot = Mtot
        self.X = X
        self.Y = Y
        self.Ltot = Ltot
        self.Rtot = Rtot
        self.Tc = Tc
        self.Z = 1-X-Y
        self.mu = 1/(2*self.X + 3*self.Y/4 + self.Z/2)

        # Fields for storing the results after integration
        self.modelo = None #It contains all the layered information of the model within a single attribute
        self.Teff = None
        self.error_rel = None
        self.tipo_energia = None #Save the energy generation process (None means --, 1 means PP and 0 means CN)
        self.fase_transp = None # None==ˆˆˆˆˆˆ, 0==RADIAT, 1==CONVEC, 2==CENTRO and 3==INICI
        self.r           = None
        self.P           = None
        self.T           = None
        self.L           =None
        self.M           = None
        self.par_rad     = None
        self.rho         = None
        self.epsilon     = None
        self.kappa       = None
        self._optimizado = None #Internal attribute to check whether the optimisation of the radio and brightness parameters has been carried out.

    def __str__(self):
        return f"Estrella de M: {self.Mtot} [10³³ g] | R: {self.Rtot} [10¹⁰ cm] | L: {self.Ltot} [10³³ erg s⁻¹]"

    def __repr__(self):
        return f"ModeloEstelar(Mtot={self.Mtot}, Rtot={self.Rtot}, Ltot={self.Ltot})"



    def parameters(self):
        """
        Returns the parameters defining the star.

        Returns
        -------
        dict
            Dictionary containing the entries 'Mtot', 'Rtot', 'Ltot', 'Tc', 'X' and 'Y'.
        """
        return {
            'Mtot': self.Mtot,
            'Rtot': float(self.Rtot),
            'Ltot': float(self.Ltot),
            'Tc':   self.Tc,
            'X':    self.X,
            'Y':    self.Y
        }

    def redefine(self, Mtot=None, X=None, Y=None, Ltot=None, Rtot=None, Tc=None, unidades_solares:bool=False):
        """
        Redefines the instance's physical parameters.
        If a parameter is None, it retains the previous value.

        Parameters
        ----------
        Mtot: float
            Total mass [10^33 g].
        X : float
            Mass fraction of hydrogen.
        Y : float
            Mass fraction of helium.
        Ltot: float
            Total brightness [10^33 erg/s].
        Rtot: float
            Total radio [10^10 cm].
        Tc: float
            Central temperature [10^7 K].
        unidades_solares: bool
            If true, the units of the mass, brightness and radius introduced will be transformed from solar units.
        """
        # We update only the parameters provided
        if Mtot is not None:
            if unidades_solares:
                self.Mtot = Mtot * 1.989
            else:
                self.Mtot = Mtot
        if Rtot is not None:
            if unidades_solares:
                self.Rtot = Rtot * 6.957
            else:
                self.Rtot = Rtot
        if Ltot is not None:
            if unidades_solares:
                self.Ltot = Ltot * 3.828
            else:
                self.Ltot = Ltot
        if Tc is not None:   self.Tc = Tc
        # If we change X or Y, we must recalculate Z and mu
        if X is not None or Y is not None:
            # We use temporary variables to check before overwriting
            nuevo_X = X if X is not None else self.X
            nuevo_Y = Y if Y is not None else self.Y
            if nuevo_X + nuevo_Y > 1:
                raise ValueError("X + Y cannot exceed 1")
            self.X=nuevo_X
            self.Y=nuevo_Y
            self.Z = 1 - self.X - self.Y
            self.mu = 1 / (2*self.X + 0.75*self.Y + 0.5*self.Z)
        #We have removed the calculated data as the parameters have been changed
        self.modelo      = None
        self.Teff        = None
        self.error_rel   = None
        self.tipo_energia = None
        self.fase_transp = None
        self.r           = None
        self.P           = None
        self.T           = None
        self.L           = None
        self.M           = None
        self.par_rad     = None
        self.rho         = None
        self.epsilon     = None
        self.kappa       = None
        print(f"Parámetros actualizados: M={self.Mtot}, R={self.Rtot}, L={self.Ltot}, Tc={self.Tc}, X={self.X}, Y={self.Y}")

    def error(self):
        """
        Returns the total relative error (E) as a percentage, calculated as the sum of the squares of the differences between the parameters at the boundary between the convective and radiative zones.
        """
        if self.error_rel == None:
            raise("The stellar model has not been calculated for the star’s initial parameters.")

        return self.error_rel


    def calculate(self, optimizar:bool =True,mallado_final:bool =True,all_mallados:bool =False,parametros_optimos:bool =True,magnitudes_extra: bool = True, **kwargs):
        """
        Calculation of the complete stellar model based on the star’s initial parameters. This includes an option (optimise True/False) to decide whether to perform a search for the most optimal radius and luminosity parameters (in both cases, this search is carried out for the core temperature).

        Parameters
        ----------
        optimizar: bool
            If false, only the core temperature is optimised (study the possible inputs of the function seccion_aurea); if true, a search is also carried out for the most optimal radius-luminosity pair (study the possible inputs of the function modelo_completo).
        mallado_final: bool
            If optimizar is true, it displays the final radio-luminance mesh showing the optimal value.
        all_mallados: bool
            If true, displays the meshes for all the refinements carried out.
        parametros_optimos: bool
            If true, displays the parameters of the optimization searched.
        magnitudes_extra: bool
            If true, the function calculates the density, opacity and epsilon.
        """
        if not optimizar:
            res=golden_section(self.Mtot,self.Ltot,self.Rtot,self.X,self.Y,**kwargs)
            self.error_rel=res[-1][-1]
            self.modelo=res[1][:-1]
            if magnitudes_extra:
                rho=[None]*len(res[1][2])
                epsilon = [None] * len(res[1][2])
                kappa = [None] * len(res[1][2])
                for i in range(len(res[1][2])):
                    rho[i]=rho_fun(self.mu,res[1][3][i],res[1][4][i])
                    epsilon[i]=gen_rate(res[1][4][i],self.X,self.Z)[0]*rho[i]
                    kappa[i]=opacity(rho[i],res[1][4][i],self.X,self.Z)
                self.modelo = list(self.modelo) + [rho, epsilon, kappa]
                self.rho     = rho
                self.epsilon = epsilon
                self.kappa   = kappa
            self.tipo_energia = self.modelo[0]
            self.fase_transp = self.modelo[1]
            self.r       = self.modelo[2]
            self.P       = self.modelo[3]
            self.T       = self.modelo[4]
            self.L       = self.modelo[5]
            self.M       = self.modelo[6]
            self.par_rad = self.modelo[7]
        else:
            res=optimized_model(self.Mtot,self.Ltot,self.Tc,self.Rtot,self.X,self.Y,debug=mallado_final,debug_opt_param=parametros_optimos,debug2=all_mallados,**kwargs)
            self.error_rel=res[-2]
            self.modelo=res[:-2]
            if magnitudes_extra:
                rho=[None]*len(res[2])
                epsilon = [None] * len(res[2])
                kappa = [None] * len(res[2])
                for i in range(len(res[2])):
                    rho[i]=rho_fun(self.mu,res[3][i],res[4][i])
                    epsilon[i]=gen_rate(res[4][i],self.X,self.Z)[0]*rho[i]
                    kappa[i]=opacity(rho[i],res[4][i],self.X,self.Z)
                self.modelo = list(self.modelo) + [rho, epsilon, kappa]
                self.rho     = rho
                self.epsilon = epsilon
                self.kappa   = kappa
            self.tipo_energia = self.modelo[0]
            self.fase_transp = self.modelo[1]
            self.r       = self.modelo[2]
            self.P       = self.modelo[3]
            self.T       = self.modelo[4]
            self.L       = self.modelo[5]
            self.M       = self.modelo[6]
            self.par_rad = self.modelo[7]
            self.Rtot=res[-1][0]
            self.Ltot=res[-1][1]
            self.Tc=res[-1][2]
            self._optimizado = True


    def layers(self, to_csv=False, name=None):
        """
        Displays the model layers and, if desired, exports them to CSV.

        Parameters
        ----------
        to_csv : bool
            If True, save the data to a .csv file.
        name : str, optional
            CSV file name.
        """

        headers = f"{'E':>3} {'fase':>8} {'i':>6} {'r':>8} {'P':>12} {'T':>12} {'L':>12} {'M':>12} {'n+1':>12} {'ρ':>12} {'ε':>12} {'κ':>12}"
        print(headers)

        for j in range(len(self.modelo[0])):
            ener_label = 'PP' if self.tipo_energia[j] == 1 else 'CN' if self.tipo_energia[j] == 0 else '--'
            fase_label = {0: 'RADIAT', 1: 'CONVEC', 2: 'CENTRO', 3: 'INICIO'}.get(self.fase_transp[j], 'ˆˆˆˆˆˆ')

            print(
                f"{ener_label:>4} {fase_label:>8} {j-11:>5} "
                f"{self.r[j]:>11.5f} {self.P[j]:>12.7f} {self.T[j]:>12.7f} "
                f"{self.L[j]:>12.6f} {self.M[j]:>12.6f} {self.par_rad[j]:>12.6f} "
                f"{self.rho[j]:>12.6f} {self.epsilon[j]:>12.6f} {self.kappa[j]:>12.6f}"
            )

        if to_csv:
            datos_csv = {
                'Generacion_energia': self.modelo[0],
                'Fase': self.modelo[1],
                'Radio': self.modelo[2],
                'Presion': self.modelo[3],
                'Temperatura': self.modelo[4],
                'Luminosidad': self.modelo[5],
                'Masa': self.modelo[6],
                'Param_Rad': self.modelo[7],
                'Densidad_ρ': self.modelo[8],
                'ritmo_generacion_energia_ε': self.modelo[9],
                'Opacidad_κ': self.modelo[10]
            }

            df = pd.DataFrame(datos_csv)
            filename = name if name else f"modelo_M{self.Mtot}.csv"
            df.to_csv(filename, index=False)
            print(f"\n Archivo CSV guardado como: {filename}")

    def effective_temp(self):
        """
        Calculation of the effective temperature using the Stefan-Boltzmann law.

        Returns
        -------
        self.Teff
            Effective temperature [K].
        """
        sigma = 5.6704e-5
        self.Teff = ((self.Ltot * 1e33) / (4 * math.pi * sigma * (self.Rtot * 1e10) ** 2)) ** (1 / 4)
        return float(self.Teff)

    def HR(self):
        """
        Plotting the HR diagram using the Hipparcos catalogue with 10.000 stars.
        """
        file_path = 'asu (1).tsv'
        df = pd.read_csv(file_path, sep='|', comment='#', skipinitialspace=True)
        df = df.iloc[2:].reset_index(drop=True)

        df.columns = df.columns.str.strip()
        for col in ['Vmag', 'Plx', 'B-V']:
            df[col] = pd.to_numeric(df[col], errors='coerce')

        df = df.dropna(subset=['Vmag', 'Plx', 'B-V'])
        df = df[df['Plx'] > 0]

        # Absolute magnitude: M = V + 5*log10(parallax_in_seconds) + 5
        df['Mv'] = df['Vmag'] + 5 * np.log10(df['Plx']) - 10

        # Conversion from B-V to Effective Temperature (Ballesteros’ formula)
        bv = df['B-V']
        df['Teff'] = 4600 * (1/(0.92*bv + 1.7) + 1/(0.92*bv + 0.62))

        # Conversion from Mv to luminosity (L/Lsun) (approximate)
        # Using M_sun = 4.83
        df['L_Lsun'] = 10**(0.4 * (4.83 - df['Mv']))


        T_star = self.Teff  # K
        L_star = self.Ltot / 3.828

        plt.figure(figsize=(10, 8))

        plt.scatter(df['Teff'], df['L_Lsun'], s=1, color='gray', alpha=0.3, label='Hipparcos Data')

        plt.scatter(T_star, L_star, color='red', s=150, marker='*', edgecolors='black', label='Mi Modelo', zorder=1000)

        plt.xscale('log')
        plt.yscale('log')
        plt.gca().invert_xaxis()

        plt.xlabel(r'Temperatura Efectiva $T_{eff}$ (K)')
        plt.ylabel(r'Luminosidad ($L/L_{\odot}$)')
        plt.title('Diagrama HR: Modelo Estelar sobre Datos de Hipparcos')
        plt.legend()
        plt.grid(True, which='both', alpha=0.2)

        plt.show()

    def Spectral(self):
        """
        Determine the star's spectral type and is plotted in the same format as the HR diagram.
        """
        HRDiagram(self.Ltot,self.Rtot)

    def charts(self, x_axis='r', which=['r','P', 'T', 'L', 'M', 'rho','epsilon','kappa'], merge:bool =True, normalise:bool =True, grid: bool=True):
        """Generate the internal structure plots against mass or radius.

        Parameters
        ----------
        x_axis : str, optional
            Magnitude to use for x-axis.
        which : list, optional
            Magnitudes to use for the y-axis.
        merge: bool
            If True, merge the plots of the different magnitudes in the y-axis into a single figure.
        normalise : bool
            If True, normalize the magnitudes of the y-axis between 0 and 1.
        grid : bool
            If True, plot the grid."""
        datos={'r':self.modelo[2],'P':self.modelo[3], 'T':self.modelo[4], 'L':self.modelo[5], 'M':self.modelo[6],'parametro_rad':self.modelo[7], 'rho':self.rho,'epsilon':self.epsilon,'kappa':self.kappa}

        parameter_charts(datos,x_axis=x_axis,which=which,merge=merge,normalise=normalise,grid=grid)

    def border_point(self):
        """Calculates and returns the physical values at the boundary between the convective and radiative zones."""

        parametro_rad = self.par_rad

        try:
            # We are looking for the transition point where n ≤ 2.5
            idx = next(i for i, v in enumerate(parametro_rad)if v <= 2.5 and i > 0 and parametro_rad[i-1] > 2.5)

            r, P, T, L, M = self.modelo[2], self.modelo[3], self.modelo[4], self.modelo[5], self.modelo[6]

            valores = border(parametro_rad, idx, r, P, T, L, M)
            _,p_f, t_f, l_f, m_f, r_f = valores
            rho_f=rho_fun(self.mu,p_f,t_f)
            eps_f=gen_rate(t_f,self.X,self.Z)[0] * rho_f
            kappa_f=opacity(rho_f,t_f,self.X,self.Z)

            print("\n" + "="*40)
            print("DATOS EN LA BASE DE LA ZONA CONVECTIVA")
            print("="*40)
            print(f"Radio (r):       {r_f:12.4f} · 10¹⁰ cm")
            print(f"Presión (P):     {p_f:12.4f} · 10¹⁵ din/cm²")
            print(f"Temperatura (T): {t_f:12.4f} · 10⁷ K")
            print(f"Luminosidad (L): {l_f:12.4f} · 10³³ erg/s")
            print(f"Masa (M):        {m_f:12.4f} · 10³³ g")
            print(f"Densidad (ρ):     {rho_f:12.4f} · g/cm³")
            print(f"Ritmo de gen. de energía (ε):        {eps_f:12.4f} · erg/g/s")
            print(f"Opacidad (κ):        {kappa_f:12.4f} · cm²/g")
            print("-" * 40)
            print(f"Punto de corte: n = 2.5")
            print("="*40 + "\n")

            return valores

        except StopIteration:
            print("No se detectó una transición de zona radiativa a convectiva.")
            return None

    def TDplane(self):
        """
        Plot all the layers of the star on the density-temperature diagram.
        """
        if self.T is None:
            raise ValueError("No se ha calculado el modelo para la estrella.")
        T_array = np.array(self.T)
        rho_array = np.array(self.rho)
        T_scale = 1e7
        rho_scale = 1.0
        parametro_rad = self.par_rad
        t_f, rho_f = None, None

        try:
            # We are looking for the transition point where n ≤ 2.5
            idx = next(i for i, v in enumerate(parametro_rad)if v <= 2.5 and i > 0 and parametro_rad[i-1] > 2.5)

            r, P, T, L, M = self.modelo[2], self.modelo[3], self.modelo[4], self.modelo[5], self.modelo[6]

            valores = border(parametro_rad, idx, r, P, T, L, M)
            _,p_f, t_f, l_f, m_f, r_f = valores
            rho_f=rho_fun(self.mu,p_f,t_f)
            t_f = np.log10(t_f * T_scale)
            rho_f = np.log10(rho_f * rho_scale)

        except StopIteration:
            print("No se detectó una transición de zona radiativa a convectiva.")

        logT_star = np.log10(T_array * T_scale)
        log_rho_star = np.log10(rho_array * rho_scale)

        plt.rcParams['font.family'] = 'serif'
        fig, ax = plt.subplots(figsize=(9, 7))

        estilo_texto = {'fontsize': 11, 'fontweight': 'bold', 'va': 'center', 'ha': 'center'}

        ax.text(5.0, 8.5, r'III: Degeneración Relativista', color='white', **estilo_texto)
        #ax.text(5.0, 7.8, r'$P \propto \rho^{4/3}$', color='white', fontsize=10, ha='center')

        ax.text(4.5, 3.5, r'II: Degeneración NR', color='white', **estilo_texto)
        #ax.text(4.5, 2.8, r'$P \propto \rho^{5/3}$', color='white', fontsize=10, ha='center')

        ax.text(5.0, -2.5, r'I: Gas Ideal', color='black', **estilo_texto)
        #ax.text(5.0, -3.2, r'$P = \frac{\mathcal{R}}{\mu} \rho T$', color='black', fontsize=10, ha='center')

        ax.text(8.5, -6.0, r'IV: Presión de Radiación', color='black', **estilo_texto)
        #ax.text(8.5, -6.7, r'$P = \frac{1}{3} a T^4$', color='black', fontsize=10, ha='center')

        logT_axis = np.linspace(3.3, 10, 500)
        # Physical items:
        log_rho_I_II = 1.5 * logT_axis - 8.0   # Gas Ideal vs You. NR
        log_rho_I_IV = 3 * logT_axis - 23.0    # Ideal Gas vs Radiation
        log_rho_limit_rel = 6.0                 # Relativistic limit

        ax.fill_between(logT_axis, log_rho_limit_rel, 10.5, color='#4d4d4d')
        ax.fill_between(logT_axis, log_rho_I_II, 6.0, where=(6.0 > log_rho_I_II), color='#7f7f7f')
        ax.fill_between(logT_axis, log_rho_I_IV, log_rho_I_II, color='#d9d9d9')
        ax.fill_between(logT_axis, -10.5, log_rho_I_IV, color='#f2f2f2')

        ax.plot(logT_axis, log_rho_I_II, color='black', lw=1, alpha=0.5)
        ax.plot(logT_axis, log_rho_I_IV, color='black', lw=1, alpha=0.5)
        ax.axhline(6.0, color='black', lw=1, xmin=0, xmax=0.85, alpha=0.5)

        ax.plot(logT_star, log_rho_star, color='orange', lw=3, zorder=10, label='Perfil Estelar')

        ax.scatter(logT_star[-1], log_rho_star[-1], marker='*', s=250,
                   color='orange', edgecolors='black', zorder=11, label='Centro')

        ax.scatter(logT_star[0], log_rho_star[0], marker='s', s=70,
                   color='orange', edgecolors='black', zorder=11, label='Superficie')

        ax.scatter(t_f, rho_f, marker='o', s=70,
                   color='red', edgecolors='black', zorder=11, label='Frontera convectiva-radiativa')

        ax.set_title("Estructura Interna en el Diagrama T-Rho", fontsize=14)
        ax.set_xlabel(r"$\log_{10} T$ (K)", fontsize=12)
        ax.set_ylabel(r"$\log_{10} \rho$ (g cm$^{-3}$)", fontsize=12)
        ax.set_xlim(3.3, 10)
        ax.set_ylim(-10.5, 10.5)
        ax.legend(loc='upper right')
        ax.grid(True, which='both', linestyle=':', alpha=0.4)

        plt.tight_layout()
        plt.show()