from matplotlib.ticker import FixedLocator, FixedFormatter
import math
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


def normalize(array):
    """Normalises an array to values between 0 and 1"""
    return (array - array.min()) / (array.max() - array.min())

def parameter_charts(data_dict, x_axis='r', which=['P', 'T', 'L', 'M'], merge=True, normalise=True, grid=True):
    """
    Plots of stellar parameters.

    data_dict : dict -> {'r': array, 'P': array, ...}
    x_axis    : str  -> Magnitude for the X-axis
    which     : list -> Magnitudes for the Y-axis
    merge     : bool -> If True, plot everything on a single graph. If False, subplots.
    normalise : bool -> If True, scale 0-1.
    grid : bool -> If true, grid
    """

    unidades = {
        'r': r'Radio ($10^{10}$ cm)',
        'P': r'Presión ($10^{15}$ din cm$^{-2}$)',
        'T': r'Temperatura ($10^{7}$ K)',
        'L': r'Luminosidad ($10^{33}$ erg s$^{-1}$)',
        'M': r'Masa ($10^{33}$ g)',
        'rho': r'Densidad ($\text{g cm}^{-3}$)',
        'epsilon': r'Generación energía ($\text{erg g}^{-1} \text{s}^{-1}$)',
        'kappa': r'Opacidad ($\text{cm}^2 \text{g}^{-1}$)'}

    x_data = np.array(data_dict[x_axis])
    x_data = normalize(x_data)
    parametro_rad=np.array(data_dict['parametro_rad'])

    v_fron = None
    if parametro_rad is not None:
        try:
            idx = next(i for i, v in enumerate(parametro_rad) if v == 0 and any(x > 0 for x in parametro_rad[:i]))
            v_fron = x_data[idx]
        except (StopIteration, IndexError):
            pass

    num_plots = len(which)
    if merge:
        fig, ax = plt.subplots(figsize=(9, 6))
        axs = [ax] * num_plots
    else:
        fig, axs = plt.subplots(num_plots, 1, figsize=(8, 3 * num_plots), sharex=True)
        if num_plots == 1: axs = [axs]
    frontera_dibujada = False
    for i, key in enumerate(which):
        y_data = np.array(data_dict[key])

        if normalise:
            y_data = normalize(y_data)

        label_text = unidades.get(key, key)
        axs[i].plot(x_data, y_data, label=label_text, linewidth=1.5)

        if v_fron is not None:
            if not frontera_dibujada:
                axs[i].axvline(v_fron, color='black', linestyle='--', alpha=0.6, label='Frontera convectivo-radiativa')
                if merge: frontera_dibujada = True
            else:
                axs[i].axvline(v_fron, color='black', linestyle='--', alpha=0.6)

        if not merge:
            #axs[i].set_ylabel('Valor normalizado' if normalize else key)
            axs[i].legend(loc='upper right', fontsize='x-small')

    if merge:
        #axs[0].set_ylabel('Valores Normalizados' if normalize else 'Valores del Modelo')
        axs[0].legend(loc='best', fontsize='small', frameon=True)
        #axs[0].set_title(f'Perfil estelar respecto a {x_axis}')

    axs[-1].set_xlabel(f'{unidades.get(x_axis, x_axis)} (normalizado)')
    axs[-1].set_xlim(0, 1)
    plt.tight_layout()
    if grid:
        plt.grid()
    plt.show()



def HRDiagram(L: float, R: float, debug_spectral: bool = True, debug_HR: bool = True):
        """
        Determine the star’s effective temperature from its luminosity and radius, and plot its position on the HR diagram.

        Parameters
        ----------
        L : float
            Brightness [10^33 erg/s].
        R : float
            Radius [10^10 cm].
        debug_spectral : bool
            Print spectral star's type information.
        debug_HR : bool
            Print HR diagram position of the star.

        Returns
        -------
        Teff : float
            Spectral star's effective temperature [K].
        """
        sigma = 5.6704e-5
        L_sol_unidades_tfg = 3.828
        Teff = ((L * 1e33) / (4 * math.pi * sigma * (R * 1e10) ** 2)) ** (1 / 4)
        if debug_HR:
            log_T_star = np.log10(Teff)
            log_L_star = np.log10(L / L_sol_unidades_tfg)

            fig, ax = plt.subplots(figsize=(11, 8))

            zonas = [
                (4.6, 4.5, 0.4, 1.2, -30, 'blue', 'O', 'Tipo O (Azul, >30k K)'),
                (4.2, 2.5, 0.5, 1.5, -35, 'deepskyblue', 'B', 'Tipo B (Blanco-Azul, 10k-30k K)'),
                (3.9, 1.1, 0.3, 0.8, -40, 'silver', 'A', 'Tipo A (Blanco, 7.5k-10k K)'),
                (3.8, 0.4, 0.2, 0.6, -40, 'yellow', 'F', 'Tipo F (Blanco-Amarillo, 6k-7.5k K)'),
                (3.75, 0.0, 0.2, 0.5, -40, 'orange', 'G', 'Tipo G (Amarillo, 5k-6k K)'),
                (3.6, -0.6, 0.3, 0.7, -40, 'orangered', 'K/M', 'Tipo K/M (Naranja/Rojo, <5k K)')
            ]

            for x, y, w, h, ang, col, short_lab, full_lab in zonas:
                ell = patches.Ellipse((x, y), w, h, angle=ang, color=col, alpha=0.2, label=full_lab)
                ax.add_patch(ell)
                ax.text(x, y, short_lab, fontsize=12, ha='center', va='center', fontweight='bold')

            ax.scatter(log_T_star, log_L_star, color='gold', marker='*', s=400,
                       edgecolor='black', zorder=10, label='Mi estrella')

            ax.set_xlim(4.9, 3.3)
            ticks_temp = [40000, 30000, 20000, 10000, 8000, 6000, 4000, 3000]
            ticks_log = [np.log10(t) for t in ticks_temp]
            ax.xaxis.set_major_locator(FixedLocator(ticks_log))
            ax.xaxis.set_major_formatter(FixedFormatter([f"{t}" for t in ticks_temp]))

            ax.set_ylim(-1.5, 5.5)
            ax.set_xlabel(r'Temperatura Efectiva $T_{eff}$ [K] (Escala Log)', fontsize=12)
            ax.set_ylabel(r'$\log_{10}(L/L_{\odot})$', fontsize=12)
            ax.set_title('Diagrama HR: Clasificación Espectral de la estrella', fontsize=14, pad=20)

            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="Clasificación", frameon=True)

            ax.grid(True, which='both', linestyle=':', alpha=0.6)
            plt.tight_layout()
            plt.show()

        if debug_spectral:
            print(f"--- RESULTADOS ---")
            print(f"Temperatura efectiva calculada: {Teff:.2f} K")
            if Teff >= 30000:
                print("Estrella clase O (Azul)")
            elif 10000 <= Teff < 30000:
                print("Estrella clase B (Blanco-Azul)")
            elif 7500 <= Teff < 10000:
                print("Estrella clase A (Blanco)")
            elif 6000 <= Teff < 7500:
                print("Estrella clase F (Blanco-Amarillo)")
            elif 5000 <= Teff < 6000:
                print("Estrella clase G (Amarillo)")
            elif 3700 <= Teff < 5000:
                print("Estrella clase K (Naranja)")
            elif 2000 <= Teff < 3700:
                print("Estrella clase M (Rojo)")
            else:
                print("Tipo desconocido (Fuera de rango estelar estándar)")