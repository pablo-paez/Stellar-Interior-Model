
# 🌟 StarModel: Stellar Internal Structure Simulation

> A Python-based numerical engine designed to model the physical and thermodynamic structure of stars based on initial mass and chemical composition parameters.

---

Python module for solving the stellar-interior equations.

 - Author: **Pablo Páez Ramos**
 - Last update: May, 2026

It provides:

* A Star object to model massive stars.
* Two different options to optimize modeling.
* A 10.000 stars data set for Hertzsprung-Russell diagram representation retrived from The Hipparcos and Tycho Catalogues ESA SP-1200, 1997, available in the catalog service VizieR (CDS), cat. I/239, http://vizier.cds.unistra.fr/.

## 🚀 Installation

---

### Prerequisites
Python 3.8 or higher is recommended. It is advised to use a virtual environment.

Clone the repository and install the package:

```bash
git clone https://github.com/pablo-paez/Stellar-Interior-Model.git
cd Stellar-Interior-Model
pip install .
```
## Star class

---

It represents a star whose mass and chemical composition allow you to calculate a variety of options. 
The first thing you should do after creating an object is to calculate the layered model so you can use the other functions.

### Atributes

---

The `StarModel` object takes the following initial parameters:
* `Mtot` (Float) : Total mass of the star [10^33 g]
*  `X` (Float) : Fraction of hydrogen
*  `Y` (Float) : Fraction of helium
* `Ltot` (Float) : Total brightness of the star [10^33 erg/s] 
* `Rtot` (Float) : Total radium of the star [10^10 cm]
* `Tc` (Float) : Central temperature [10^7 K]
*  `unidades_solares` (Bool) : If true, the units of the mass, brightness and radius introduced will be transformed from solar units.

Once the object with these attributes is initialized, it calculates directly:
* `Z=1-X-Y` : Fraction of metalicity
*  `mu=1/(2*X+3*Y/4+Z/2)` : Average molecular weight

In addition, there are other attributes that start with the value `None` and will be updated when methods are called. These are:
* `modelo` : Save all the values for each layer of the star.
*  `Teff` : Effective temperature of the star
*  `error_rel` : Percentage of the star's relative error.
*  `tipo_energia` : Save the power generation process (None means --, 1 means PP and 0 means CN).
*  `fase_transp` : It retains the energy transfer mechanism (None==ˆˆˆˆˆˆ, 0==RADIAT, 1==CONVEC, 2==CENTRO and 3==INICIO).
* `r` : radium of each layer [10^10 cm]
* `P` : pressure of each layer [10^15 din/cm^2]
* `T` : temperature of each layer [10^7 K]
* `L` : brightness of each layer [10^33 erg/s]
* `M` : mass of each layer [10^33 g]
*  `par_rad` : radiative parameter value
*  `rho` : density [g/cm^3]
*  `epsilon` : rate of power generation [erg/g/s]
*  `kappa` : opacity [cm^2/g]
*  `_optimizado` : Internal attribute to check whether the optimisation of the radio and brightness parameters has been carried out.

### Methods

---

There are several methods for each object:
* `parameters`

  Returns the parameters defined.
* `redefine`

  Redefines the instance's physical parameters (the parameters that are not introduced remains the same). Arguments
  * `Mtot` (Float) : Total mass [10^33 g]
  * `X` (Float) : Hydrogen fraction
  * `Y` (Float) : Helium fraction
  * `Ltot` (Float) : Total brightness of the star [10^33 erg/s] 
  * `Rtot` (Float) : Total radium of the star [10^10 cm]
  * `Tc` (Float) : Central temperature [10^7 K]
  *  `unidades_solares` (Bool) : If true, the units of the mass, brightness and radius introduced will be transformed from solar units.
* `error` 

  Returns the total relative error (E) as a percentage, calculated as the sum of the squares of the differences between the parameters at the boundary between the convective and radiative zones.

* `calculate`
  
  Calculation of the complete stellar model based on the star’s initial parameters. This includes an option (optimise True/False) to decide whether to perform a search for the most optimal radius and luminosity parameters (in both cases, this search is carried out for the core temperature). Arguments:
  * `optimizar` (Bool) :  If false, only the core temperature is optimised (study the possible inputs of the function seccion_aurea); if true, a search is also carried out for the most optimal radius-luminosity pair (study the possible inputs of the function modelo_completo).
  * `mallado_final` (Bool) : If optimizar is true, it displays the final radio-luminance mesh showing the optimal value.
  * `all_mallados` (Bool) : If true, displays the meshes for all the refinements carried out.
  * `parametros_optimos` (Bool) : If true, displays the parameters of the optimization searched.
  * `magnitudes_extra` (Bool) : If true, the function calculates the density, opacity and epsilon.
* `layers`

  Displays the model layers and, if desired, exports them to CSV. Arguments:
  * `to_csv` (Bool) :  If True, save the data to a .csv file.
  * `name` (str) : CSV file name
* `effective_temp`

  Calculation of the effective temperature using the Stefan-Boltzmann law.
* `HR`

  Plotting the HR diagram using the Hipparcos catalogue with 10.000 stars in asu(1).tsv.
* `Spectral`

  Determine the star's spectral type and is plotted in the same format as the HR diagram.
* `charts`

  Generate the internal structure plots against mass or radius. Arguments:
  * `x_axis` (str) : Magnitude to use for x-axis
  * `which` (list) : List of magnitudes using for the y-axis
  * `merge` (Bool) : If True, merge the plots of the different magnitudes in the y-axis into a single figure.
  * `normalise` (Bool) : If True, normalize the magnitudes of the y-axis between 0 and 1.
  * `grid` (Bool) : If true, plot the grid.
* `border_point`

  Calculates and returns the physical values at the boundary between the convective and radiative zones.
* `TDplane`

  Plot all the layers of the star on the density-temperature diagram.
  
## Units

---

In order to properly estimate the variables of the star, the unit system adopted for internal calculations of the model varies with respect to CGS. However, both input and output values of the model can be expressed in any unit system using `Quantity` objects from astropy.

    radius (r)                         ->   1e10 cm
    pressure (P)                       ->   1e15 dyn cm^-2
    temperature (T)                    ->   1e7 K
    mass (m)                           ->   1e33 g
    luminosity (l)                     ->   1e33 erg s^-1
    density (rho)                      ->   1 g cm^-3
    energy generation rate (epsilon)   ->   1 erg g^-1 s^-1
    opacity (kappa)                    ->   1 cm^2 g^-1


## 📂 Repository Structure

---

* `Star_Model/`: Core package containing the logic and the `StarModel` class.
* `Example_Usage.ipynb`: Interactive Jupyter Notebook with simulation examples and plots.
* `setup.py`: Configuration file for package installation.
* `requirements.txt`: List of required Python dependencies.


## Usage example

---

For an example of how it works, go to the file example.ipynb

## License

---

This repository is licensed under the MIT License. See the LICENSE file for further details.
