

Guide to the X-ray Light Curve Generator, XGEN.

_________________________________________________
XGEN was written by Abygail R. Waggoner(1*) and L. Ilsedore Cleeves(1,2)

(1) University of Virginia, Department of Chemistry
(2) University of Virginia, Department of Astronomy


This model is published in Waggoner & Cleeves (2022), and any use of XGEN should cite this work.


_________________________________________________
Purpose of XGEN:
_________________________________________________
XGEN is capable of producing a light curve based on the power law distribution

dN/dE = beta * logE ** -alpha

where N is the number of flares that over the energy range dE,
        E     = the total energy of the flare
        beta  = a normalization constant dependent on average flare frequency
        logE  = the total integrated energy of the flaring event
        alpha = the power law index
The script xray_generator.py and Waggoner & Cleeves (2022) describe how the
model produces the light curve in further detail.

__________________________________________________
Included files
___________________________________________________
XGEN itself is a single python script, xray_generator.py, but this directory
includes two additional scripts, an observed energy distribution for TTauri stars,
and an schematic of the model. All python scripts are executable in Python 3,
and last updated for Python 3.7.11.

Included with model:
- xray_generator.py
  This script generates the light curves, as described above and in the script.
- plot_curve.py
  This script plots individual light curves simulated by xray_generator.py.
  Modeled flare peaks can also be plotted.
- run_stats.py
  This script determines the *observed* energy distribution and *observable* flare frequency
  of the generated light curve(s). This script is optimized for comparing the simulated
  light curve to observed flare statistics.
- wolk05_edistdata.csv
  This csv file contains the observed energy distribution of solar mass TTauri stars
  from Wolk et al. (2005). This csv file is read into run_stats.py to compare
  simulated and observed flare statistics.
- xgen_visualized.png
  This image demonstrates how xray_generator.py creates the light curve and
  how run_stats.py determines the energy distribution of observable flares.

________________________________________________
Additional Notes:
________________________________________________

- The default input parameters used for XGEN are found to yeild the best
  fit to the energy distribution and flare frequency observed for a TTauri star
  as discussed in Waggoner & Cleeves (2022)

- Flares of different energy values can occur simultaneously

- XGEN is currently written where all flares have a uniform rise and decay time,
  but is written such that a non-uniform rise and decay time can be introduced

- XGEN distinguishes between `observable' and `modeled' flares, as defined below:
  modeled: any and all flares modeled. This includes any nano or microflares that can
           occur within the set energy range
  observable: these are *only* flares that would be considered distinguishable if the flare
           were observed. The user defines the minimum flare energy considered observable,
           and any flares with Etot < Emin (observed) are considered non-observable.
  Additionally: only flares with a distinguishable flare peak in the
           integrated light curve are considered observable. The final
           flare frequency and energy distribution takes overlap of modeled flares
           into account when determining observable flares.       
  See run_stats.py and Waggoner & Cleeves (2022) for further details.

- xray_generator considers both a modeled (L_characteristics) and an observed
  characteristic luminosity (Lchar_obs). This distinction is necessary, because
  micro and nano-flaring events can raise the modeled characterstic luminosity.
  If you find that XGEN is yielding a higher baseline luminosity than desired,
  (i.e. deltaLXR always > 1), then decrease Lcharacteristic.

______________________________________________
Acknowledgements
______________________________________________
A.R.W. would like to thank Dana Anderson, Rachel Gross,
and Korash Assani for testing XGEN.

__________________________________________
end of readme.txt
_____________________________________________
