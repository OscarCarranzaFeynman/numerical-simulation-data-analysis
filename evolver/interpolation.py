# interpolation.py
# Interpolación del ancho total de h2 (hadrons + leptones)
# a partir de la curva digitalizada en plot-data.csv

import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d

# Ruta al csv en la MISMA carpeta que este archivo
_datafile = Path(__file__).with_name("plot-data.csv")

# Carga de datos desde plot-data.csv
_data = np.genfromtxt(_datafile, delimiter=",", names=True)

# Tomamos las dos primeras columnas: x = mh2 [GeV], y = Gamma_tot [GeV]
_colnames = _data.dtype.names
_mh2_raw   = _data[_colnames[0]]
_Gamma_raw = _data[_colnames[1]]

_mask   = _Gamma_raw > 0.0
_mh2    = _mh2_raw[_mask]
_Gamma  = _Gamma_raw[_mask]

# Interpolación en log10(Gamma) vs mh2 (robusta a mh2 repetidos / desordenados)
_logGamma = np.log10(_Gamma)

# 1) ordenar por mh2
order = np.argsort(_mh2)
_mh2 = _mh2[order]
_logGamma = _logGamma[order]

# 2) colapsar mh2 duplicados (promediando logGamma en cada mh2 repetido)
mh2_unique, inv, counts = np.unique(_mh2, return_inverse=True, return_counts=True)
logGamma_sum = np.bincount(inv, weights=_logGamma)
logGamma_avg = logGamma_sum / counts

_mh2 = mh2_unique
_logGamma = logGamma_avg

_logGamma_interp = interp1d(
    _mh2,
    _logGamma,
    kind="linear",
    fill_value="extrapolate"
)


def Gamma_toth2(mh2):
    """
    Ancho total Gamma_tot(mh2) en GeV,
    obtenido de la curva digitalizada (hadrons + leptones).
    mh2 puede ser escalar o array.
    """
    mh2_arr = np.asarray(mh2, dtype=float)
    return 10.0**_logGamma_interp(mh2_arr)

__all__ = ["Gamma_toth2"]
