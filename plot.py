#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot.py

Makes a plot of a run
"""
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pickle
import math
from evolver.utilities import unpack
from evolver.model import Model
from math import pi
from matplotlib.backends.backend_pdf import PdfPages
from enum import Enum
from evolver.eoms import eoms_terms
from cycler import cycler 
from matplotlib import colors as mcolors
# === Seminar preset (simple y robusto) ===
plt.rcParams.update({
    "figure.dpi": 140,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,   # texto editable al exportar PDF
    "ps.fonttype": 42,
    "font.family": "serif",
    "font.size": 16,
    "axes.labelsize": 18,
    "legend.fontsize": 14,
    "axes.linewidth": 1.3,
    "xtick.major.width": 1.1,
    "ytick.major.width": 1.1,
    "xtick.major.size": 6,
    "ytick.major.size": 6,
})

plt.rcParams["axes.prop_cycle"] = cycler(color=[
    "#377EB8",  # azul brillante
    "#E41A1C",  # rojo fuerte
    "#4DAF4A",  # verde saturado
    "#984EA3",  # morado intenso
    "#FF7F00",  # naranja fuerte
    "#A65628",  # café/naranja oscuro
    "#F781BF",  # rosa saturado
])


####################################
# Deal with command line arguments #
####################################
parser = argparse.ArgumentParser(description="Plot data from a run")
parser.add_argument("filename", help="Base of the filename to read data in from")
parser.add_argument("outfilename", help="Filename to output to (include .pdf)")
args = parser.parse_args()


#################
# Load the data #
#################
# Parameters
model = Model.load(args.filename + ".params")
params = model.eomparams

# File 1
with open(args.filename + ".dat") as f:
    data = f.readlines()
results = np.array([list(map(float, line.split(", "))) for line in data]).transpose()
t = results[0]
(Ydm) = unpack(results[1:])

# File 2
with open(args.filename + ".dat2") as f:
    data2 = f.readlines()
results2 = np.array([list(map(float, line.split(", "))) for line in data2]).transpose()
(Rate3to2, RateXiXitoh2h2, RateXih2toXjXk, Rate3toh2, RateHtoInv, RateDMtoSM, Hubble, sent, neqh1, neqh2, neqDM) = results2[1:]

# File 3
with open(args.filename + ".quick", 'rb') as f:
    quickdata = pickle.load(f)

Omega_h2 = quickdata.get("Omega_h2", None)
Yinf_final = quickdata.get("Yinf", None)

    
    
######################
# Plotting Functions #
######################

class PlotStyle(Enum):
    LINEAR = 1
    LOG10 = 2
    LOG10_XY = 3

def create_cover_sheet(fig):
    fig.patch.set_alpha(0.0)

    left_lines = [
        rf"$m_X$ = {params.model.mx:g} GeV",
        rf"$m_{{h_1}}$ = {params.model.mh1:g} GeV",
        rf"$m_{{h_2}}$ = {params.model.mh2:g} GeV",
        rf"$g_X$ = {params.model.gX:g}",
        rf"$\lambda_H$ = {params.model.lambdaH:.3g}",
        rf"$\lambda_\Phi$ = {params.model.lambdaPhi:.3g}",
        rf"$\lambda_M$ = {params.model.lambdaM:.3g}",
        rf"$\beta$ = {params.model.beta:g}",
    ]


    right_lines = [
        rf"$Y(0)$ = {Ydm[0]:.5g}",
        rf"$Y_\mathrm{{eq}}(0)$ = {(neqDM[0]/sent[0]):.5g}",
        rf"$H(0)$ = {Hubble[0]:.2e}",
        rf"Rate$_{{3\to2}}(0)$ = {Rate3to2[0]:.2e}",
        rf"Rate$_{{3\to2}}(0)/H(0)$ = {Rate3to2[0]/Hubble[0]:.2e}",
        rf"Rate$_{{h\to\mathrm{{inv}}}}(0)$ = {RateHtoInv[0]:.2e}",
        rf"Rate$_{{h\to\mathrm{{inv}}}}(0)/H(0)$ = {RateHtoInv[0]/Hubble[0]:.2e}",
        rf"Rate$_{{\mathrm{{DM}}\to\mathrm{{SM}}}}(0)$ = {RateDMtoSM[0]:.2e}",
        rf"Rate$_{{\mathrm{{DM}}\to\mathrm{{SM}}}}(0)/H(0)$ = {RateDMtoSM[0]/Hubble[0]:.2e}",
        rf"$Y_\infty$ = {Yinf_final:.3e}",
        rf"$\Omega_\chi h^2$ = {Omega_h2:.3f}",
        ("Stiff solver ON" if model.parameters["stiff_solver"] else "Stiff solver OFF"),
        f"Filename: {model.basefilename}",
    ]

    fig.text(0.07, 0.95, "\n".join(left_lines),
             ha="left", va="top", family="monospace", fontsize=14, linespacing=1.35)
    fig.text(0.55, 0.95, "\n".join(right_lines),
             ha="left", va="top", family="monospace", fontsize=14, linespacing=1.35)

        
        
def create_cover_sheet2(canvas):
    # Create a plot on the canvas
    ax = canvas.add_subplot(1, 1, 1)

    # Add the text we want
    ax.text(0.05, 0.65, r'inflationstart='+str(quickdata['inflationstart']))

    # Hide the ticks (this is an empty plot!)
    ax.tick_params(
        axis='both',
        which='both',
        bottom=False,
        top=False,
        left=False,
        right=False,
        labelleft=False,
        labelbottom=False)
    
    
def make_pdf(pages, filename):
    # Create the PDF file
    print(f"Creating {filename}")
    pdf_pages = PdfPages(filename)

    # Write the cover sheet
    plt.rcParams["font.family"] = "serif"
    canvas = plt.figure(figsize=(8.0, 8.0), dpi=70)
    create_cover_sheet(canvas)
    pdf_pages.savefig(canvas)

    # Write the second cover sheet
    plt.rcParams["font.family"] = "serif"
    canvas2 = plt.figure(figsize=(8.0, 8.0), dpi=70)
    create_cover_sheet2(canvas2)
    pdf_pages.savefig(canvas2)

    # Make the dangerplot
    #create_danger_plot(pdf_pages)

    # Create the plots
    for idx, page in enumerate(pages):
        # Create the canvas
        canvas = plt.figure(figsize=(14.0, 14.0), dpi=100)

        # Sort out the page configuration
        numfigs = len(page)
        if numfigs in [1, 2]:
            rows = 2
            cols = 1
        elif numfigs in [3, 4]:
            rows = 2
            cols = 2
        else:
            print(f"Page {idx+1} has too many figures; only 4 will be produced")
            rows = 2
            cols = 2
            page = page[0:4]

        # Create the figures
        for fig, definition in enumerate(page):
            # Define the plotting location
            plt.subplot(rows, cols, fig + 1)

            # Specify when to use scientific notation
            plt.ticklabel_format(style='scientific', axis='y', scilimits=(-2, 2))

            # Determine the plotting function
            if definition['y_type'] == PlotStyle.LOG10:
                plotter = plt.semilogy
            elif definition['y_type'] == PlotStyle.LINEAR:
                plotter = plt.plot
            elif definition['y_type'] == PlotStyle.LOG10_XY:
                plotter = plt.loglog
            else:
                print(f"Bad plotting instruction on page {idx+1}, figure {fig+1}")

            # Do we have one data series, or a list?
            data = definition['y']
            if not isinstance(data, list):
                data = [data]

            
            # Create the plot
            labels = definition.get('labels', None)
            line_colors = definition.get('colors', None)
            
            for i, y_series in enumerate(data):
                y_data = np.real(y_series)
            
                plotter(
                    definition['x'], y_data,
                    label=(labels[i] if labels and i < len(labels) and labels[i] else None),
                    color=(line_colors[i] if line_colors and i < len(line_colors) else None),
                    linewidth=2.6,
                    alpha=0.95
                )

         


            if labels:
                plt.legend(frameon=True, framealpha=0.18, loc="best", fontsize=26)


            # Set the plot range
            if definition['x_range']:
                plt.xlim(*definition['x_range'])
            else:
                plt.xlim((definition['x'][0], definition['x'][-1]))
            if definition['y_range']:
                plt.ylim(*definition['y_range'])

            # Apply labels
            plt.xlabel(definition['x_label'],fontsize=30)
            plt.ylabel(definition['y_label'],fontsize=30)
            plt.tick_params(axis="both", which="major", labelsize=30)
            
            plt.title(r"$\bf{FORBIDDEN\ regime}$", fontsize=33, pad=18)

            
            # Call the extra plotting function if desired
            # Call the extra plotting function if desired
            if definition['extra_plotting']:
                definition['extra_plotting']()
            
            # --- Toques estéticos por subplot ---
            ax = plt.gca()
            ax.grid(True, which="major", alpha=0.55)  # cuadrícula suave
            ax.minorticks_on()
            ax.spines["top"].set_visible(False)       # limpia bordes
            ax.spines["right"].set_visible(False)


        # Save the page
        # Save the page
        plt.tight_layout()
        pdf_pages.savefig(canvas)
        canvas.savefig(
    f"{filename.replace('.pdf','')}_page{idx+1}.png",
    dpi=300,
    bbox_inches="tight"
)


    # Finish the file
    pdf_pages.close()
    print("Finished!")
    
    
    
    
def define_fig(x_data, y_data,
               x_label=r"$x \equiv m_{DM}/T$", y_label=None,
               x_range=None, y_range=None,
               y_type=PlotStyle.LINEAR,
               extra_plotting=None,
               labels=None,
               colors=None):  ##This option is for labels
    """Constructs data for a figure"""
    return {
        'x': x_data,
        'y': y_data,
        'x_label': x_label,
        'y_label': y_label,
        'x_range': x_range,
        'y_range': y_range,
        'y_type': y_type,
        'extra_plotting': extra_plotting,
        'labels': labels,
        'colors': colors
    }
    
    
    
    
    
def early(data, range=(0, 6)):
    """Restricts the plotting range in x to the given range"""
    return {**data, 'x_range': range}

def decorate_regime(ax, x_left, x_right, params):
    """
    Reproduce exactamente el estilo visual de la imagen:
    - tres regiones coloreadas
    - dos líneas verticales punteadas
    - tres cajas de texto
    - etiqueta x_f
    """

    xmin, xmax = ax.get_xlim()

    # --- regiones ---
    ax.axvspan(xmin, x_left, color="#52559A", alpha=0.95)
    ax.axvspan(x_left, x_right, color="#85A7C8", alpha=0.90)
    ax.axvspan(x_right, xmax, color="#B6D7CB", alpha=0.95)  

    # --- líneas verticales ---
    ax.axvline(x_left, linestyle="--", color="k", linewidth=1.6)
    ax.axvline(x_right, linestyle="--", color="k", linewidth=1.6)

    # --- cajas de texto ---
    box = dict(boxstyle="round", fc="white", ec="0.4", alpha=0.95)

    ax.text(0.08, 0.57, r"$\Gamma \gg H$", transform=ax.transAxes, bbox=box, fontsize=26)
    ax.text(0.3, 0.57, r"$\Gamma \gtrsim H$", transform=ax.transAxes, bbox=box, fontsize=26)
    ax.text(0.70, 0.57, r"$\Gamma < H$", transform=ax.transAxes, bbox=box, fontsize=26)

    # --- etiqueta x_f ---
    ax.text(x_right-1, 0.02, r"$x_f$", transform=ax.get_xaxis_transform(),
            ha="center", rotation=90, fontsize=28)

        # --- caja con parámetros (arriba-izquierda) ---
    pbox = dict(
    boxstyle="round,pad=0.45,rounding_size=0.2",
    fc="white",
    ec="0.35",
    lw=1.5,
    alpha=0.97
)

    txt = (
        rf"$m_X = {params.model.mx:g}\ \mathrm{{GeV}}$" "\n"
        rf"$m_{{h_2}} = {params.model.mh2:g}\ \mathrm{{GeV}}$" "\n"
        rf"$g_X = {params.model.gX:g}$" "\n"
        rf"$\beta = {params.model.beta:g}$"
    )

    ax.text(0.70, 0.36, txt, transform=ax.transAxes,
            ha="left", va="top", bbox=pbox, fontsize=30)
    
    
####################
# Plot Definitions #
####################

# Note that it is relatively quick to make all these definitions
# The slow part is the plotting of whatever is actually included in the PDF
# So, it is convenient to define everything we could ever want to plot here!

#x-axis
#t = results[0]
#(Ydm) = unpack(results[1:]

#Equilibrium abundance
Yeqphi = neqDM/sent

#Abundance plots
YPlot = define_fig(
    x_data=t,
    y_data=[Ydm, Yeqphi],
    y_label=r'$Y_{\chi}$',
    y_type=PlotStyle.LOG10_XY,
    y_range=(1e-14, 1e-2),
    labels=[r"$Y_{\chi}$", r"$Y_{\chi}^{0}$"],
    colors=["#1A1A1A", "#FFFFFF"]
)

YPlot['extra_plotting'] = lambda: decorate_regime(
    plt.gca(),
    x_left=3,      # primera línea vertical
    x_right=20,
    params=params# segunda línea vertical
)
                    
#Rate plots
RatesPlot = define_fig(x_data=t,
                        y_data=[Hubble,
                                RateDMtoSM,
                                RateXiXitoh2h2,    
                                RateXih2toXjXk,    
                                Rate3to2,          
                                Rate3toh2 ],
                        y_label=r'Rates per-DM particle',
                        y_type=PlotStyle.LOG10,
                        labels=[
    r"$H$",
    r"$\Gamma_{X_i X_i\rightarrow \mathrm{SM\,SM}}$",
    r"$\Gamma_{X_i X_i \rightarrow h_2 h_2}$",
    r"$\Gamma_{X_i h_2\rightarrow X_j X_k}$",
    r"$\Gamma_{3\rightarrow2}$",
    r"$\Gamma_{3\rightarrow h_2}$"
],
                       x_range=(10, 50),
                       y_range=(1e-30, 1e-10) )

# === Boltzmann Terms
term_names = ["XiXi_to_h2h2_raw", "Xih2_to_XjXk_raw", "3to2_to_XX_raw", "3to2_to_Xh2_raw", "DM_to_SM_raw"]
term_curves = {k: np.zeros_like(t) for k in term_names}

for i in range(len(t)):
    tr = eoms_terms(Ydm[i], params, t[i])
    for k in term_names:
        term_curves[k][i] = tr[k]

neqH_DM = neqDM * Hubble

TermsRawPlot = define_fig(
    x_data=t,
    y_data=[
        neqH_DM,
        term_curves["XiXi_to_h2h2_raw"],
        term_curves["Xih2_to_XjXk_raw"],
        term_curves["3to2_to_XX_raw"],
        term_curves["3to2_to_Xh2_raw"],
        term_curves["DM_to_SM_raw"]
    ],
    y_label=r'RHS termss',
    y_type=PlotStyle.LINEAR,  
    labels=[
        r'$n_{\rm eq}^{\rm DM}\,H$',      
        r'$XX\leftrightarrow h_2 h_2$',
        r'$X_i\,h_2\to X_j X_k$',
        r'$3\to 2:\,XXX\to XX$',
        r'$3\to 2:\,XXX\to X\,h_2$',
        r'$XX\leftrightarrow \mathrm{SM\,SM}$'
    ],
    x_range=(-10, 50)   
)
TermsRawPlot['extra_plotting'] = lambda: plt.yscale('symlog', linthresh=1e-30)
                       
                        



###########################
# PDF Layout and Creation #
###########################

# Lay out the figures in pages
# We recommend commenting out pages that you don't want, rather than deleting them
pages = [
    [YPlot],
    [RatesPlot],
    [TermsRawPlot]
]

# Construct the PDF
make_pdf(pages, args.outfilename)
