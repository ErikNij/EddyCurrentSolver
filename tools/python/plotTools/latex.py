#!/usr/bin/python
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

from __future__ import nested_scopes, generators, division, absolute_import
from __future__ import with_statement, print_function, unicode_literals

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

__name__
__path__ = os.path.realpath(__file__)
__base__ = os.path.basename(__path__)
__dir__ = os.path.dirname(__path__)
__head__ = os.path.splitext(__base__)[0]

import locale as lc
import matplotlib as mpl
from matplotlib.backends.backend_pgf import FigureCanvasPgf
mpl.backend_bases.register_backend("pdf", FigureCanvasPgf)
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------- #
# --- Function -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def pdflatexify(
        fig_width=None, fig_height=None, dpi=None,
        fontsize=None, fontfamily=None, locale=None, preamble=None
    ):
    """Set up matplotlib"s RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inch
    fig_height : float,  optional, inch
    dpi : int, optional,
    fontsize : int, optional
    fontfamily : str, optional
    """

    inch = 25.4

    ## Adjust size and dpi
    if fig_width and fig_height:

        mpl.rcParams["figure.figsize"] = [fig_width,fig_height]

    if dpi:

        mpl.rcParams["figure.dpi"] = dpi
        mpl.rcParams["savefig.dpi"] = dpi

    # Adjust locale
    if locale:

        lc.setlocale(lc.LC_ALL, locale)

    else:

        localeSet = lc.getlocale()
        locale = localeSet[0] + "." + localeSet[1]

    mpl.rcParams["axes.formatter.use_locale"] = True

    # Adjust font sizes
    if fontsize:

        rc_font_sizes = {
            "font.size": fontsize,
            "axes.labelsize": fontsize,
            "axes.titlesize": fontsize,
            "legend.fontsize": fontsize,
            "xtick.labelsize": fontsize,
            "ytick.labelsize": fontsize
        }

        mpl.rcParams.update(rc_font_sizes)

    # Adjust font family (important for clabels!!!
    if fontfamily:

        mpl.rcParams["font.family"] = fontfamily

    # Adjust pdflatex settings
    if not preamble:

        preamble = [
            "\\usepackage[utf8x]{inputenc}",
            "\\usepackage[T1]{fontenc}",
            "\\usepackage{amsmath,amssymb,amsthm,amsfonts,mathrsfs}",
            "\\usepackage{mathtools}",
            "\\usepackage[locale=" + locale[3:5] + "]{siunitx}"]

    rc_pgf_with_pdflatex = {
        "backend": "pgf",
        "font.serif": [],
        "font.sans-serif": [],
        "font.monospace": [],
        "text.usetex": True,
        "text.latex.unicode": True,
        "pgf.texsystem": "pdflatex",
        "pgf.rcfonts": False,
        "pgf.preamble": preamble
    }

    mpl.rcParams.update(rc_pgf_with_pdflatex)


# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

