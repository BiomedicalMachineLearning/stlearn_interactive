# -*- coding: utf-8 -*-
"""Helper utilities and decorators."""
from flask import flash
import matplotlib.pyplot as plt


def flash_errors(form, category="warning"):
    """Flash all errors for a form."""
    for field, errors in form.errors.items():
        for error in errors:
            flash(getattr(form, field).label.text + " - " + error + ", category")


def savePlot(plot_name, dpi=100, tight_layout=True, folder="temp_plots/"):
    """Deals with the current matplotlib.pyplot."""

    if tight_layout:
        plt.tight_layout()

    plt.savefig(folder + plot_name, dpi=dpi, format=plot_name.split(".")[-1])
