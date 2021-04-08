""" This is more a general views focussed on defining functions which are \
	called by other views for specify pages. This way different pages can be \
	used to display different data, but in a consistent way.
"""

import sys
from flask import flash
from source.forms import forms

from source.forms.utils import flash_errors
from flask import render_template

import scanpy as sc

# Creating the forms using a class generator #
PreprocessForm = forms.getPreprocessForm()

def getVal(preprocess_form, element):
	return getattr(preprocess_form, element).data

def run_preprocessing(request, adata, step_log):
	""" Performs the scanpy pre-processing steps based on the inputted data.
	"""

	form = PreprocessForm(request.form)

	if not form.validate_on_submit():
		flash_errors(form)

	elif type(adata) == type(None):
		flash("Need to load data first!")

	else:
		# Logging params used #
		form_elements = form.elements
		form_fields = form.element_fields
		for i, element in enumerate(form_elements):
			if form_fields[i] != 'Title':
				data = getVal(form, element)
				step_log['preprocessed_params'][element] = data

		# QC filtering #
		sc.pp.filter_cells(adata,
						   min_genes=getVal(form, 'Minimum genes per spot'))
		sc.pp.filter_cells(adata,
						   min_counts=getVal(form, 'Minimum counts per spot'))
		sc.pp.filter_genes(adata,
						   min_cells=getVal(form, 'Minimum spots per gene'))
		sc.pp.filter_genes(adata,
						   min_counts=getVal(form, 'Minimum counts per gene'))

		# Pre-processing #
		if getVal(form, 'Normalize total'):
			sc.pp.normalize_total(adata, target_sum=1e4)
		if getVal(form, 'Log 1P'):
			sc.pp.log1p(adata)
		if getVal(form, 'Scale'):
			sc.pp.scale(adata, max_value=10)

		# Setting pre-process to true #
		step_log['preprocessed'][0] = True

	if step_log['preprocessed'][0]:
		flash("Preprocessing completed!")

	updated_page = render_template("preprocessing.html",
								   title="Preprocessing",
									preprocess_form=form,
								   preprocessed=step_log['preprocessed'][0],
                                   step_log=step_log)

	return updated_page

def run_cci(request, adata, step_log):
	""" Performs CCI analysis & alters the CCI form depending on the optional
		aspects of the analysis chosen by the user as specified in step_log.
	"""

	CCIForm = forms.getCCIForm()
	form = CCIForm(request.form)

	if not form.validate_on_submit():
		flash_errors(form)

	elif type(adata) == type(None):
		flash("Need to load data first!")

	elif type(step_log['cci_het']) != type(None):
		""" If we have been through the steps of:
			1. Choosing with/without cell heterogeneity information.
			2. Choosing with/without permutation testing.
			3. Filling out information.
			We are now ready to take the form input & perform CCI analysis.
		"""
		# TODO impliment the CCI analyses, saving plots to temp_plots/

		print("We are using the cell heterogeneity!", file=sys.stdout)

	updated_page = render_template("cci.html",
								   title="Cell-Cell Interaction",
								   cci_form=form
								   )

	return updated_page



