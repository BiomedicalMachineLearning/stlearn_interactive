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

# Creating the preprocessing form using a class generator #
PreprocessForm = forms.getPreprocessForm()

def getVal(preprocess_form, element):
	return getattr(preprocess_form, element).data

def run_preprocessing(request, adata, step_log):
	""" Performs the scanpy pre-processing steps based on the inputted data.
	"""

	form = PreprocessForm(request.form)

	print(adata, step_log, file=sys.stdout)

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
				print(data, file=sys.stdout)
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
		step_log['preprocessed'] = True

	if step_log['preprocessed']:
		flash("Preprocessing completed!")

	updated_page = render_template("superform.html",
								   title="Preprocessing",
									superForm=form,
								   preprocessed=True,#step_log['preprocessed']
									)

	return updated_page

