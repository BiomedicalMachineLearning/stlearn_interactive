""" This is more a general views focussed on defining functions which are \
	called by other views for specify pages. This way different pages can be \
	used to display different data, but in a consistent way.
"""

import sys
import numpy
from flask import flash
from source.forms import forms

from source.forms.utils import flash_errors, savePlot
import source.forms.view_helpers as vhs
import traceback

from flask import render_template

import scanpy as sc
import stlearn as st

# Creating the forms using a class generator #
PreprocessForm = forms.getPreprocessForm()
CCIForm = forms.getCCIForm()

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
		step_log['preprocessed_params'] = vhs.getData(form)
		print(step_log['preprocessed_params'], file=sys.stdout)

		# QC filtering #
		sc.pp.filter_cells(adata,
						   min_genes=vhs.getVal(form, 'Minimum genes per spot'))
		sc.pp.filter_cells(adata,
						 min_counts=vhs.getVal(form, 'Minimum counts per spot'))
		sc.pp.filter_genes(adata,
						   min_cells=vhs.getVal(form, 'Minimum spots per gene'))
		sc.pp.filter_genes(adata,
						 min_counts=vhs.getVal(form, 'Minimum counts per gene'))

		# Pre-processing #
		if vhs.getVal(form, 'Normalize total'):
			sc.pp.normalize_total(adata, target_sum=1e4)
		if vhs.getVal(form, 'Log 1P'):
			sc.pp.log1p(adata)
		if vhs.getVal(form, 'Scale'):
			sc.pp.scale(adata, max_value=10)

		# Setting pre-process to true #
		step_log['preprocessed'][0] = True

	if step_log['preprocessed'][0]:
		flash("Preprocessing completed!")

	updated_page = render_template("preprocessing.html",
								   title="Preprocessing",
									preprocess_form=form,
								   flash_bool=step_log['preprocessed'][0],
								   step_log=step_log)

	return updated_page

def run_cci(request, adata, step_log):
	""" Performs CCI analysis & alters the CCI form depending on the optional
		aspects of the analysis chosen by the user as specified in step_log.
	"""

	form = CCIForm(request.form)

	step_log['cci_params'] = vhs.getData(form)
	print(step_log['cci_params'], file=sys.stdout)
	elements = numpy.array(list(step_log['cci_params'].keys()))
	#order: cell_het file, neighbour_dist, L-R pairs, permutations
	element_values = numpy.array(list(step_log['cci_params'].values()))

	cell_het = type(element_values[0]) != type(None)
	lrs, msg = vhs.getLR(element_values[2], adata.var_names)
	print(lrs, file=sys.stdout)

	if not form.validate_on_submit():
		flash_errors(form)

	elif type(adata) == type(None):
		flash("Need to load data first!")

	elif msg != '':
		flash(msg)

	else:
		""" If we have been through the steps of:
			1. Choosing with/without cell heterogeneity information.
			2. Choosing with/without permutation testing.
			3. Filling out information.
			We are now ready to take the form input & perform CCI analysis.
		"""
		try:
			adata.uns['lr'] = lrs
			st.tl.cci.lr(adata=adata, distance=element_values[1])
			st.pl.het_plot(adata, use_het='cci_lr', image_alpha=0.7)
			savePlot('cell_lr.png')

			dist, n_pairs = element_values[1], element_values[-1]
			if cell_het: # Conduct with cell heterogeneity information #
				print("We are using the cell heterogeneity!", file=sys.stdout)
				# Adding the label transfer information #
				st.add.labels(adata, element_values[0], sep='\t')
				st.pl.cluster_plot(adata, use_label="predictions")
				savePlot('label_transfer.png') # saves to temp_pots

				# Calculating cell heterogeneity #
				st.tl.cci.het.count(adata, distance=dist,
									use_label='label_transfer')
				st.pl.het_plot(adata, use_het='cci_het')
				savePlot('cell_het.png')

				# Merging with the lR values #
				st.tl.cci.merge(adata, use_lr='cci_lr', use_het='cci_het')
				st.pl.het_plot(adata, use_het='merged', cell_alpha=0.7)
				savePlot('merged.png')

			if n_pairs != 0: # Permutation testing #
				st.tl.cci.permutation(adata,
									  use_het='cci_het' if cell_het else None,
												 n_pairs=n_pairs, distance=dist)
				st.pl.het_plot(adata, cell_alpha=0.7,
						 use_het='merged_pvalues' if cell_het else "lr_pvalues")
				savePlot('cci-log10pvalues.png')

				st.pl.het_plot(adata, cell_alpha=0.7,
								use_het='merged_sign'if cell_het else "lr_sign")
				savePlot('cci-sig-log10pvalues.png')

			step_log['cci'][0] = True

			flash('CCI analysis completed!')

		except Exception as msg:
			traceback.print_exc(file=sys.stdout)
			flash('Analysis ERROR: '+str(msg))
			print(msg)

	updated_page = render_template("cci.html",
								   title="Cell-Cell Interaction",
								   cci_form=form,
								   flash_bool=True,
								   step_log=step_log
								   )

	return updated_page



