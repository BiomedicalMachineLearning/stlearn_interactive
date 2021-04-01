"""Purpose of this script is to create general forms that are programmable with
	particular input. Will impliment forms for subsetting the data and
	visualisation options in a general way so can be used with any
	SingleCellAnalysis dataset.
"""

import sys
from flask_wtf import FlaskForm
from wtforms import SelectMultipleField, SelectField, StringField, \
					IntegerField, BooleanField

def createSuperForm(elements, element_fields, element_values,
					validators):
	""" Creates a general form; goal is to create a fully programmable form \
	that essentially governs all the options the user will select.

	Args:
		elements (list<str>): Element names to be rendered on the page, in \
							  order of how they will appear on the page.

		element_fields (list<str>): The names of the fields to be rendered. \
								Each field is in same order as 'elements'. \
								Currently supported are: \
								'Title', 'SelectMultipleField', 'SelectField', \
								'StringField', 'Text', 'List'.

		element_values (list<object>): The information which will be put into \
									the field. Changes depending on field: \

									'Title' and 'Text': 'object' is a string
									containing the title which will be added as \
									a heading when rendered on the page.

									'SelectMultipleField' and 'SelectField':
									'object' is list of options to select from.

									'StringField':
									The example values to display within the \
									fields text area. The 'placeholder' option.

									'List':
									A list of objects which will be attached \
									to the form.

		validators (list<FunctionHandles>): A list of functions which take the \
						form as input, used to construct the form validator. \
						Form validator constructed by calling these \
						sequentially with form 'self' as input.

	Args:
		form (list<WTForm>): A WTForm which has attached as variable all the \
		fields mentioned, so then when rendered as input to
		'SuperDataDisplay.html' shows the form.
	"""

	class SuperForm(FlaskForm):
		""" A base form on which all of the fields will be added.
		"""

		def validate(self):
			# Return True if no validators attached #
			if type(self.validators) == type(None):
				return True

			# Only return True if all the validators return true #
			for function in self.validators:
				if not function(self):
					return False

			return True

	# Add the information #
	SuperForm.elements = elements
	SuperForm.element_fields = element_fields
	SuperForm.validators = validators

	multiSelectLeft = True # Places multi-select field to left, alternatives
							# if many multi-selects in row
	for i, element in enumerate( elements ):
		fieldName = element_fields[i]

		# Adding each element as the appropriate field to the form #
		if fieldName=='SelectMultipleField':
			setattr(SuperForm, element,
					SelectMultipleField(element, choices=element_values[i]))
			# The point of this number is to give an order for the attributes,
			# so that odd numbers get rendered to right of page, even numbers
			# left.
			setattr(SuperForm, element + '_number', int(multiSelectLeft))
			# inverts, so if left, goes right for the next multiSelectField
			multiSelectLeft = multiSelectLeft == False

		else:
			multiSelectLeft = True # Reset the MultiSelectField position

			if fieldName in ['Title', 'List']:
				setattr(SuperForm, element, element_values[i])

			elif fieldName == 'SelectField':
				setattr(SuperForm, element,
						SelectField(element, choices=element_values[i]))

			elif fieldName == 'StringField':
				setattr(SuperForm, element,
						StringField(element))
				setattr(SuperForm, element+'_placeholder', #Setting default
						element_values[i])

			elif fieldName == 'IntegerField':
				setattr(SuperForm, element,
						IntegerField(element))
				setattr(SuperForm, element + '_placeholder',  # Setting default
						element_values[i])

			elif fieldName == 'BooleanField':
				setattr(SuperForm, element,
						BooleanField(element))
				setattr(SuperForm, element + '_placeholder',  # Setting default
						element_values[i])

	return SuperForm

def getPreprocessForm():
	""" Gets the preprocessing form generated from the superform above.

	Args:
		adata (AnnData): AnnData containing the spatial RNA-seq data.
	Returns:
		FlaskForm: With attributes that allow for inputs that are related to
					pre-processing.
	"""
	elements = ['Spot Quality Control Filtering', # Title
				'Minimum genes per spot', 'Minimum counts per spot',
				'Gene Quality Control Filtering', # Title
				'Minimum spots per gene', 'Minimum counts per gene',
				'Normalisation, Log-transform, & Scaling', # Title
				'Normalize total', 'Log 1P', 'Scale'
				]
	element_fields = ['Title', 'IntegerField', 'IntegerField',
					  'Title', 'IntegerField', 'IntegerField',
					  'Title', 'BooleanField', 'BooleanField', 'BooleanField'
					  ]
	element_values = ['', 200, 300, '', 3, 5, '', True, True, True]
	return createSuperForm(elements, element_fields, element_values, None)


