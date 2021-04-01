# Pre-processing Form Design Notes
Flow of data:

    app.py/preprocessing() 
       -> source/forms/views.py/run_preprocessing(request, adata, step_log)
          -> source/forms/forms.py/getPreprocessForm()
          -> templates/superform.html

Notes:    
                  
    * step_log defined in app.py, keeps track whether pre-processing was run.
        -> If not run, then run_preprocessing shows just the form.
        -> If has run, shows banner that preprocessing complete.
    
    * If attempt to run_preprocessing() when adata not yet loaded, shows banner
        indicating need to upload the data first. 
        
    * source/forms/forms.py/getPreprocessForm() generates a WTForm class using a 
        general WTForm generator, as defined in:
        source/forms/forms.py/createSuperForm()
        
    * templates/superform.html renders a general WTForm that was created using 
        the source/forms/forms.py/createSuperForm() function, thereby allowing
        easy generation of new forms if need to add extra information. 
    
