from flask import Flask, render_template, request, flash, url_for, redirect
from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.resources import INLINE
from werkzeug.utils import secure_filename

import os
import stlearn

import asyncio
from bokeh.server.server import BaseServer
from bokeh.server.tornado import BokehTornado
from tornado.httpserver import HTTPServer
from tornado.ioloop import IOLoop
from bokeh.application import Application
from bokeh.application.handlers import FunctionHandler
from bokeh.server.server import Server
from bokeh.embed import server_document

from bokeh.layouts import column, row

app = Flask(__name__)
app.secret_key = b'_5#y2L"F4Q8z\n\xec]/'

UPLOAD_FOLDER = 'uploads/'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER 
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0



@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')
    

@app.route('/upload')
def upload():
    return render_template('upload.html')

allow_files = [
    "filtered_feature_bc_matrix.h5",
    "tissue_hires_image.png",
    "tissue_lowres_image.png",
    "tissue_positions_list.csv",
    "scalefactors_json.json"
]

@app.route('/uploader', methods = ['GET', 'POST'])
def uploader_file():
    if request.method == 'POST':
        # Clean uploads folder before upload a new data
        import shutil
        shutil.rmtree(app.config['UPLOAD_FOLDER'])
        os.makedirs(app.config['UPLOAD_FOLDER'])

        # Get list of files from selected folder
        files = request.files.getlist("file")
        os.mkdir(os.path.join(app.config['UPLOAD_FOLDER'], "spatial"))
        for file in files:
            filename = secure_filename(file.filename)
            if allow_files[0] in filename:
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                os.rename(os.path.join(app.config['UPLOAD_FOLDER'], filename), 
                    os.path.join(app.config['UPLOAD_FOLDER'], allow_files[0]))
            for allow_file in allow_files[1:]:
                if allow_file in filename:
                    file.save(os.path.join(app.config['UPLOAD_FOLDER'] + "/spatial", filename))
                    os.rename(os.path.join(app.config['UPLOAD_FOLDER'] + "/spatial", filename), 
                        os.path.join(app.config['UPLOAD_FOLDER'] + "/spatial", allow_file))

    flash('File uploaded successfully')
    global adata
    adata = stlearn.Read10X(app.config['UPLOAD_FOLDER'])

    return redirect(url_for('index')) 

@app.route('/gene_plot')
def gene_plot():
    script = server_document('http://localhost:5006/bokeh_gene_plot')
    return render_template('gene_plot.html', script=script, template="Flask", 
        relative_urls=False)



# adata = stlearn.Read10X("../UQ/10X/BCBA")


def modify_doc(doc):
    from stlearn.plotting.classes_bokeh import BokehGenePlot
    gp_object = BokehGenePlot(adata)
    doc.add_root(row(gp_object.layout, width=800))
                   
    gp_object.data_alpha.on_change("value", gp_object.update_data)
    gp_object.tissue_alpha.on_change("value", gp_object.update_data)
    gp_object.spot_size.on_change("value", gp_object.update_data)  
    gp_object.gene_select.on_change("value", gp_object.update_data)


bkapp = Application(FunctionHandler(modify_doc))
bkapp2 = Application(FunctionHandler(modify_doc))


def bk_worker():
    asyncio.set_event_loop(asyncio.new_event_loop())

    server = Server({'/bokeh_gene_plot': bkapp,
                    '/bokeh_cluster_plot': bkapp2}, io_loop=IOLoop(), allow_websocket_origin=["localhost:5000"])
    server.start()
    server.io_loop.start()



from threading import Thread
Thread(target=bk_worker).start()

if __name__ == '__main__':
   app.run(debug = True)