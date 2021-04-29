from __future__ import division
import numpy as np
from PIL import Image
from bokeh.plotting import (
    figure,
    show,
    output_file,
    output_notebook,
    ColumnDataSource,
    curdoc,
)
from bokeh.models import (
    BoxSelectTool,
    LassoSelectTool,
    CustomJS,
    CheckboxGroup,
    LinearColorMapper,
    Slider,
    Panel,
    Select,
    AutocompleteInput,
    IndexFilter,
    CDSView,
    Div,
    Paragraph,
    Span,
)
from bokeh.palettes import Spectral11
from bokeh.layouts import column, row
import matplotlib
from matplotlib import pyplot as plt
import scanpy
import stlearn as st
import io
import os


def get_img_from_fig(fig, dpi=300):
    buf = io.BytesIO()
    from io import BytesIO

    fig.savefig(
        buf, format="jpg", dpi=dpi, bbox_inches="tight", pad_inches=0, transparent=False
    )
    buf.seek(0)
    img_arr = np.frombuffer(buf.getvalue(), dtype=np.uint8)
    img = np.asarray(Image.open(BytesIO(img_arr)))
    buf.close()
    # img = cv2.imdecode(img_arr, 1)
    # img = cv2.cvtColor(img, cv2.COLOR_BGRA2RGB)

    return img


def make_paga_plot(data, use_label):

    # img_pillow =Image.fromarray(get_img_from_fig(fig)).convert("RGBA")

    # xdim, ydim = img_pillow.size

    # print("Dimensions: ({xdim}, {ydim})".format(**locals()))
    # Create an array representation for the image `img`, and an 8-bit "4
    # layer/RGBA" version of it `view`.
    # img = np.empty((ydim, xdim), dtype=np.uint32)
    # view = img.view(dtype=np.uint8).reshape((ydim, xdim, 4))
    # Copy the RGBA image into view, flipping it so it comes right-side up
    # with a lower-left origin
    # view[:,:,:] = np.flipud(np.asarray(img_pillow))

    # Display the 32-bit RGBA image
    # dim = max(xdim, ydim)

    header = Div(
        text="""<h2>PAGA plot: </h2>""", width=400, height=30, sizing_mode="fixed"
    )

    fig, a = plt.subplots()

    try:
        del data.obsm["X_diffmap"]
    except:
        pass
    try:
        del data.obsm["X_draw_graph_fr"]
    except:
        pass

    data.uns[use_label + "_colors"] = data.uns["tmp_color"]
    data.uns["iroot"] = np.flatnonzero(data.obs[use_label] == str(0))[0]
    st.spatial.trajectory.pseudotime(data, eps=100, use_rep="X_pca")

    scanpy.pl.paga(
        data,
        color=use_label,
        ax=a,
        show=False,
        edge_width_scale=3.5,
        node_size_scale=3.5,
    )
    a.axis("off")
    fig.set_size_inches(5, 5)

    filelist = [
        f for f in os.listdir("stlearn_interactive/static/") if f.endswith(".png")
    ]
    for f in filelist:
        os.remove(os.path.join("stlearn_interactive/static/", f))

    fig.savefig(
        "stlearn_interactive/static/tmp.png", dpi=150, bbox_inches="tight", pad_inches=0
    )

    paga = Div(
        text="""<img src='stlearn_interactive/static/tmp.png' style='float: left;'>"""
    )

    fig2, a2 = plt.subplots()

    scanpy.tl.draw_graph(data, init_pos="paga")
    data.uns["louvain_colors"] = data.uns["tmp_color"]

    scanpy.pl.draw_graph(
        data,
        color="louvain",
        legend_loc="on data",
        ax=a2,
        show=False,
    )
    a2.axis("off")
    fig2.set_size_inches(5, 5)

    fig2.savefig(
        "stlearn_interactive/static/tmp2.png",
        dpi=150,
        bbox_inches="tight",
        pad_inches=0,
    )

    fa2 = Div(text="""<img src='stlearn_interactive/static/tmp2.png'>""")

    # if len(inputs.children) == 2:
    #    inputs.children.pop()

    figures = row(paga, fa2)

    inputs = column(header, figures)

    return inputs
