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
    HoverTool,
)
from bokeh.palettes import Spectral11
from bokeh.layouts import column, row
import matplotlib
from matplotlib import pyplot as plt
from collections import OrderedDict


def make_cluster_plot(data, use_label):

    library_id = list(data.uns["spatial"].keys())[0]
    # Open image, and make sure it's RGB*A*
    image = (
        data.uns["spatial"][library_id]["images"][data.uns["spatial"]["use_quality"]]
        * 255
    ).astype(np.uint8)

    img_pillow = Image.fromarray(image).convert("RGBA")

    xdim, ydim = img_pillow.size
    print("Dimensions: ({xdim}, {ydim})".format(**locals()))
    # Create an array representation for the image `img`, and an 8-bit "4
    # layer/RGBA" version of it `view`.
    img = np.empty((ydim, xdim), dtype=np.uint32)
    view = img.view(dtype=np.uint8).reshape((ydim, xdim, 4))
    # Copy the RGBA image into view, flipping it so it comes right-side up
    # with a lower-left origin
    view[:, :, :] = np.flipud(np.asarray(img_pillow))

    # Display the 32-bit RGBA image
    dim = max(xdim, ydim)

    list_index = {}
    tmp = list(data.obs[use_label])

    for cl in np.sort(data.obs[use_label].unique().astype(int)).astype(str):
        list_index[cl] = [i for i, x in enumerate(tmp) if x == cl]

    def make_fig(list_index):
        fig = figure(
            title="Cluster plot",
            x_range=(0, dim),
            y_range=(dim, 0),
            # Specifying xdim/ydim isn't quire right :-(
            # width=xdim, height=ydim,
        )

        fig.image_rgba(
            image=[img], x=0, y=xdim, dw=ydim, dh=xdim, global_alpha=tissue_alpha.value
        )

        # Get query clusters
        command = []
        for i in list_cluster.active:
            command.append(use_label + ' == "' + str(i) + '"')
        tmp = data.obs.query(" or ".join(command))

        coordinate = tmp[["imagecol", "imagerow"]].values
        x = coordinate[:, 0]
        y = coordinate[:, 1]

        category_items = np.sort(data.obs[use_label].unique().astype(int)).astype(str)
        palette = data.uns["tmp_color"]
        colormap = dict(zip(category_items, palette))
        color = list(tmp[use_label].map(colormap))
        cluster = list(tmp[use_label])

        s1 = ColumnDataSource(data=dict(x=x, y=y, color=color, cluster=cluster))

        fig.scatter(
            x="x",
            y="y",
            source=s1,
            size=5,
            color="color",
            legend_group="cluster",
            fill_alpha=data_alpha.value,
            line_alpha=data_alpha.value,
        )

        fig.toolbar.logo = None
        fig.xaxis.visible = False
        fig.yaxis.visible = False
        fig.xgrid.grid_line_color = None
        fig.ygrid.grid_line_color = None
        fig.outline_line_alpha = 0
        fig.add_tools(LassoSelectTool())
        fig.add_tools(BoxSelectTool())
        fig.add_tools(HoverTool())

        hover = fig.select(dict(type=HoverTool))
        hover.tooltips = OrderedDict(
            [
                ("Spot", "$index"),
                ("X location", "@x{1.11}"),
                ("Y location", "@y{1.11}"),
                ("Cluster", "@cluster"),
            ]
        )

        return fig

    def update_data(attrname, old, new):
        layout.children[1] = make_fig(list_index)

    header = Div(
        text="""<h2>Clustering result plot: </h2>""",
        width=400,
        height=30,
        sizing_mode="fixed",
    )

    data_alpha = Slider(title="Spot alpha", value=1.0, start=0, end=1.0, step=0.1)
    data_alpha.on_change("value", update_data)

    tissue_alpha = Slider(title="Tissue alpha", value=1.0, start=0, end=1.0, step=0.1)
    tissue_alpha.on_change("value", update_data)

    p = Paragraph(text="""Choose clusters:""", width=400, height=20)
    list_cluster = CheckboxGroup(
        labels=list(np.sort(data.obs["louvain"].unique().astype(int)).astype(str)),
        active=list(np.sort(data.obs["louvain"].unique().astype(int))),
    )
    list_cluster.on_change("active", update_data)

    inputs = column(header, data_alpha, tissue_alpha, p, list_cluster)

    layout = row(inputs, make_fig(list_index))

    return layout


def centroidpython(data):
    x, y = zip(*data)
    l = len(x)
    return sum(x) / l, sum(y) / l
