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
    Div,
    LinearColorMapper,
    Slider,
    Panel,
    Select,
    AutocompleteInput,
    ColorBar,
    TextInput,
    BasicTicker,
    CrosshairTool,
    HoverTool,
)
from bokeh.palettes import Spectral11
from bokeh.layouts import column, row
from collections import OrderedDict


def make_gene_plot(data):

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

    gene_list = list(data.var_names)

    def make_fig():

        fig = figure(
            title="Gene plot",
            x_range=(0, dim - 150),
            y_range=(dim, 0),
        )

        colors = _gene_plot(data, "CumSum", [gene_select.value])

        if threshold.value != "":
            colors = colors[colors > float(threshold.value)]

        index_filter = colors.index

        filter_obs = data.obs.loc[index_filter]
        imagecol = filter_obs["imagecol"].values
        imagerow = filter_obs["imagerow"].values

        # coordinate = data.obs[["imagecol","imagerow"]].values
        # x = coordinate[:,0]
        # y = coordinate[:,1]

        # color = _gene_plot(data,"CumSum",[gene_list[0]])
        s1 = ColumnDataSource(data=dict(x=imagecol, y=imagerow, color=colors))

        fig.image_rgba(
            image=[img], x=0, y=xdim, dw=ydim, dh=xdim, global_alpha=tissue_alpha.value
        )

        color_mapper = LinearColorMapper(
            palette=Spectral11, low=min(s1.data["color"]), high=max(s1.data["color"])
        )

        fig.circle(
            "x",
            "y",
            color={"field": "color", "transform": color_mapper},
            size=5,
            source=s1,
            fill_alpha=data_alpha.value,
            line_alpha=data_alpha.value,
        )

        color_bar = ColorBar(
            color_mapper=color_mapper, ticker=BasicTicker(), location=(10, 0)
        )
        fig.add_layout(color_bar, "right")

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
                ("Gene expression", "@color{1.11}"),
            ]
        )

        return fig

    def update_data(attrname, old, new):
        layout.children[1] = make_fig()

    header = Div(
        text="""<h2>Gene plot: </h2>""", width=400, height=50, sizing_mode="fixed"
    )

    data_alpha = Slider(title="Spot alpha", value=1.0, start=0, end=1.0, step=0.1)
    data_alpha.on_change("value", update_data)

    tissue_alpha = Slider(title="Tissue alpha", value=1.0, start=0, end=1.0, step=0.1)
    tissue_alpha.on_change("value", update_data)

    threshold = TextInput(title="Threshold:", value="")
    threshold.on_change("value", update_data)

    gene_select = AutocompleteInput(
        title="Gene:", value=gene_list[0], completions=gene_list, min_characters=1
    )
    gene_select.on_change("value", update_data)

    inputs = column(header, gene_select, data_alpha, tissue_alpha, threshold)

    layout = row(inputs, make_fig())

    # Make a tab with the layout
    tab = Panel(child=layout, title="Gene plot")

    return tab


def _gene_plot(adata, method, genes):

    # Gene plot option

    if len(genes) == 0:
        raise ValueError("Genes shoule be provided, please input genes")

    elif len(genes) == 1:

        if genes[0] not in adata.var.index:
            raise ValueError(
                genes[0] + " is not exist in the data, please try another gene"
            )

        colors = adata[:, genes].to_df().iloc[:, -1]

        return colors
    else:

        for gene in genes:
            if gene not in adata.var.index:
                genes.remove(gene)
                warnings.warn(
                    "We removed " + gene + " because they not exist in the data"
                )

            if len(genes) == 0:
                raise ValueError("All provided genes are not exist in the data")

        count_gene = adata[:, genes].to_df()

        if method is None:
            raise ValueError(
                "Please provide method to combine genes by NaiveMean/NaiveSum/CumSum"
            )

        if method == "NaiveMean":
            present_genes = (count_gene > 0).sum(axis=1) / len(genes)

            count_gene = (count_gene.mean(axis=1)) * present_genes
        elif method == "NaiveSum":
            present_genes = (count_gene > 0).sum(axis=1) / len(genes)

            count_gene = (count_gene.sum(axis=1)) * present_genes

        elif method == "CumSum":
            count_gene = count_gene.cumsum(axis=1).iloc[:, -1]

        colors = count_gene

        return colors
