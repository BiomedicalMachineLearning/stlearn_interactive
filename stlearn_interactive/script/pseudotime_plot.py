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
    Button,
    RadioGroup,
)
from bokeh.palettes import Spectral11
from bokeh.layouts import column, row
import matplotlib
from matplotlib import pyplot as plt
import stlearn as st
import os
import random


def make_ptvis_plot(data):

    list_index = {}
    tmp = list(data.obs["louvain"])

    for cl in np.sort(data.obs["louvain"].unique().astype(int)).astype(str):
        list_index[cl] = [i for i, x in enumerate(tmp) if x == cl]

    def doing_ps(root, spot_index, eps, list_cluster, data_alpha, tissue_alpha):

        fig, a = plt.subplots()

        try:
            del data.obsm["X_diffmap"]
        except:
            pass
        try:
            del data.obsm["X_draw_graph_fr"]
        except:
            pass
        st.pl.cluster_plot(data, use_label="louvain", show_plot=False)
        data.uns["iroot"] = np.flatnonzero(data.obs["louvain"] == str(root))[spot_index]
        st.spatial.trajectory.pseudotime(data, eps=eps, use_rep="X_pca")

        if len(inputs.children) == 10:
            inputs.children.pop()
            final_layout.children.pop()
            final_layout.children.pop()

        id_img = random.randint(3, 99999999)

        # filelist = [ f for f in os.listdir("stlearn_interactive/static/") if f.endswith(".png") ]
        # for f in filelist:
        # os.remove(os.path.join("stlearn_interactive/static/", f))

        st.pl.trajectory.pseudotime_plot(
            data,
            list_cluster=list_cluster,
            show_graph=False,
            node_alpha=1,
            data_alpha=data_alpha,
            tissue_alpha=tissue_alpha,
            edge_alpha=0.1,
            node_size=3,
            show_plot=False,
            dpi=150,
            output="stlearn_interactive/static/",
            name="tmp" + str(id_img),
        )
        fig = Div(
            text="<img src='stlearn_interactive/static/tmp"
            + str(id_img)
            + ".png' height='482' width='568'>",
            sizing_mode="fixed",
        )

        inputs.children.insert(10, fig)

        final_layout.children.append(pst_layout)
        final_layout.children.append(pstg_layout)

        log_clustering.text = """Done!"""

    def change_ps():
        log_clustering.width = 400
        log_clustering.height = 50
        log_clustering.text = """Processing data ..."""
        log_clustering.style = {"color": "black"}

        curdoc().add_next_tick_callback(
            doing_ps(
                root.value,
                spot_index.value,
                eps.value,
                list_cluster.active,
                data_alpha.value,
                tissue_alpha.value,
            )
        )

    def update_spot(attrname, old, new):

        len_index = len(data.obs[data.obs["louvain"] == root.value])
        spot_index = Slider(
            start=0,
            end=len_index,
            value=0,
            step=1,
            title="Spot index",
            sizing_mode="fixed",
        )

        # insert_at = inputs.children.index(root) + 1
        inputs.children.insert(2, spot_index)

        if spot_index == inputs.children[2]:
            del inputs.children[3]

    def create_figures(data, cluster, reverse, tissue_alpha, data_alpha):
        fig, a = plt.subplots()

        st.spatial.trajectory.local_level(data, use_label="louvain", cluster=cluster)

        # if len(pst_layout.children) == 8:
        #    del pst_layout.children[7]

        id_img = random.randint(3, 99999999)
        id_img2 = random.randint(3, 99999999)

        # filelist = [ f for f in os.listdir("stlearn_interactive/static/") if f.endswith(".png") ]
        # for f in filelist:
        # os.remove(os.path.join("stlearn_interactive/static/", f))

        if reverse == 0:
            reverse_value = True
        else:
            reverse_value = False

        st.pl.subcluster_plot(
            data,
            use_label="louvain",
            cluster=cluster,
            tissue_alpha=tissue_alpha,
            data_alpha=data_alpha,
            dpi=150,
            output="stlearn_interactive/static/",
            name="tmp" + str(id_img),
            show_plot=False,
        )

        st.pl.trajectory.local_plot(
            data,
            use_cluster=cluster,
            branch_alpha=0.2,
            dpi=150,
            reverse=reverse_value,
            output="stlearn_interactive/static/",
            name="tmp" + str(id_img2),
            show_plot=False,
        )

        fig = Div(
            text="<img src='stlearn_interactive/static/tmp"
            + str(id_img)
            + ".png' height='554' width='554' style='float: left;'>",
            sizing_mode="fixed",
        )

        fig2 = Div(
            text="<img src='stlearn_interactive/static/tmp"
            + str(id_img2)
            + ".png' height='577' width='600' style='padding-left: 300px;'>",
            sizing_mode="fixed",
        )

        return [fig, fig2]

    def change_pst(attrname, old, new):

        if len(pst_layout.children) == 7:
            pst_layout.children.pop()

        fig, fig2 = create_figures(
            data,
            cluster_pst.value,
            reverse.active,
            tissue_alpha_pst.value,
            data_alpha_pst.value,
        )

        figures = row([fig, fig2])
        pst_layout.children.insert(6, figures)

    def doing_pstg(list_cluster, data_alpha, tissue_alpha):

        fig, a = plt.subplots()

        if len(pstg_layout.children) == 8:
            pstg_layout.children.pop()

        id_img = random.randint(3, 99999999)
        id_img2 = random.randint(3, 99999999)

        filelist = [
            f for f in os.listdir("stlearn_interactive/static/") if f.endswith(".png")
        ]
        for f in filelist:
            os.remove(os.path.join("stlearn_interactive/static/", f))

        st.spatial.trajectory.pseudotimespace_global(
            data, use_label="louvain", list_cluster=list_cluster
        )

        st.pl.cluster_plot(
            data,
            use_label="louvain",
            show_trajectory=True,
            list_cluster=list_cluster,
            show_plot=False,
            show_subcluster=True,
            show_legend=False,
            dpi=150,
            output="stlearn_interactive/static/",
            name="tmp" + str(id_img),
            data_alpha=data_alpha,
            tissue_alpha=tissue_alpha,
        )

        st.pl.trajectory.tree_plot(
            data,
            dpi=150,
            output="stlearn_interactive/static/",
            name="tmp" + str(id_img2),
            show_plot=False,
        )

        fig = Div(
            text="<img src='stlearn_interactive/static/tmp"
            + str(id_img)
            + ".png' height='500' width='500' style='float: left;'>",
            sizing_mode="fixed",
        )

        fig2 = Div(
            text="<img src='stlearn_interactive/static/tmp"
            + str(id_img2)
            + ".png' height='472' width='1200' style='padding-left: 300px;'>",
            sizing_mode="fixed",
        )

        figures = row([fig, fig2])

        pstg_layout.children.insert(8, figures)

        log_clustering.text = """Done!"""

    def change_pstg():
        log_clustering.width = 400
        log_clustering.height = 50
        log_clustering.text = """Processing data ..."""
        log_clustering.style = {"color": "black"}

        curdoc().add_next_tick_callback(
            doing_pstg(
                list_cluster_pstg.active,
                data_alpha_pstg.value,
                tissue_alpha_pstg.value,
            )
        )

    log_clustering = Paragraph(text="""""", width=0, height=0, sizing_mode="fixed")

    header = Div(
        text="""<h2>Step 4.2: Pseudotime analysis </h2>""",
        width=400,
        height=30,
        sizing_mode="fixed",
    )

    data_alpha = Slider(
        title="Spot alpha", value=1.0, start=0, end=1.0, step=0.1, sizing_mode="fixed"
    )
    # data_alpha.on_change('value',update_data)

    tissue_alpha = Slider(
        title="Tissue alpha", value=1.0, start=0, end=1.0, step=0.1, sizing_mode="fixed"
    )
    # tissue_alpha.on_change('value',update_data)

    root = Select(
        title="Choose root:",
        value="0",
        options=list(np.sort(data.obs["louvain"].unique().astype(int)).astype(str)),
    )
    root.on_change("value", update_spot)

    spot_index = Slider(
        start=0,
        end=len(data.obs[data.obs["louvain"] == "0"]),
        value=0,
        step=1,
        title="Spot index",
        sizing_mode="fixed",
    )

    eps = Slider(
        title="Eps (DBSCAN)", value=100, start=0, end=500, step=1, sizing_mode="fixed"
    )

    p = Paragraph(
        text="""Choose clusters to display:""",
        width=400,
        height=20,
        sizing_mode="fixed",
    )

    list_cluster = CheckboxGroup(
        labels=list(np.sort(data.obs["louvain"].unique().astype(int)).astype(str)),
        active=list(np.sort(data.obs["louvain"].unique().astype(int))),
    )
    # list_cluster.on_change('active',update_data)

    ps_bt = Button(
        label="Do pseudotime analysis", button_type="success", sizing_mode="fixed"
    )

    ps_bt.on_click(change_ps)

    inputs = column(
        header, root, spot_index, eps, p, list_cluster, data_alpha, tissue_alpha, ps_bt
    )

    ###############################################

    header_pst = Div(
        text="""<h2>Pseudo-space-time. Local analysis </h2>""",
        width=400,
        height=30,
        sizing_mode="fixed",
    )

    cluster_pst = Select(
        title="Choose clusters:",
        value="0",
        options=list(np.sort(data.obs["louvain"].unique().astype(int)).astype(str)),
        sizing_mode="fixed",
    )
    cluster_pst.on_change("value", change_pst)

    reverse = RadioGroup(labels=["Reverse", "Inverse"], active=0, sizing_mode="fixed")
    reverse.on_change("active", change_pst)

    data_alpha_pst = Slider(
        title="Spot alpha", value=1.0, start=0, end=1.0, step=0.1, sizing_mode="fixed"
    )
    data_alpha_pst.on_change("value", change_pst)

    tissue_alpha_pst = Slider(
        title="Tissue alpha", value=1.0, start=0, end=1.0, step=0.1, sizing_mode="fixed"
    )
    tissue_alpha_pst.on_change("value", change_pst)

    pst_layout = column(
        [header_pst, cluster_pst, reverse, data_alpha_pst, tissue_alpha_pst]
    )

    figures = row(
        create_figures(
            data,
            cluster_pst.value,
            reverse.active,
            tissue_alpha_pst.value,
            data_alpha_pst.value,
        )
    )

    pst_layout.children.append(figures)

    ###############################################

    header_pstg = Div(
        text="""<h2>Pseudo-space-time. Global analysis </h2>""",
        width=400,
        height=30,
        sizing_mode="fixed",
    )

    p = Paragraph(text="""Choose clusters:""", width=400, height=20)

    list_cluster_pstg = CheckboxGroup(
        labels=list(np.sort(data.obs["louvain"].unique().astype(int)).astype(str)),
        sizing_mode="fixed",
    )
    # list_cluster_pstg.on_change('value',change_pstg)

    data_alpha_pstg = Slider(
        title="Spot alpha", value=1.0, start=0, end=1.0, step=0.1, sizing_mode="fixed"
    )
    # data_alpha_pstg.on_change('value',change_pstg)

    tissue_alpha_pstg = Slider(
        title="Tissue alpha", value=1.0, start=0, end=1.0, step=0.1, sizing_mode="fixed"
    )
    # tissue_alpha_pstg.on_change('value',change_pstg)

    ps_bt = Button(
        label="Do global pseudo-space-time analysis",
        button_type="success",
        sizing_mode="fixed",
    )
    ps_bt.on_click(change_pstg)

    pstg_layout = column(
        [
            header_pstg,
            p,
            list_cluster_pstg,
            data_alpha_pstg,
            tissue_alpha_pstg,
            ps_bt,
        ]
    )

    ###############################################

    final_layout = column([inputs])

    tab = Panel(child=final_layout, title="Pseudotime analysis")

    return tab
