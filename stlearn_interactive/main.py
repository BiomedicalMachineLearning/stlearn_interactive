from bokeh.plotting import show, output_file, output_notebook, curdoc
from bokeh.models import (
    TextInput,
    Button,
    Paragraph,
    Select,
    Div,
    CheckboxGroup,
    Tabs,
    Panel,
    AutocompleteInput,
    RadioGroup,
    Slider,
)
from bokeh.layouts import column, row
import stlearn as st
import scanpy as sc

from script.gene_plot import make_gene_plot
from script.cluster_plot import make_cluster_plot
from script.paga_plot import make_paga_plot
from script.pseudotime_plot import make_ptvis_plot
from os import walk

import matplotlib
from matplotlib import pyplot as plt

##################### Tab 1 #############################
header_step1 = Div(text="""<h2>Step 1: Reading data </h2>""", width=400, height=50)

header_step2 = Div(text="""<h2>Step 2: Preprocessing </h2>""", width=400, height=50)

# Setup input
file_input = TextInput(
    value="/home/d.pham/10X/BCBA", title="10X Visium folder path (default):"
)


h5_input = TextInput(
    value="filtered_feature_bc_matrix.h5",
    title="Count matrix .h5 file path (optional):",
)

# Setup button
submit_bt = Button(label="Read data", button_type="success")
preprocessing_bt = Button(label="Preprocessing", button_type="success")

# Setup error
read_error = Paragraph(text="""""", width=0, height=0)

# Setup checkbox preprocessing
preprocessing = CheckboxGroup(
    labels=["Filter genes", "Normalize total", "Log 1P", "Scale"], active=[0, 1, 2, 3]
)

##################### Tab 3 #############################

# Setup setting for clustering
header_step3 = Div(
    text="""<h2>Step 3: SME Clustering </h2>""",
    width=400,
    height=50,
    sizing_mode="fixed",
)


pca_slider = Slider(
    title="PCA - Number of component",
    value=50,
    start=1,
    end=100,
    step=1,
    sizing_mode="fixed",
)

sme_clust = RadioGroup(
    labels=["Using SMEClust", "Don't use"], active=0, sizing_mode="fixed"
)

tiling_slider = Slider(
    title="Crop size for tiling image",
    value=40,
    start=1,
    end=100,
    step=1,
    sizing_mode="fixed",
)

cl_method = Select(
    title="Choosing clustering method:",
    value="Louvain clustering",
    options=["Louvain clustering", "K-means clustering"],
    sizing_mode="fixed",
)

# Setup for louvain clustering

knn_slider = Slider(
    title="Number of neighbors for KNN",
    value=25,
    start=2,
    end=100,
    step=1,
    sizing_mode="fixed",
)

morpho_radius = Slider(
    title="Radius for adjusting spots",
    value=50,
    start=2,
    end=100,
    step=1,
    sizing_mode="fixed",
)

# Setup for k-means clustering

n_clusters = Slider(
    title="Number of clusters", value=10, start=2, end=100, step=1, sizing_mode="fixed"
)

sme_bt = Button(label="Do clustering", button_type="success", sizing_mode="fixed")

log_clustering = Paragraph(text="""""", width=0, height=0, sizing_mode="fixed")


##################### Tab 4 #############################

# Setup setting for pseudotime analysis
header_step4 = Div(
    text="""<h2>Step 4.1: Pseudotime analysis </h2>""",
    width=400,
    height=50,
    sizing_mode="fixed",
)

##################### Tab 5 #############################

# Setup setting for pseudotime visualization
header_step5 = Div(
    text="""<h2>Step 4.2: Pseudotime analysis ST visualization </h2>""",
    width=400,
    height=50,
    sizing_mode="fixed",
)


def reading_data():
    try:
        global data
        data = st.Read10X(path=file_input.value)

        layout.children.append(header_step2)
        layout.children.append(preprocessing)
        layout.children.append(preprocessing_bt)

        read_error.text = """Done!"""

    except:
        read_error.text = """Error! Your input path is not correct. Please try again."""
        read_error.style = {"color": "red"}
        read_error.width = 400
        read_error.height = 25


def change_click():

    read_error.width = 400
    read_error.height = 25
    read_error.text = """Reading data..."""
    read_error.style = {"color": "black"}

    # layout.update(tab1_layout)
    # tabs.update(tabs=tab_list)
    curdoc().add_next_tick_callback(reading_data)


def change_preprocessing():

    if 0 in preprocessing.active:
        sc.pp.filter_genes(data, min_cells=3)
        print("Filter genes!")
    if 1 in preprocessing.active:
        sc.pp.normalize_total(data)
        print("Normalized genes!")
    if 2 in preprocessing.active:
        sc.pp.log1p(data)
        print("Log transformed genes!")
    if 3 in preprocessing.active:
        sc.pp.scale(data)
        print("Scaled genes!")

    if len(tab_list) == 1:
        tab2 = make_gene_plot(data)
        tab_list.append(tab2)
        tab_list.append(tab3)
        tabs.update(tabs=tab_list)
        tabs.active += 1
    else:
        tabs.active += 1


####### Step 3 #######


def doing_smeclust():

    # Do PCA
    st.em.run_pca(data, n_comps=pca_slider.value, random_state=0)

    # Using SMEClust
    if sme_clust.active == 0:

        st.pp.tiling(data, out_path="../tiling", crop_size=tiling_slider.value)
        st.pp.extract_feature(data)
        st.spatial.morphology.adjust(
            data, use_data="X_pca", radius=morpho_radius.value, method="mean"
        )

        if cl_method.value == "Louvain clustering":
            st.pp.neighbors(
                data,
                n_neighbors=knn_slider.value,
                use_rep="X_pca_morphology",
                random_state=0,
            )
            st.tl.clustering.louvain(data, random_state=0)
            st.pl.cluster_plot(data, use_label="louvain", show_plot=False)
            get_cluster_color(data, use_label="louvain")

            plot_cluster = make_cluster_plot(data, "louvain")

        else:
            st.tl.clustering.kmeans(
                data, n_clusters=n_clusters.value, use_data="X_pca_morphology"
            )
            st.pl.cluster_plot(data, use_label="kmeans", show_plot=False)
            get_cluster_color(data, use_label="kmeans")
            plot_cluster = make_cluster_plot(data, "kmeans")
    else:
        print("Choose non-SME")
        if cl_method.value == "Louvain clustering":
            st.pp.neighbors(
                data, n_neighbors=knn_slider.value, use_rep="X_pca", random_state=0
            )
            st.tl.clustering.louvain(data, random_state=0)
            st.pl.cluster_plot(data, use_label="louvain", show_plot=False)
            get_cluster_color(data, use_label="louvain")
            plot_cluster = make_cluster_plot(data, "louvain")

        else:
            st.tl.clustering.kmeans(data, n_clusters=n_clusters.value, use_data="X_pca")
            st.pl.cluster_plot(data, use_label="kmeans", show_plot=False)
            get_cluster_color(data, use_label="kmeans")
            plot_cluster = make_cluster_plot(data, "kmeans")

    log_clustering.text = """Done!"""

    if len(sme_layout.children) == 8:
        sme_layout.children.pop()

    sme_layout.children.append(plot_cluster)

    if len(tab_list) == 3:

        plot_paga = make_paga_plot(data, "louvain")

        if len(pt_layout.children) == 2:
            pt_layout.children.pop()

        pt_layout.children.append(plot_paga)

        tab_list.append(tab4)
        tabs.update(tabs=tab_list)

    if len(tab_list) == 4:
        tab5 = make_ptvis_plot(data)
        tab_list.append(tab5)
        tabs.update(tabs=tab_list)


def change_smeclust():

    log_clustering.width = 400
    log_clustering.height = 50
    log_clustering.text = """Waiting for few minutes ..."""
    log_clustering.style = {"color": "black"}

    # layout.update(tab1_layout)
    # tabs.update(tabs=tab_list)
    curdoc().add_next_tick_callback(doing_smeclust)


def choose_sme(attrname, old, new):

    if sme_clust.active == 0:
        insert_at = sme_layout.children.index(sme_clust) + 1
        sme_layout.children[insert_at:insert_at] = [tiling_slider]
        insert_at += 1
        sme_layout.children[insert_at:insert_at] = [morpho_radius]
    else:
        if tiling_slider in sme_layout.children:
            del_index = sme_layout.children.index(tiling_slider)
            del sme_layout.children[del_index]
            del sme_layout.children[del_index]


def choose_clustering(attrname, old, new):
    if cl_method.value == "Louvain clustering":

        # Louvain
        insert_at = sme_layout.children.index(cl_method) + 1
        sme_layout.children[insert_at:insert_at] = [knn_slider]

        # Kmeans
        if n_clusters in sme_layout.children:
            del_index = sme_layout.children.index(n_clusters)
            del sme_layout.children[del_index]

    else:
        # Louvain
        if knn_slider in sme_layout.children:
            del_index = sme_layout.children.index(knn_slider)
            del sme_layout.children[del_index]

        # Kmeans
        insert_at = sme_layout.children.index(cl_method) + 1
        sme_layout.children[insert_at:insert_at] = [n_clusters]


submit_bt.on_click(change_click)

sme_clust.on_change("active", choose_sme)
cl_method.on_change("value", choose_clustering)
sme_bt.on_click(change_smeclust)

preprocessing_bt.on_click(change_preprocessing)

# pt_bt.on_click(change_click_pt)
#######################################

tab1_layout = [header_step1, file_input, h5_input, submit_bt, read_error]
layout = column(tab1_layout)
tab1 = Panel(child=layout, title="Reading data")

tab3_layout = [
    header_step3,
    pca_slider,
    sme_clust,
    tiling_slider,
    morpho_radius,
    cl_method,
    knn_slider,
    sme_bt,
    log_clustering,
]

sme_layout = column(tab3_layout)
tab3 = Panel(child=sme_layout, title="SME Clustering")

tab4_layout = [header_step4]
pt_layout = column(tab4_layout)
tab4 = Panel(child=pt_layout, title="Pseudotime preview")

tab_list = [tab1]

# Put all the tabs into one application
tabs = Tabs(tabs=tab_list)

curdoc().add_root(tabs)
curdoc().title = "stlearn_interactive"


def get_cluster_color(data, use_label="louvain", cmap="vega_20_scanpy"):

    from scanpy.plotting import palettes

    if cmap == "vega_10_scanpy":
        cmap = palettes.vega_10_scanpy
    elif cmap == "vega_20_scanpy":
        cmap = palettes.vega_20_scanpy
    elif cmap == "default_102":
        cmap = palettes.default_102
    elif cmap == "default_28":
        cmap = palettes.default_28
    else:
        raise ValueError(
            "We only support vega_10_scanpy, vega_20_scanpy, default_28, default_102"
        )

    cmaps = matplotlib.colors.LinearSegmentedColormap.from_list("", cmap)

    cmap_ = plt.cm.get_cmap(cmaps)

    # Plot scatter plot based on pixel of spots
    data.uns["tmp_color"] = []

    for i, cluster in enumerate(data.obs.groupby(use_label)):

        data.uns["tmp_color"].append(matplotlib.colors.to_hex(cmap_(int(i) / 19)))
