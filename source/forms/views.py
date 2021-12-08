""" This is more a general views focussed on defining functions which are \
	called by other views for specify pages. This way different pages can be \
	used to display different data, but in a consistent way.
"""

import sys
import numpy
import numpy as np
from flask import flash
from source.forms import forms

from source.forms.utils import flash_errors, savePlot
import source.forms.view_helpers as vhs
import traceback

from flask import render_template

import scanpy as sc
import stlearn as st

from scipy.spatial.distance import cosine

# Creating the forms using a class generator #
PreprocessForm = forms.getPreprocessForm()
#CCIForm = forms.getCCIForm() #OLD
ClusterForm = forms.getClusterForm()
LRForm = forms.getLRForm()

def run_preprocessing(request, adata, step_log):
    """Performs the scanpy pre-processing steps based on the inputted data."""

    form = PreprocessForm(request.form)

    if not form.validate_on_submit():
        flash_errors(form)

    elif type(adata) == type(None):
        flash("Need to load data first!")

    else:
        # Logging params used #
        step_log["preprocessed_params"] = vhs.getData(form)
        print(step_log["preprocessed_params"], file=sys.stdout)

        # QC filtering #
        sc.pp.filter_cells(adata, min_genes=vhs.getVal(form, "Minimum genes per spot"))
        sc.pp.filter_cells(
            adata, min_counts=vhs.getVal(form, "Minimum counts per spot")
        )
        sc.pp.filter_genes(adata, min_cells=vhs.getVal(form, "Minimum spots per gene"))
        sc.pp.filter_genes(
            adata, min_counts=vhs.getVal(form, "Minimum counts per gene")
        )

        # Pre-processing #
        if vhs.getVal(form, "Normalize total"):
            sc.pp.normalize_total(adata, target_sum=1e4)
        if vhs.getVal(form, "Log 1P"):
            sc.pp.log1p(adata)
        adata.raw = adata
        if vhs.getVal(form, "Scale"):
            sc.pp.scale(adata, max_value=10)

        # Setting pre-process to true #
        step_log["preprocessed"][0] = True

    if step_log["preprocessed"][0]:
        flash("Preprocessing is completed!")

    updated_page = render_template(
        "preprocessing.html",
        title=step_log["preprocessed"][1],
        preprocess_form=form,
        flash_bool=step_log["preprocessed"][0],
        step_log=step_log,
    )

    return updated_page

def run_lr(request, adata, step_log):
    """ Runs LR analysis.
    """

    form = LRForm(request.form)

    if not form.validate_on_submit():
        flash_errors(form)

    elif type(adata) == type(None):
        flash("Need to load data first!")

    else:
        step_log["lr_params"] = vhs.getData(form)
        print(step_log["lr_params"], file=sys.stdout)
        elements = numpy.array(list(step_log["lr_params"].keys()))
        # order: Species, Spot neighbourhood, min_spots, n_pairs, CPUs
        element_values = list(step_log["lr_params"].values())
        dist = element_values[1]
        dist = dist if dist!=-1 else None

        # Loading the LR databases available within stlearn (from NATMI)
        lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'],
                                 species=element_values[0])

        # Running the analysis #
        st.tl.cci.run(adata, lrs, min_spots=element_values[2],
                      distance=dist, n_pairs=element_values[3],
                      n_cpus=element_values[-1])

        flash("LR analysis is completed!")

    step_log["lr"][0] = 'lr_summary' in adata.uns

    updated_page = render_template(
        "lr.html",
        title=step_log["lr"][1],
        lr_form=form,
        flash_bool=True,
        step_log=step_log,
    )
    return updated_page

def run_cci(request, adata, step_log):
    """Performs CCI analysis."""

    CCIForm = forms.getCCIForm(adata)
    form = CCIForm(request.form)

    if not form.validate_on_submit():
        flash_errors(form)

    elif type(adata) == type(None):
        flash("Need to load data first!")

    else:
        step_log["cci_params"] = vhs.getData(form)
        print(step_log["cci_params"], file=sys.stdout)
        elements = numpy.array(list(step_log["cci_params"].keys()))
        # order: cell_type, min_spots, spot_mixtures, cell_prop_cutoff, sig_spots
        #           n_perms
        element_values = list(step_log["cci_params"].values())

        if not form.validate_on_submit():
            flash_errors(form)

        else:
            try:
                # Running the counting of co-occurence of cell types and LR expression #
                st.tl.cci.run_cci(adata, element_values[0],
                                  min_spots=element_values[1],
                                  spot_mixtures=element_values[2],
                                  cell_prop_cutoff=element_values[3],
                                  sig_spots=True,#Should make this not optional..
                                  n_perms=element_values[4]
                                  )

                flash("CCI analysis is completed!")

            except Exception as msg:
                traceback.print_exc(file=sys.stdout)
                flash("Analysis ERROR: " + str(msg))
                print(msg)

    step_log["cci"][0] = np.any(['lr_cci_' in key for key in adata.uns])

    updated_page = render_template(
        "cci.html",
        title=step_log["cci"][1],
        cci_form=form,
        flash_bool=True,
        step_log=step_log,
    )

    return updated_page


def run_clustering(request, adata, step_log):
    """Performs clustering analysis."""

    form = ClusterForm(request.form)

    step_log["cluster_params"] = vhs.getData(form)
    print(step_log["cluster_params"], file=sys.stdout)
    elements = list(step_log["cluster_params"].keys())
    # order: pca_comps, SME bool, method, method_param
    element_values = list(step_log["cluster_params"].values())

    if not form.validate_on_submit():
        flash_errors(form)

    elif type(adata) == type(None):
        flash("Need to load data first!")

    else:
        try:
            # Running PCA, performs scaling internally #
            n_comps = element_values[0]
            st.em.run_pca(adata, n_comps=n_comps)

            print(element_values[1], file=sys.stdout, flush=True)
            if element_values[1]:  # Performing SME clustering #
                # Image feature extraction #
                st.pp.tiling(adata)
                st.pp.extract_feature(adata)

                # apply stSME to data (format of data depending on preprocess)
                st.spatial.SME.SME_normalize(adata, use_data="raw")
                adata.X = adata.obsm["raw_SME_normalized"]
                st.em.run_pca(adata, n_comps=n_comps)

            # Performing the clustering on the PCA #
            if element_values[2] == "KMeans":  # KMeans
                param = int(element_values[3])
                st.tl.clustering.kmeans(adata, n_clusters=param, use_data="X_pca")

                st.pp.neighbors(adata, n_neighbors=element_values[5], use_rep="X_pca")
                sc.tl.paga(adata, groups="kmeans")
                st.pl.cluster_plot(adata, use_label="kmeans")

            elif element_values[2] == "Louvain":  # Louvain
                param = element_values[4]
                st.pp.neighbors(adata, n_neighbors=element_values[5], use_rep="X_pca")
                st.tl.clustering.louvain(adata, resolution=param)
                sc.tl.paga(adata, groups="louvain")
                st.pl.cluster_plot(adata, use_label="louvain")

            else:  # Leiden
                param = element_values[4]
                st.pp.neighbors(adata, n_neighbors=element_values[5], use_rep="X_pca")
                sc.tl.leiden(adata, resolution=param)
                sc.tl.paga(adata, groups="leiden")
                st.pl.cluster_plot(adata, use_label="leiden")

            savePlot("clustering.png")

            step_log["clustering"][0] = True
            flash("Clustering is completed!")

        except Exception as msg:
            traceback.print_exc(file=sys.stdout)
            flash("Analysis ERROR: " + str(msg))
            print(msg)

    updated_page = render_template(
        "clustering.html",
        title=step_log["clustering"][1],
        clustering_form=form,
        flash_bool=True,
        step_log=step_log,
    )

    return updated_page


def run_psts(request, adata, step_log):
    """Performs psts analysis; must have performed clustering first."""
    # Creating the form with the clustering information #
    cluster_set = numpy.unique(adata.obs["clusters"].values)
    order = numpy.argsort([int(cluster) for cluster in cluster_set])
    cluster_set = cluster_set[order]

    options = ["Auto", "Spatial distance only", "Gene expression distance only"]

    from .utils import get_all_paths

    trajectory_set = get_all_paths(adata)

    PSTSForm = forms.getPSTSForm(trajectory_set, cluster_set, options)
    form = PSTSForm(request.form)

    step_log["psts_params"] = vhs.getData(form)
    print(step_log["psts_params"], file=sys.stdout)
    elements = list(step_log["psts_params"].keys())
    # order: pca_comps, SME bool, method, method_param
    element_values = list(step_log["psts_params"].values())

    if element_values[4] == "Auto":
        w = None
    elif element_values[4] == "Spatial distance only":
        w = 0
    else:
        w = 1

    if not form.validate_on_submit():
        flash_errors(form)

    elif type(adata) == type(None):
        flash("Need to load data first!")

    else:
        try:
            from stlearn.spatials.trajectory import set_root

            root_index = set_root(
                adata, use_label="clusters", cluster=str(element_values[0])
            )

            adata.uns["iroot"] = root_index

            print(root_index, file=sys.stdout, flush=True)

            # Performing the TI #
            print(element_values[3], file=sys.stdout, flush=True)

            node_order = element_values[3].split(" - ")

            st.spatial.trajectory.pseudotime(
                adata,
                eps=element_values[2],
                use_rep="X_pca",
                use_label="clusters",
                reverse=element_values[1],
            )
            print(node_order)
            st.spatial.trajectory.pseudotimespace_global(
                adata, use_label="clusters", list_clusters=node_order, w=w
            )

            st.pl.cluster_plot(
                adata,
                use_label="clusters",
                show_trajectories=True,
                list_clusters=node_order,
                show_subcluster=True,
            )
            savePlot("trajectory_inference.png")

            step_log["psts"][0] = True
            flash("Trajectory inference is completed!")

        except Exception as msg:
            traceback.print_exc(file=sys.stdout)
            flash("Analysis ERROR: " + str(msg))
            print(msg)

    updated_page = render_template(
        "psts.html",
        title=step_log["psts"][1],
        psts_form=form,
        flash_bool=True,
        step_log=step_log,
    )

    return updated_page


def run_dea(request, adata, step_log):

    list_labels = []

    for col in adata.obs.columns:
        if adata.obs[col].dtype.name == "category":
            if col != "sub_cluster_labels":
                list_labels.append(col)

    list_labels = numpy.array(list_labels)

    methods = numpy.array(["t-test", "t-test_overestim_var", "logreg", "wilcoxon"])

    DEAForm = forms.getDEAForm(list_labels, methods)
    form = DEAForm(request.form)

    step_log["dea_params"] = vhs.getData(form)
    print(step_log["dea_params"], file=sys.stdout)
    elements = list(step_log["dea_params"].keys())
    element_values = list(step_log["dea_params"].values())

    if not form.validate_on_submit():
        flash_errors(form)

    elif type(adata) == type(None):
        flash("Need to load data first!")

    else:
        try:

            sc.tl.rank_genes_groups(adata, element_values[0], method=element_values[1])

            step_log["dea"][0] = True
            flash("Differential expression analysis is completed!")

        except Exception as msg:
            traceback.print_exc(file=sys.stdout)
            flash("Analysis ERROR: " + str(msg))
            print(msg)

    updated_page = render_template(
        "dea.html",
        title=step_log["dea"][1],
        dea_form=form,
        flash_bool=True,
        step_log=step_log,
    )

    return updated_page

""" Junk code:

def run_cci(request, adata, step_log):
    #Performs CCI analysis.#

    form = CCIForm(request.form)

    if not form.validate_on_submit():
        flash_errors(form)

    elif type(adata) == type(None):
        flash("Need to load data first!")

    else:
        step_log["cci_params"] = vhs.getData(form)
        print(step_log["cci_params"], file=sys.stdout)
        elements = numpy.array(list(step_log["cci_params"].keys()))
        # order: cell_het file, neighbour_dist, L-R pairs, permutations
        element_values = list(step_log["cci_params"].values())

        cell_het = type(element_values[0]) != type(None)
        lrs, msg = vhs.getLR(element_values[2], adata.var_names)
        print(lrs, file=sys.stdout)

        if not form.validate_on_submit():
            flash_errors(form)

        elif type(adata) == type(None):
            flash("Need to load data first!")

        elif msg != "":
            flash(msg)

        else:
            # If we have been through the steps of:
            # 1. Choosing with/without cell heterogeneity information.
            # 2. Choosing with/without permutation testing.
            # 3. Filling out information.
            # We are now ready to take the form input & perform CCI analysis.
            
            try:
                adata.uns["lr"] = lrs
                st.tl.cci.lr(adata=adata, distance=element_values[1])
                st.pl.het_plot(adata, use_het="cci_lr", image_alpha=0.7)
                savePlot("cell_lr.png")

                dist, n_pairs = element_values[1], element_values[-1]
                if cell_het:  # Conduct with cell heterogeneity information #
                    print("We are using the cell heterogeneity!", file=sys.stdout)
                    # Adding the label transfer information #
                    st.add.labels(adata, element_values[0], sep="\t")
                    st.pl.cluster_plot(adata, use_label="predictions")
                    savePlot("label_transfer.png")  # saves to temp_pots

                    # Calculating cell heterogeneity #
                    st.tl.cci.het.count(adata, distance=dist, use_label="label_transfer")
                    st.pl.het_plot(adata, use_het="cci_het")
                    savePlot("cell_het.png")

                    # Merging with the lR values #
                    st.tl.cci.merge(adata, use_lr="cci_lr", use_het="cci_het")
                    st.pl.het_plot(adata, use_het="merged", cell_alpha=0.7)
                    savePlot("merged.png")

                if n_pairs != 0:  # Permutation testing #
                    st.tl.cci.permutation(
                        adata,
                        use_het="cci_het" if cell_het else None,
                        n_pairs=n_pairs,
                        distance=dist,
                    )
                    st.pl.het_plot(
                        adata,
                        cell_alpha=0.7,
                        use_het="merged_pvalues" if cell_het else "lr_pvalues",
                    )
                    savePlot("cci_rank-log10pvalues.png")

                    st.pl.het_plot(
                        adata,
                        cell_alpha=0.7,
                        use_het="merged_sign" if cell_het else "lr_sign",
                    )
                    savePlot("cci_rank-sig-log10pvalues.png")

                step_log["cci_rank"][0] = True

                flash("CCI analysis is completed!")

            except Exception as msg:
                traceback.print_exc(file=sys.stdout)
                flash("Analysis ERROR: " + str(msg))
                print(msg)

    updated_page = render_template(
        "cci.html",
        title=step_log["cci_rank"][1],
        cci_form=form,
        flash_bool=True,
        step_log=step_log,
    )

    return updated_page
"""

