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

from scipy.spatial.distance import cosine

# Creating the forms using a class generator #
PreprocessForm = forms.getPreprocessForm()
CCIForm = forms.getCCIForm()
ClusterForm = forms.getClusterForm()


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
        flash("Preprocessing completed!")

    updated_page = render_template(
        "preprocessing.html",
        title=step_log["preprocessed"][1],
        preprocess_form=form,
        flash_bool=step_log["preprocessed"][0],
        step_log=step_log,
    )

    return updated_page


def run_cci(request, adata, step_log):
    """Performs CCI analysis."""

    form = CCIForm(request.form)

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
        """If we have been through the steps of:
        1. Choosing with/without cell heterogeneity information.
        2. Choosing with/without permutation testing.
        3. Filling out information.
        We are now ready to take the form input & perform CCI analysis.
        """
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
                savePlot("cci-log10pvalues.png")

                st.pl.het_plot(
                    adata,
                    cell_alpha=0.7,
                    use_het="merged_sign" if cell_het else "lr_sign",
                )
                savePlot("cci-sig-log10pvalues.png")

            step_log["cci"][0] = True

            flash("CCI analysis completed!")

        except Exception as msg:
            traceback.print_exc(file=sys.stdout)
            flash("Analysis ERROR: " + str(msg))
            print(msg)

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
                st.tl.clustering.kmeans(
                    adata, n_clusters=param, use_data="X_pca", key_added="clusters"
                )

            else:  # Louvain
                param = element_values[4]
                st.pp.neighbors(adata, n_neighbors=element_values[5], use_rep="X_pca")
                st.tl.clustering.louvain(adata, resolution=param, key_added="clusters")

            st.pp.neighbors(adata, n_neighbors=element_values[5], use_rep="X_pca")
            sc.tl.paga(adata, groups="clusters")

            st.pl.cluster_plot(adata, use_label="clusters")
            savePlot("clustering.png")

            step_log["clustering"][0] = True
            flash("Clustering completed!")

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

    from .utils import get_all_paths

    trajectory_set = get_all_paths(adata)

    PSTSForm = forms.getPSTSForm(trajectory_set, cluster_set)
    form = PSTSForm(request.form)

    step_log["psts_params"] = vhs.getData(form)
    print(step_log["psts_params"], file=sys.stdout)
    elements = list(step_log["psts_params"].keys())
    # order: pca_comps, SME bool, method, method_param
    element_values = list(step_log["psts_params"].values())

    if not form.validate_on_submit():
        flash_errors(form)

    elif type(adata) == type(None):
        flash("Need to load data first!")

    else:
        try:
            # Getting root spot position as one closest to median location for
            # cluster #
            # cluster_spots = adata.obs['clusters'].values == element_values[0]
            # spot_locs = adata.obs.iloc[:,1:3].loc[cluster_spots,:]
            # median = [numpy.median(spot_locs.values[:,0]),
            # 		  numpy.median(spot_locs.values[:,1])]
            # med_dists = numpy.apply_along_axis(cosine, 1,
            # 								   spot_locs.values, median)
            # min_dist = min(med_dists)
            # root_name = spot_locs.index.values[med_dists==min_dist][0]
            from stlearn.spatials.trajectory import set_root

            root_index = set_root(
                adata, use_label="clusters", cluster=str(element_values[0])
            )

            adata.uns["iroot"] = root_index

            print(root_index, file=sys.stdout, flush=True)

            # Performing the TI #
            print(element_values[-1], file=sys.stdout, flush=True)

            node_order = element_values[-1].split(" - ")

            st.spatial.trajectory.pseudotime(
                adata,
                eps=element_values[2],
                use_rep="X_pca",
                use_label="clusters",
                reverse=element_values[1],
            )
            print(node_order)
            st.spatial.trajectory.pseudotimespace_global(
                adata, use_label="clusters", list_clusters=node_order
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
            flash("Trajectory inference completed!")

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
