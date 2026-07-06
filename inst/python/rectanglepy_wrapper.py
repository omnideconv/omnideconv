"""Wrapper script for rectanglepy integration with omnideconv.

Called via system2() from R with a dedicated conda environment (r-omnideconv-rectangle).
Data is exchanged through temporary files:
  - Single-cell data: .h5ad (AnnData, obs["label"] = cell type)
  - Bulk data: CSV, genes x samples (transposed to samples x genes on read)
  - Signature: pickle (RectangleSignatureResult)
  - Results: CSV, samples x cell_types
"""

import argparse
import pickle
import sys

import anndata as ad
import pandas as pd


CELL_TYPE_COL = "label"


def build_model(args):
    import rectanglepy

    adata = ad.read_h5ad(args.sc_h5ad)

    bulks = None
    if args.bulk_csv is not None:
        # omnideconv convention: genes x samples -> transpose to samples x genes
        bulks = pd.read_csv(args.bulk_csv, index_col=0).T

    kwargs = dict(
        optimize_cutoffs=args.optimize_cutoffs,
        p=args.p,
        lfc=args.lfc,
        n_cpus=args.n_cpus,
        gene_expression_threshold=args.gene_expression_threshold,
    )

    signature_result = rectanglepy.pp.build_rectangle_signatures(
        adata, CELL_TYPE_COL, bulks=bulks, **kwargs
    )

    with open(args.output_pickle, "wb") as f:
        pickle.dump(signature_result, f)


def deconvolute(args):
    import rectanglepy

    bulks = pd.read_csv(args.bulk_csv, index_col=0).T  # genes x samples -> samples x genes

    with open(args.signature_pickle, "rb") as f:
        signature_result = pickle.load(f)

    estimations, _ = rectanglepy.tl.deconvolution(
        signature_result, bulks,
        correct_mrna_bias=args.correct_mrna_bias,
        n_cpus=args.n_cpus,
    )

    estimations.to_csv(args.output_csv)


def get_signature_matrix(args):
    with open(args.signature_pickle, "rb") as f:
        signature_result = pickle.load(f)

    # Returns genes x cell_types DataFrame; omnideconv convention: rows=genes, cols=cell_types
    sig_matrix = signature_result.get_signature_matrix(include_mrna_bias=args.include_mrna_bias)
    sig_matrix.to_csv(args.output_csv)


def rectangle_all(args):
    import rectanglepy

    adata = ad.read_h5ad(args.sc_h5ad)
    bulks = pd.read_csv(args.bulk_csv, index_col=0).T  # genes x samples -> samples x genes

    estimations, _ = rectanglepy.rectangle(
        adata, bulks,
        cell_type_col=CELL_TYPE_COL,
        correct_mrna_bias=args.correct_mrna_bias,
        optimize_cutoffs=args.optimize_cutoffs,
        p=args.p,
        lfc=args.lfc,
        n_cpus=args.n_cpus,
        gene_expression_threshold=args.gene_expression_threshold,
    )

    estimations.to_csv(args.output_csv)


def main():
    parser = argparse.ArgumentParser(description="rectanglepy wrapper for omnideconv")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- build_model ---
    p_build = subparsers.add_parser("build_model", help="Build Rectangle signature model")
    p_build.add_argument("--sc_h5ad", required=True)
    p_build.add_argument("--output_pickle", required=True)
    p_build.add_argument("--bulk_csv", default=None)
    p_build.add_argument("--optimize_cutoffs", type=lambda x: x.lower() == "true", default=True)
    p_build.add_argument("--p", type=float, default=0.015)
    p_build.add_argument("--lfc", type=float, default=1.5)
    p_build.add_argument("--n_cpus", type=int, default=None)
    p_build.add_argument("--gene_expression_threshold", type=float, default=0.5)

    # --- deconvolute ---
    p_deconv = subparsers.add_parser("deconvolute", help="Deconvolute bulk using pre-built signature")
    p_deconv.add_argument("--bulk_csv", required=True)
    p_deconv.add_argument("--signature_pickle", required=True)
    p_deconv.add_argument("--output_csv", required=True)
    p_deconv.add_argument("--correct_mrna_bias", type=lambda x: x.lower() == "true", default=True)
    p_deconv.add_argument("--n_cpus", type=int, default=None)

    # --- get_signature_matrix ---
    p_sig = subparsers.add_parser("get_signature_matrix", help="Extract signature matrix from pickle")
    p_sig.add_argument("--signature_pickle", required=True)
    p_sig.add_argument("--output_csv", required=True)
    p_sig.add_argument("--include_mrna_bias", type=lambda x: x.lower() == "true", default=True)

    # --- rectangle (all-in-one) ---
    p_rect = subparsers.add_parser("rectangle", help="Build signature and deconvolute in one step")
    p_rect.add_argument("--sc_h5ad", required=True)
    p_rect.add_argument("--bulk_csv", required=True)
    p_rect.add_argument("--output_csv", required=True)
    p_rect.add_argument("--correct_mrna_bias", type=lambda x: x.lower() == "true", default=True)
    p_rect.add_argument("--optimize_cutoffs", type=lambda x: x.lower() == "true", default=True)
    p_rect.add_argument("--p", type=float, default=0.015)
    p_rect.add_argument("--lfc", type=float, default=1.5)
    p_rect.add_argument("--n_cpus", type=int, default=None)
    p_rect.add_argument("--gene_expression_threshold", type=float, default=0.5)

    args = parser.parse_args()

    if args.command == "build_model":
        build_model(args)
    elif args.command == "deconvolute":
        deconvolute(args)
    elif args.command == "get_signature_matrix":
        get_signature_matrix(args)
    elif args.command == "rectangle":
        rectangle_all(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
