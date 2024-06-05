#!/usr/bin/env python3

import plotly.graph_objects as go
import pandas as pd
import argparse


def load_data(data_filename, karyotype_ref_file_name, karyotype_target_filename):
    data = pd.read_csv(data_filename, sep="\t", header=None)
    data.columns = [
        "ref",
        "ref_start",
        "ref_stop",
        "target",
        "target_start",
        "target_stop",
        "color",
    ]

    karyotype_ref = pd.read_csv(karyotype_ref_file_name, sep="\t", header=None)
    karyotype_ref.columns = ["chr", "name", "name1", "zero", "size", "color"]
    karyotype_target = pd.read_csv(karyotype_target_filename, sep="\t", header=None)
    karyotype_target.columns = ["chr", "name", "name1", "zero", "size", "color"]
    return data, karyotype_ref, karyotype_target


def insert_offsets_and_return_dict(karyotype_df):
    karyotype_df["offset"] = karyotype_df["size"].cumsum().shift(fill_value=0)
    return dict(zip(karyotype_df["name"], karyotype_df["offset"]))


def insert_midpoints_and_return_dict(karyotype_df):
    karyotype_df["midpoint"] = (karyotype_df["size"] / 2.0) + karyotype_df["offset"]
    return dict(zip(karyotype_df["name"], karyotype_df["midpoint"]))


def insert_figure_data(
    data, offsets_ref, offsets_target, midpoints_ref, midpoints_target, fig
):
    for index, row in data.iterrows():
        x_ = [i + offsets_ref[row["ref"]] for i in [row["ref_start"], row["ref_stop"]]]

        y_ = [
            j + offsets_target[row["target"]]
            for j in [row["target_start"], row["target_stop"]]
        ]

        midpoint_x = midpoints_ref[row["ref"]]
        midpoint_y = midpoints_target[row["target"]]

        fig.add_trace(
            go.Scatter(
                x=x_,
                y=y_,
                mode="lines",
                line=dict(
                    color=f"rgba{(row['color'].replace('color=', '').replace('0.5)', '1)'))}",
                    width=2,
                ),
                name=f"{index}: {row['target']}:{row['ref']}",
                legendgroup=f"{row['target']}:{row['ref']}",
                legendgrouptitle=dict(text=f"{row['target']}:{row['ref']}"),
            )
        )

        fig.add_trace(
            go.Scatter(
                x=[midpoint_x],
                y=[0],
                xaxis="x2",
                line=dict(color="#ffffff"),
                showlegend=False,
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[0],
                y=[midpoint_y],
                yaxis="y2",
                line=dict(color="#ffffff"),
                showlegend=False,
            )
        )


def format_figure(
    karyotype_ref,
    karyotype_target,
    offsets_ref,
    offsets_target,
    midpoints_ref,
    midpoints_target,
    fig,
):

    xaxis_range = [0, list(offsets_ref.values())[-1] + list(karyotype_ref["size"])[-1]]
    yaxis_range = [
        0,
        list(offsets_target.values())[-1] + list(karyotype_target["size"])[-1],
    ]

    fig.update_layout(
        xaxis=dict(
            range=xaxis_range,
        ),
        yaxis=dict(range=yaxis_range),
        hovermode="closest",
        xaxis2=dict(
            overlaying="x",
            scaleanchor="x1",
            range=xaxis_range,
            tickmode="array",
            tickvals=list(offsets_ref.values()),
            ticktext=list(offsets_ref.keys()),
            side="top",
        ),
        yaxis2=dict(
            overlaying="y",
            scaleanchor="y1",
            range=yaxis_range,
            tickmode="array",
            tickvals=list(offsets_target.values()),
            ticktext=list(offsets_target.keys()),
            autoshift=True,
            anchor="free",
        ),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        showlegend=False,
    )

    fig.update_xaxes(
        showgrid=True,
        zeroline=False,
        gridcolor="rgba(0, 0, 0, 0.1)",
        griddash="dashdot",
    )
    fig.update_yaxes(
        showgrid=True,
        zeroline=False,
        gridcolor="rgba(0, 0, 0, 0.1)",
        griddash="dashdot",
    )
    fig["layout"]["yaxis1"]["showgrid"] = False
    fig["layout"]["xaxis1"]["showgrid"] = False


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Create a linear synteny plot from `nucmer/dnadiff/circos bundlelinks` bundles",
        epilog="Author: Usman Rashid",
    )
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.1")
    parser.add_argument(
        "--output",
        type=str,
        default="synteny_plot.html",
        required=False,
        help="Output filename",
    )
    parser.add_argument("bundlelinks", type=argparse.FileType("r"))
    parser.add_argument("karyotype_ref", type=argparse.FileType("r"))
    parser.add_argument("karyotype_target", type=argparse.FileType("r"))

    args = parser.parse_args()

    data_filename = args.bundlelinks
    karyotype_ref_file_name = args.karyotype_ref
    karyotype_target_filename = args.karyotype_target
    output_filename = args.output

    data, karyotype_ref, karyotype_target = load_data(
        data_filename, karyotype_ref_file_name, karyotype_target_filename
    )

    offsets_ref = insert_offsets_and_return_dict(karyotype_ref)
    offsets_target = insert_offsets_and_return_dict(karyotype_target)

    midpoints_ref = insert_midpoints_and_return_dict(karyotype_ref)
    midpoints_target = insert_midpoints_and_return_dict(karyotype_target)

    fig = go.Figure()

    insert_figure_data(
        data, offsets_ref, offsets_target, midpoints_ref, midpoints_target, fig
    )

    format_figure(
        karyotype_ref,
        karyotype_target,
        offsets_ref,
        offsets_target,
        midpoints_ref,
        midpoints_target,
        fig,
    )

    fig.write_html(output_filename)
