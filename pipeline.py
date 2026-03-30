import os
import subprocess
from typing import List
import pandas as pd
from Bio import SeqIO
import plotly.graph_objects as go
import plotly.express as px


# -----------------------------
# CONFIG
# -----------------------------
class Config:
    def __init__(
        self,
        ref_proteins,
        ref_gff,
        test_proteins,
        test_gff,
        gene_name,
        window_size,
        outdir="output"
    ):
        self.ref_proteins = ref_proteins
        self.ref_gff = ref_gff
        self.test_proteins = test_proteins
        self.test_gff = test_gff
        self.gene_name = gene_name
        self.window_size = window_size
        self.outdir = outdir
        os.makedirs(outdir, exist_ok=True)


# -----------------------------
# GFF PARSING (CDS → protein)
# -----------------------------
def parse_gff(gff_file: str) -> pd.DataFrame:
    rows = []

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            seqid, _, feature, start, end, _, strand, _, attr = parts
            if feature != "CDS":
                continue

            attrs = dict(
                field.split("=", 1) for field in attr.split(";") if "=" in field
            )

            protein_id = attrs.get("protein_id") or attrs.get("Name")
            if protein_id is None:
                continue

            rows.append({
                "seqid": seqid,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "gene_id": protein_id
            })

    df = pd.DataFrame(rows)

    df = (
        df.groupby(["seqid", "gene_id", "strand"], as_index=False)
        .agg(start=("start", "min"), end=("end", "max"))
        .sort_values(["seqid", "start"])
        .reset_index(drop=True)
    )

    return df


# -----------------------------
# WINDOW SELECTION
# -----------------------------
def get_gene_window(df, gene_name, window):
    idx = df[df["gene_id"] == gene_name].index
    if idx.empty:
        raise ValueError(f"{gene_name} not found")

    i = idx[0]
    return df.iloc[max(0, i - window): i + window + 1]


# -----------------------------
# FASTA EXTRACTION
# -----------------------------
def extract_proteins(fasta, ids, out):
    SeqIO.write(
        [r for r in SeqIO.parse(fasta, "fasta") if r.id in ids],
        out,
        "fasta"
    )


# -----------------------------
# DIAMOND
# -----------------------------
def run_diamond(query, db_fa, out):
    db = os.path.join(os.path.dirname(out), "diamond_db")
    subprocess.run(["diamond", "makedb", "--in", db_fa, "-d", db], check=True)
    subprocess.run([
        "diamond", "blastp",
        "-q", query,
        "-d", db,
        "-o", out,
        "--max-target-seqs", "1",
        "--evalue", "1e-50"
    ], check=True)


def parse_hits(hit_file):
    cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore"
    ]
    return pd.read_csv(hit_file, sep="\t", names=cols)


# -----------------------------
# SYNTENY PLOT
# -----------------------------
def plot_synteny(
    ref_df,
    test_df,
    hits_df,
    focal_gene,
    focal_test_gene,
    outprefix
):
    chromosomes = test_df["seqid"].unique()
    chr_to_y = {c: i for i, c in enumerate(chromosomes)}
    ref_y = len(chromosomes)

    colors = px.colors.qualitative.Dark24
    chr_color = {c: colors[i % len(colors)] for i, c in enumerate(chromosomes)}

    fig = go.Figure()

    # -----------------------------
    # Gene spans (rectangles)
    # -----------------------------
    def draw_gene(row, y, color, highlight=False):
        fig.add_shape(
            type="rect",
            x0=row.start,
            x1=row.end,
            y0=y - 0.2,
            y1=y + 0.2,
            line=dict(width=2 if highlight else 1),
            fillcolor=color,
            opacity=0.9
        )

    # Reference genes
    for _, r in ref_df.iterrows():
        draw_gene(
            r,
            ref_y,
            "gold" if r.gene_id == focal_gene else "black",
            highlight=r.gene_id == focal_gene
        )

        fig.add_trace(go.Scatter(
            x=[(r.start + r.end) / 2],
            y=[ref_y],
            mode="markers",
            marker=dict(
                symbol="triangle-right" if r.strand == "+" else "triangle-left",
                size=10
            ),
            hovertext=(
                f"Gene: {r.gene_id}<br>"
                f"Chr: {r.seqid}<br>"
                f"Start: {r.start}<br>"
                f"End: {r.end}<br>"
                f"Strand: {r.strand}"
            ),
            showlegend=False
        ))

    # Test genes
    for _, r in test_df.iterrows():
        y = chr_to_y[r.seqid]
        color = "red" if r.gene_id == focal_test_gene else chr_color[r.seqid]

        draw_gene(
            r,
            y,
            color,
            highlight=r.gene_id == focal_test_gene
        )

        fig.add_trace(go.Scatter(
            x=[(r.start + r.end) / 2],
            y=[y],
            mode="markers",
            marker=dict(
                symbol="triangle-right" if r.strand == "+" else "triangle-left",
                size=9
            ),
            hovertext=(
                f"Gene: {r.gene_id}<br>"
                f"Chr: {r.seqid}<br>"
                f"Start: {r.start}<br>"
                f"End: {r.end}<br>"
                f"Strand: {r.strand}"
            ),
            showlegend=False
        ))

    # -----------------------------
    # Ribbons (inversions highlighted)
    # -----------------------------
    max_bitscore = hits_df["bitscore"].max()

    for _, h in hits_df.iterrows():
        if h.qseqid not in ref_df.gene_id.values:
            continue
        if h.sseqid not in test_df.gene_id.values:
            continue

        r = ref_df[ref_df.gene_id == h.qseqid].iloc[0]
        t = test_df[test_df.gene_id == h.sseqid].iloc[0]

        inverted = r.strand != t.strand
        is_focal = (r.gene_id == focal_gene and t.gene_id == focal_test_gene)

        fig.add_trace(go.Scatter(
            x=[(r.start + r.end) / 2, (t.start + t.end) / 2],
            y=[ref_y, chr_to_y[t.seqid]],
            mode="lines",
            line=dict(
                width=2 + 4 * h.bitscore / max_bitscore,
                color="red" if is_focal else (
                    "purple" if inverted else chr_color[t.seqid]
                ),
                dash="dash" if inverted else "solid"
            ),
            opacity=0.7,
            hovertext=(
                f"Query: {h.qseqid}<br>"
                f"Subject: {h.sseqid}<br>"
                f"Bitscore: {h.bitscore}<br>"
                f"Inversion: {inverted}"
            ),
            showlegend=False
        ))

    fig.update_layout(
        title="Synteny Plot (gene spans, inversions, bp scale)",
        xaxis_title="Genomic position (bp)",
        yaxis=dict(
            tickvals=list(chr_to_y.values()) + [ref_y],
            ticktext=list(chr_to_y.keys()) + ["Reference"]
        ),
        template="simple_white"
    )

    fig.write_html(f"{outprefix}.html")
    fig.write_image(f"{outprefix}.svg")
    fig.write_image(f"{outprefix}.pdf")


# -----------------------------
# MAIN PIPELINE
# -----------------------------
def run_pipeline(cfg: Config):
    ref_gff = parse_gff(cfg.ref_gff)
    test_gff = parse_gff(cfg.test_gff)

    ref_window = get_gene_window(ref_gff, cfg.gene_name, cfg.window_size)

    query_fa = os.path.join(cfg.outdir, "query.fa")
    extract_proteins(cfg.ref_proteins, ref_window.gene_id.tolist(), query_fa)

    hits_file = os.path.join(cfg.outdir, "diamond.tsv")
    run_diamond(query_fa, cfg.test_proteins, hits_file)
    hits_df = parse_hits(hits_file)

    focal_hit = hits_df[hits_df.qseqid == cfg.gene_name]
    focal_test_gene = focal_hit.sseqid.iloc[0] if not focal_hit.empty else None

    test_subset = test_gff[test_gff.gene_id.isin(hits_df.sseqid)]

    plot_synteny(
        ref_window,
        test_subset,
        hits_df,
        cfg.gene_name,
        focal_test_gene,
        os.path.join(cfg.outdir, "synteny")
    )


if __name__ == "__main__":
    cfg = Config(
        ref_proteins="ref.fa",
        ref_gff="ref.gff",
        test_proteins="test.fa",
        test_gff="test.gff",
        gene_name="AGB95346.1",
        window_size=20
    )
    run_pipeline(cfg)