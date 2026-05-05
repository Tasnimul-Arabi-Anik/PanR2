import math
import os

import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
from scipy.stats import fisher_exact

from panr2.viz import crowding_warning_row, label_map, label_warning_rows, plot_style_config


DATABASE_PREFIXES = {
    "amr": "AMR",
    "vfdb": "VFDB",
    "plasmidfinder": "PLASMID",
    "mobileelementfinder": "MGE",
    "isfinder": "IS",
    "integronfinder": "INTEGRON",
    "iceberg": "ICE",
    "mlst": "MLST",
    "defensefinder": "DEFENSE",
    "prophage": "PROPHAGE",
}

DATABASE_COLORS = {
    "AMR": "#C43B3B",
    "VFDB": "#7A4EB3",
    "PLASMID": "#2F6CBE",
    "MGE": "#3E8E41",
    "IS": "#67A953",
    "INTEGRON": "#008C8C",
    "ICE": "#8A6A3D",
    "MLST": "#4B5563",
    "DEFENSE": "#9467BD",
    "PROPHAGE": "#D97706",
}

METADATA_COLUMNS = [
    "Geographic Location",
    "Continent",
    "Subcontinent",
    "Collection Date",
    "Collection Year",
    "Host",
    "Host_SD",
    "Host_Rank",
    "Host_Genus",
    "Host_Species",
    "Host_Class",
    "Host_Order",
    "Host_Anatomical_Site_SD",
    "Host_Age_Group_SD",
    "Host_Production_Context_SD",
    "Isolation Source",
    "Isolation_Source_SD",
    "Sample_Type_SD",
    "Environment_Medium_SD",
    "Environment_Broad_Scale_SD",
    "Environment_Local_Scale_SD",
    "Organism Name",
    "Genus",
    "Species",
    "TaxID",
]


def _safe_feature(value):
    text = str(value).strip()
    return text if text and text.lower() not in {"nan", "none", "0"} else "unknown_feature"


def _feature_database(feature_name):
    return str(feature_name).split(":", 1)[0] if ":" in str(feature_name) else "UNKNOWN"


def _fdr_bh(p_values):
    """Benjamini-Hochberg correction implemented without an extra dependency."""
    indexed = [(i, p) for i, p in enumerate(p_values) if pd.notna(p)]
    q_values = [float("nan")] * len(p_values)
    if not indexed:
        return q_values
    indexed.sort(key=lambda item: item[1])
    m = len(indexed)
    running = 1.0
    for rank, (original_index, p_value) in reversed(list(enumerate(indexed, start=1))):
        running = min(running, (p_value * m) / rank)
        q_values[original_index] = min(running, 1.0)
    return q_values


def _sample_metadata(df):
    columns = [col for col in ["Assembly Accession", "Assembly BioSample Accession"] + METADATA_COLUMNS if col in df.columns]
    if "Assembly Accession" not in columns:
        return pd.DataFrame()
    meta = df[columns].drop_duplicates(subset=["Assembly Accession"]).copy()
    if "Collection Date" in meta.columns:
        meta["Collection Year"] = pd.to_numeric(meta["Collection Date"], errors="coerce").astype("Int64").astype(str).replace("<NA>", pd.NA)
    return meta


def _matrix_from_amr_tidy(tidy_df):
    if tidy_df is None or tidy_df.empty or "Assembly Accession" not in tidy_df.columns or "Gene" not in tidy_df.columns:
        return pd.DataFrame(), {}, pd.DataFrame()
    df = tidy_df.copy()
    presence_col = "Presence" if "Presence" in df.columns else "presence"
    if presence_col not in df.columns:
        return pd.DataFrame(), {}, _sample_metadata(df)
    df[presence_col] = pd.to_numeric(df[presence_col], errors="coerce").fillna(0).astype(int)
    df = df[df[presence_col] > 0].copy()
    if df.empty:
        return pd.DataFrame(), {}, _sample_metadata(tidy_df)
    df["feature_name"] = "AMR:" + df["Gene"].map(_safe_feature)
    matrix = (
        df[["Assembly Accession", "feature_name"]]
        .dropna()
        .drop_duplicates()
        .assign(value=1)
        .pivot_table(index="Assembly Accession", columns="feature_name", values="value", fill_value=0, aggfunc="max")
        .astype(int)
    )
    feature_map = {feature: "AMR" for feature in matrix.columns}
    return matrix, feature_map, _sample_metadata(tidy_df)


def _matrix_from_optional_tidy(feature_type, tidy_path):
    if not tidy_path or not os.path.exists(tidy_path):
        return pd.DataFrame(), {}, pd.DataFrame()
    df = pd.read_csv(tidy_path)
    if df.empty or "Assembly Accession" not in df.columns or "feature_id" not in df.columns:
        return pd.DataFrame(), {}, _sample_metadata(df)
    if "presence" not in df.columns:
        return pd.DataFrame(), {}, _sample_metadata(df)
    prefix = DATABASE_PREFIXES.get(feature_type, feature_type.upper())
    df["presence"] = pd.to_numeric(df["presence"], errors="coerce").fillna(0).astype(int)
    df = df[df["presence"] > 0].copy()
    if df.empty:
        return pd.DataFrame(), {}, _sample_metadata(pd.read_csv(tidy_path))
    df["feature_name"] = prefix + ":" + df["feature_id"].map(_safe_feature)
    matrix = (
        df[["Assembly Accession", "feature_name"]]
        .dropna()
        .drop_duplicates()
        .assign(value=1)
        .pivot_table(index="Assembly Accession", columns="feature_name", values="value", fill_value=0, aggfunc="max")
        .astype(int)
    )
    feature_map = {feature: prefix for feature in matrix.columns}
    return matrix, feature_map, _sample_metadata(pd.read_csv(tidy_path))


def build_unified_feature_matrix(amr_tidy_df, feature_outputs):
    matrices = []
    feature_map = {}
    metadata_frames = []

    amr_matrix, amr_map, amr_meta = _matrix_from_amr_tidy(amr_tidy_df)
    if not amr_matrix.empty:
        matrices.append(amr_matrix)
        feature_map.update(amr_map)
    if not amr_meta.empty:
        metadata_frames.append(amr_meta)

    for feature_type, outputs in sorted((feature_outputs or {}).items()):
        matrix, db_map, meta = _matrix_from_optional_tidy(feature_type, outputs.get("tidy"))
        if not matrix.empty:
            matrices.append(matrix)
            feature_map.update(db_map)
        if not meta.empty:
            metadata_frames.append(meta)

    if matrices:
        feature_matrix = pd.concat(matrices, axis=1).fillna(0).astype(int)
        feature_matrix = feature_matrix.groupby(level=0).max()
    else:
        feature_matrix = pd.DataFrame()

    metadata = pd.DataFrame()
    if metadata_frames:
        metadata = pd.concat(metadata_frames, ignore_index=True)
        metadata = metadata.drop_duplicates(subset=["Assembly Accession"], keep="first")
    if not feature_matrix.empty and not metadata.empty and "Assembly Accession" in metadata.columns:
        all_samples = metadata["Assembly Accession"].dropna().astype(str).drop_duplicates().tolist()
        feature_matrix = feature_matrix.reindex(all_samples, fill_value=0).astype(int)
    return feature_matrix, feature_map, metadata


def _select_features_for_pairwise(matrix, max_features):
    if matrix.empty:
        return matrix, []
    prevalence = matrix.sum(axis=0).sort_values(ascending=False)
    if max_features and max_features > 0 and len(prevalence) > max_features:
        selected = prevalence.head(max_features).index.tolist()
        return matrix[selected], selected
    return matrix, prevalence.index.tolist()


def _association_rows(matrix, feature_map):
    rows = []
    features = list(matrix.columns)
    total = int(matrix.shape[0])
    for i, first in enumerate(features):
        first_values = matrix[first].astype(int)
        first_db = feature_map.get(first, _feature_database(first))
        for second in features[i + 1:]:
            second_values = matrix[second].astype(int)
            second_db = feature_map.get(second, _feature_database(second))
            n11 = int(((first_values == 1) & (second_values == 1)).sum())
            n10 = int(((first_values == 1) & (second_values == 0)).sum())
            n01 = int(((first_values == 0) & (second_values == 1)).sum())
            n00 = int(((first_values == 0) & (second_values == 0)).sum())
            denom = math.sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
            phi = ((n11 * n00) - (n10 * n01)) / denom if denom else 0.0
            jaccard_denominator = n11 + n10 + n01
            try:
                odds_ratio, p_value = fisher_exact([[n11, n10], [n01, n00]])
            except Exception:
                odds_ratio, p_value = float("nan"), float("nan")
            rows.append({
                "feature_1": first,
                "feature_2": second,
                "database_1": first_db,
                "database_2": second_db,
                "association_scope": "cross_database" if first_db != second_db else "within_database",
                "n11": n11,
                "n10": n10,
                "n01": n01,
                "n00": n00,
                "cooccurring_samples": n11,
                "cooccurrence_percentage": (n11 / max(total, 1)) * 100,
                "feature_1_samples": n11 + n10,
                "feature_2_samples": n11 + n01,
                "jaccard_index": n11 / jaccard_denominator if jaccard_denominator else 0.0,
                "phi_coefficient": phi,
                "odds_ratio": odds_ratio,
                "p_value": p_value,
            })
    pair_df = pd.DataFrame(rows)
    if pair_df.empty:
        return pair_df
    pair_df["q_value"] = _fdr_bh(pair_df["p_value"].tolist())
    return pair_df.sort_values(
        ["association_scope", "q_value", "cooccurring_samples", "jaccard_index"],
        ascending=[True, True, False, False],
    )


def _write_pair_matrices(matrix, analysis_dir):
    cooc_path = os.path.join(analysis_dir, "cross_database_cooccurrence_matrix.csv")
    jaccard_path = os.path.join(analysis_dir, "cross_database_jaccard_matrix.csv")
    phi_path = os.path.join(analysis_dir, "cross_database_phi_correlation_matrix.csv")
    if matrix.empty:
        pd.DataFrame().to_csv(cooc_path)
        pd.DataFrame().to_csv(jaccard_path)
        pd.DataFrame().to_csv(phi_path)
        return cooc_path, jaccard_path, phi_path

    cooc = matrix.T.dot(matrix).astype(int)
    cooc.to_csv(cooc_path)

    features = list(matrix.columns)
    jaccard = pd.DataFrame(0.0, index=features, columns=features)
    phi = pd.DataFrame(0.0, index=features, columns=features)
    for first in features:
        first_values = matrix[first].astype(int)
        for second in features:
            second_values = matrix[second].astype(int)
            n11 = int(((first_values == 1) & (second_values == 1)).sum())
            n10 = int(((first_values == 1) & (second_values == 0)).sum())
            n01 = int(((first_values == 0) & (second_values == 1)).sum())
            n00 = int(((first_values == 0) & (second_values == 0)).sum())
            jaccard.loc[first, second] = n11 / (n11 + n10 + n01) if (n11 + n10 + n01) else 0.0
            denom = math.sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
            phi.loc[first, second] = ((n11 * n00) - (n10 * n01)) / denom if denom else 0.0
    jaccard.to_csv(jaccard_path)
    phi.to_csv(phi_path)
    return cooc_path, jaccard_path, phi_path


def _write_specific_associations(pair_df, analysis_dir):
    specs = {
        "amr_mge_associations": ({"AMR"}, {"MGE", "IS", "INTEGRON", "ICE"}),
        "amr_plasmid_associations": ({"AMR"}, {"PLASMID"}),
        "amr_integron_associations": ({"AMR"}, {"INTEGRON"}),
        "amr_virulence_associations": ({"AMR"}, {"VFDB"}),
        "amr_defense_associations": ({"AMR"}, {"DEFENSE"}),
        "amr_prophage_associations": ({"AMR"}, {"PROPHAGE"}),
        "plasmid_mge_associations": ({"PLASMID"}, {"MGE", "IS", "INTEGRON", "ICE"}),
        "defense_mge_associations": ({"DEFENSE"}, {"MGE", "IS", "INTEGRON", "ICE", "PROPHAGE"}),
        "prophage_mge_associations": ({"PROPHAGE"}, {"MGE", "IS", "INTEGRON", "ICE"}),
    }
    paths = {}
    for name, (left, right) in specs.items():
        path = os.path.join(analysis_dir, f"{name}.csv")
        if pair_df.empty:
            pd.DataFrame().to_csv(path, index=False)
        else:
            mask = (
                pair_df["database_1"].isin(left) & pair_df["database_2"].isin(right)
            ) | (
                pair_df["database_1"].isin(right) & pair_df["database_2"].isin(left)
            )
            pair_df[mask].sort_values(["q_value", "cooccurring_samples", "jaccard_index"], ascending=[True, False, False]).to_csv(path, index=False)
        paths[name] = path
    return paths


def _write_integrated_burden(matrix, metadata, analysis_dir, figures_dir, fig_format):
    path = os.path.join(analysis_dir, "sample_integrated_feature_burden.csv")
    if matrix.empty:
        pd.DataFrame().to_csv(path, index=False)
        return {"sample_integrated_feature_burden": path}

    rows = pd.DataFrame({"Assembly Accession": matrix.index})
    for db in sorted(set(_feature_database(col) for col in matrix.columns)):
        cols = [col for col in matrix.columns if _feature_database(col) == db]
        rows[f"{db.lower()}_feature_count"] = matrix[cols].sum(axis=1).astype(int).values
    mobile_cols = [col for col in ["mge_feature_count", "is_feature_count", "integron_feature_count", "ice_feature_count", "prophage_feature_count"] if col in rows.columns]
    rows["total_mobileome_count"] = rows[mobile_cols].sum(axis=1).astype(int) if mobile_cols else 0
    count_cols = [col for col in rows.columns if col.endswith("_feature_count")]
    rows["total_feature_count"] = rows[count_cols].sum(axis=1).astype(int) if count_cols else 0
    if "amr_feature_count" in rows.columns:
        rows["mobility_associated_amr_score"] = rows["amr_feature_count"] * (rows["total_mobileome_count"] > 0).astype(int)
    else:
        rows["mobility_associated_amr_score"] = 0
    if metadata is not None and not metadata.empty:
        rows = metadata.merge(rows, on="Assembly Accession", how="right")
    rows.to_csv(path, index=False)

    plot_paths = {"sample_integrated_feature_burden": path}
    html_dir = os.path.join(figures_dir, "html_files")
    os.makedirs(html_dir, exist_ok=True)
    scatter_specs = [
        ("amr_feature_count", "total_mobileome_count", "AMR_vs_mobileome_burden"),
        ("amr_feature_count", "plasmid_feature_count", "AMR_vs_plasmid_burden"),
        ("amr_feature_count", "vfdb_feature_count", "AMR_vs_virulence_burden"),
        ("amr_feature_count", "defense_feature_count", "AMR_vs_defense_burden"),
        ("amr_feature_count", "prophage_feature_count", "AMR_vs_prophage_burden"),
    ]
    for x_col, y_col, name in scatter_specs:
        if x_col not in rows.columns or y_col not in rows.columns:
            continue
        plt.figure(figsize=(7, 5))
        hue = "Continent" if "Continent" in rows.columns else None
        sns.scatterplot(data=rows, x=x_col, y=y_col, hue=hue, s=55)
        plt.xlabel(x_col.replace("_", " ").title())
        plt.ylabel(y_col.replace("_", " ").title())
        plt.title(name.replace("_", " "))
        plt.tight_layout()
        png_path = os.path.join(figures_dir, f"{name}.{fig_format}")
        plt.savefig(png_path, dpi=300, bbox_inches="tight", format=fig_format)
        plt.close()
        html_path = os.path.join(html_dir, f"{name}.html")
        hover_cols = [col for col in ["Assembly Accession", "Geographic Location", "Continent", "Collection Date", "Organism Name"] if col in rows.columns]
        fig = px.scatter(rows, x=x_col, y=y_col, color="Continent" if "Continent" in rows.columns else None,
                         hover_data=hover_cols, title=name.replace("_", " "))
        fig.write_html(html_path)
        plot_paths[f"{name}_plot"] = png_path
        plot_paths[f"{name}_html"] = html_path
    return plot_paths


def _write_heatmaps(matrix, phi_path, figures_dir, fig_format, top_n, style_config):
    plot_paths = {}
    warnings = []
    if matrix.empty:
        return plot_paths, warnings
    prevalence = matrix.sum(axis=0).sort_values(ascending=False)
    heatmap_limit = min(top_n, style_config.get("heatmap_max_features", top_n))
    warning = crowding_warning_row("integrated_feature_presence_heatmap", len(prevalence), heatmap_limit, item="features")
    if warning:
        warnings.append(warning)
    top_features = prevalence.head(heatmap_limit).index.tolist()
    warnings.extend(label_warning_rows(top_features, style_config.get("label_max_length", 55), "cross_database_heatmaps"))
    html_dir = os.path.join(figures_dir, "html_files")
    os.makedirs(html_dir, exist_ok=True)
    labels = label_map(
        top_features,
        max_length=style_config.get("label_max_length", 55),
        wrap_width=style_config.get("wrap_width", 18),
    )

    if len(top_features) >= 2 and os.path.exists(phi_path):
        phi = pd.read_csv(phi_path, index_col=0)
        plot_matrix = phi.loc[top_features, top_features].rename(index=labels, columns=labels)
        size = max(6, min(18, 0.36 * len(top_features) + 4))
        plt.figure(figsize=(size, size))
        sns.heatmap(plot_matrix, cmap="vlag", center=0, vmin=-1, vmax=1, linewidths=0.2, linecolor="white", cbar_kws={"label": "Phi coefficient"})
        plt.title("Global cross-database feature association")
        plt.xticks(rotation=45, ha="right")
        plt.yticks(rotation=0)
        plt.tight_layout()
        path = os.path.join(figures_dir, f"global_feature_association_heatmap.{fig_format}")
        plt.savefig(path, dpi=300, bbox_inches="tight", format=fig_format)
        plt.close()
        html_path = os.path.join(html_dir, "global_feature_association_heatmap.html")
        fig = px.imshow(plot_matrix, color_continuous_scale="RdBu_r", zmin=-1, zmax=1, aspect="auto",
                        labels={"color": "Phi coefficient"}, title="Global cross-database feature association")
        fig.write_html(html_path)
        plot_paths["global_feature_association_heatmap"] = path
        plot_paths["global_feature_association_heatmap_html"] = html_path

    presence = matrix[top_features].rename(columns=labels)
    if not presence.empty:
        width = max(8, min(18, 0.35 * presence.shape[1] + 4))
        height = max(5, min(18, 0.22 * presence.shape[0] + 4))
        plt.figure(figsize=(width, height))
        sns.heatmap(presence, cmap="Blues", linewidths=0.1, linecolor="white", cbar_kws={"label": "Presence"})
        plt.title("Integrated feature presence")
        plt.xlabel("Feature")
        plt.ylabel("Sample")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        path = os.path.join(figures_dir, f"integrated_feature_presence_heatmap.{fig_format}")
        plt.savefig(path, dpi=300, bbox_inches="tight", format=fig_format)
        plt.close()
        html_path = os.path.join(html_dir, "integrated_feature_presence_heatmap.html")
        fig = px.imshow(presence, color_continuous_scale="Blues", aspect="auto",
                        labels={"color": "Presence"}, title="Integrated feature presence")
        fig.write_html(html_path)
        plot_paths["integrated_feature_presence_heatmap"] = path
        plot_paths["integrated_feature_presence_heatmap_html"] = html_path
    return plot_paths, warnings


def _write_network(pair_df, matrix, analysis_dir, figures_dir, top_n, style_config):
    edges_path = os.path.join(analysis_dir, "cross_database_feature_network_edges.csv")
    nodes_path = os.path.join(analysis_dir, "cross_database_feature_network_nodes.csv")
    html_path = os.path.join(figures_dir, "html_files", "cross_database_feature_network.html")
    warnings = []
    if pair_df.empty or matrix.empty:
        pd.DataFrame().to_csv(edges_path, index=False)
        pd.DataFrame().to_csv(nodes_path, index=False)
        return {"network_edges": edges_path, "network_nodes": nodes_path}, warnings

    edges = pair_df[(pair_df["association_scope"] == "cross_database") & (pair_df["cooccurring_samples"] > 0)].copy()
    if edges.empty:
        edges.to_csv(edges_path, index=False)
        pd.DataFrame().to_csv(nodes_path, index=False)
        return {"network_edges": edges_path, "network_nodes": nodes_path}, warnings
    edges["abs_phi"] = edges["phi_coefficient"].abs()
    edge_limit = min(max(top_n * 2, 25), style_config.get("network_max_edges", max(top_n * 2, 25)))
    warning = crowding_warning_row("cross_database_feature_network", len(edges), edge_limit, item="edges")
    if warning:
        warnings.append(warning)
    edges = edges.sort_values(["q_value", "abs_phi", "jaccard_index"], ascending=[True, False, False]).head(edge_limit)
    edges.to_csv(edges_path, index=False)

    features = sorted(set(edges["feature_1"]).union(set(edges["feature_2"])), key=lambda value: (_feature_database(value), value))
    nodes = pd.DataFrame({
        "feature": features,
        "database": [_feature_database(feature) for feature in features],
        "sample_count": [int(matrix[feature].sum()) if feature in matrix.columns else 0 for feature in features],
    })
    nodes.to_csv(nodes_path, index=False)
    warnings.extend(label_warning_rows(nodes["feature"].tolist(), style_config.get("label_max_length", 55), "cross_database_feature_network"))
    display_labels = label_map(
        nodes["feature"].tolist(),
        max_length=style_config.get("label_max_length", 55),
        wrap_width=0,
    )

    os.makedirs(os.path.dirname(html_path), exist_ok=True)
    if len(nodes) >= 2:
        coords = {}
        for index, feature in enumerate(nodes["feature"]):
            angle = (2 * math.pi * index) / len(nodes)
            coords[feature] = (math.cos(angle), math.sin(angle))
        edge_x, edge_y = [], []
        for _, row in edges.iterrows():
            x0, y0 = coords[row["feature_1"]]
            x1, y1 = coords[row["feature_2"]]
            edge_x += [x0, x1, None]
            edge_y += [y0, y1, None]
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines", line={"width": 1, "color": "#9AA6B2"}, hoverinfo="skip"))
        for database, sub_df in nodes.groupby("database"):
            xs = [coords[feature][0] for feature in sub_df["feature"]]
            ys = [coords[feature][1] for feature in sub_df["feature"]]
            fig.add_trace(go.Scatter(
                x=xs,
                y=ys,
                mode="markers+text",
                text=[display_labels.get(feature, feature) for feature in sub_df["feature"]],
                textposition="top center",
                marker={"size": 10 + sub_df["sample_count"].clip(upper=30), "color": DATABASE_COLORS.get(database, "#6B7280")},
                name=database,
                hovertext=sub_df["feature"] + "<br>Samples: " + sub_df["sample_count"].astype(str),
                hoverinfo="text",
            ))
        fig.update_layout(title="Cross-database feature association network", showlegend=True, xaxis={"visible": False}, yaxis={"visible": False})
        fig.write_html(html_path)
    return {"network_edges": edges_path, "network_nodes": nodes_path, "cross_database_feature_network_html": html_path}, warnings


def _write_metadata_enrichment(matrix, metadata, analysis_dir, min_group_size):
    combined_path = os.path.join(analysis_dir, "cross_database_feature_enrichment_by_metadata.csv")
    if matrix.empty or metadata is None or metadata.empty or "Assembly Accession" not in metadata.columns:
        pd.DataFrame().to_csv(combined_path, index=False)
        return {"feature_enrichment_by_metadata": combined_path}

    work = metadata.drop_duplicates(subset=["Assembly Accession"]).set_index("Assembly Accession")
    work = work.reindex(matrix.index)
    rows = []
    for metadata_col in [col for col in METADATA_COLUMNS if col in work.columns]:
        values = work[metadata_col].fillna("Unknown").astype(str).replace({"": "Unknown", "nan": "Unknown", "None": "Unknown"})
        groups = [group for group, sub in values.groupby(values) if group != "Unknown" and len(sub) >= min_group_size]
        if not groups:
            continue
        for group in groups:
            in_group = values == group
            if int((~in_group).sum()) < min_group_size:
                continue
            for feature in matrix.columns:
                present = matrix[feature].astype(int) == 1
                n11 = int((in_group & present).sum())
                n01 = int((in_group & ~present).sum())
                n10 = int((~in_group & present).sum())
                n00 = int((~in_group & ~present).sum())
                if n11 == 0 and n10 == 0:
                    continue
                try:
                    odds_ratio, p_value = fisher_exact([[n11, n01], [n10, n00]])
                except Exception:
                    odds_ratio, p_value = float("nan"), float("nan")
                rows.append({
                    "feature": feature,
                    "database": _feature_database(feature),
                    "metadata_variable": metadata_col,
                    "group": group,
                    "feature_present_in_group": n11,
                    "feature_absent_in_group": n01,
                    "feature_present_outside_group": n10,
                    "feature_absent_outside_group": n00,
                    "group_size": int(in_group.sum()),
                    "outside_group_size": int((~in_group).sum()),
                    "odds_ratio": odds_ratio,
                    "p_value": p_value,
                })
    enrichment = pd.DataFrame(rows)
    if not enrichment.empty:
        enrichment["q_value"] = _fdr_bh(enrichment["p_value"].tolist())
        enrichment["direction"] = "not_significant"
        enrichment.loc[(enrichment["q_value"] < 0.05) & (enrichment["odds_ratio"] > 1), "direction"] = "enriched"
        enrichment.loc[(enrichment["q_value"] < 0.05) & (enrichment["odds_ratio"] < 1), "direction"] = "depleted"
        enrichment = enrichment.sort_values(["q_value", "odds_ratio", "feature"], ascending=[True, False, True])
    enrichment.to_csv(combined_path, index=False)

    paths = {"feature_enrichment_by_metadata": combined_path}
    if not enrichment.empty:
        for metadata_col, sub_df in enrichment.groupby("metadata_variable"):
            safe_name = metadata_col.lower().replace(" ", "_").replace("/", "_")
            path = os.path.join(analysis_dir, f"feature_enrichment_by_{safe_name}.csv")
            sub_df.to_csv(path, index=False)
            paths[f"feature_enrichment_by_{safe_name}"] = path
    return paths


def _write_figure_manifest(figures_dir, output_paths):
    rows = []
    descriptions = {
        "global_feature_association_heatmap": "Phi-correlation heatmap across selected AMR, virulence, plasmid, and mobileome features.",
        "integrated_feature_presence_heatmap": "Sample-by-feature binary presence heatmap for selected integrated features.",
        "cross_database_feature_network_html": "Interactive network of strongest cross-database feature associations.",
        "AMR_vs_mobileome_burden_plot": "Burden-level scatter comparing AMR count and total mobileome count.",
        "AMR_vs_plasmid_burden_plot": "Burden-level scatter comparing AMR count and plasmid replicon count.",
        "AMR_vs_virulence_burden_plot": "Burden-level scatter comparing AMR count and virulence feature count.",
        "AMR_vs_defense_burden_plot": "Burden-level scatter comparing AMR count and defense-system count.",
        "AMR_vs_prophage_burden_plot": "Burden-level scatter comparing AMR count and prophage or viral-region count.",
    }
    for key, path in sorted(output_paths.items()):
        if path and isinstance(path, str) and os.path.exists(path) and (path.endswith((".png", ".pdf", ".svg", ".tiff", ".html"))):
            rows.append({
                "figure_file": os.path.relpath(path, figures_dir),
                "figure_key": key,
                "figure_type": "interactive_html" if path.endswith(".html") else "static",
                "database": "cross_database",
                "description": descriptions.get(key, key.replace("_", " ")),
                "recommended_use": "dashboard" if path.endswith(".html") else "publication_or_review",
                "publication_ready": bool(not path.endswith(".html")),
            })
    path = os.path.join(figures_dir, "figure_manifest.csv")
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def generate_cross_database_associations(
    output_dir,
    base_name,
    amr_tidy_df,
    feature_outputs,
    fig_format="png",
    top_n=25,
    min_prevalence=0.0,
    max_features=300,
    min_group_size=5,
    plot_style="publication",
    label_max_length=None,
):
    """Generate cross-database comparative genomics outputs."""
    cross_dir = os.path.join(output_dir, "cross_database")
    analysis_dir = os.path.join(cross_dir, "analysis")
    figures_dir = os.path.join(cross_dir, "figures")
    os.makedirs(analysis_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(os.path.join(figures_dir, "html_files"), exist_ok=True)
    style_config = plot_style_config(plot_style, label_max_length=label_max_length)
    plot_warnings = []

    matrix, feature_map, metadata = build_unified_feature_matrix(amr_tidy_df, feature_outputs)
    matrix_path = os.path.join(analysis_dir, "cross_database_feature_matrix.csv")
    matrix.to_csv(matrix_path)

    prevalence_path = os.path.join(analysis_dir, "cross_database_feature_prevalence.csv")
    if matrix.empty:
        prevalence = pd.DataFrame(columns=["feature", "database", "present_samples", "prevalence_percentage"])
    else:
        prevalence = pd.DataFrame({
            "feature": matrix.columns,
            "database": [_feature_database(feature) for feature in matrix.columns],
            "present_samples": matrix.sum(axis=0).astype(int).values,
            "prevalence_percentage": (matrix.sum(axis=0).astype(float).values / max(matrix.shape[0], 1)) * 100,
        }).sort_values(["prevalence_percentage", "present_samples", "feature"], ascending=[False, False, True])
    prevalence.to_csv(prevalence_path, index=False)

    pairwise_matrix, selected_features = _select_features_for_pairwise(matrix, max_features)
    if min_prevalence > 0 and not pairwise_matrix.empty:
        prevalence_pct = (pairwise_matrix.sum(axis=0) / max(pairwise_matrix.shape[0], 1)) * 100
        kept = prevalence_pct[prevalence_pct >= min_prevalence].index.tolist()
        pairwise_matrix = pairwise_matrix[kept]
    cooc_path, jaccard_path, phi_path = _write_pair_matrices(pairwise_matrix, analysis_dir)
    pair_df = _association_rows(pairwise_matrix, {feature: feature_map.get(feature, _feature_database(feature)) for feature in pairwise_matrix.columns})
    top_assoc_path = os.path.join(analysis_dir, "cross_database_top_associations.csv")
    pair_df.to_csv(top_assoc_path, index=False)

    outputs = {
        "cross_database_feature_matrix": matrix_path,
        "cross_database_feature_prevalence": prevalence_path,
        "cross_database_cooccurrence_matrix": cooc_path,
        "cross_database_jaccard_matrix": jaccard_path,
        "cross_database_phi_correlation_matrix": phi_path,
        "cross_database_top_associations": top_assoc_path,
        "pairwise_feature_count": len(selected_features),
        "base_name": base_name,
    }
    outputs.update(_write_specific_associations(pair_df, analysis_dir))
    outputs.update(_write_integrated_burden(matrix, metadata, analysis_dir, figures_dir, fig_format))
    heatmap_outputs, heatmap_warnings = _write_heatmaps(pairwise_matrix, phi_path, figures_dir, fig_format, top_n, style_config)
    outputs.update(heatmap_outputs)
    plot_warnings.extend(heatmap_warnings)
    network_outputs, network_warnings = _write_network(pair_df, pairwise_matrix, analysis_dir, figures_dir, top_n, style_config)
    outputs.update(network_outputs)
    plot_warnings.extend(network_warnings)
    outputs.update(_write_metadata_enrichment(pairwise_matrix, metadata, analysis_dir, min_group_size))
    warning_path = os.path.join(figures_dir, "plot_readability_warnings.csv")
    pd.DataFrame(plot_warnings, columns=["figure", "warning_type", "detail", "value", "limit"]).to_csv(warning_path, index=False)
    outputs["plot_readability_warnings"] = warning_path
    outputs["figure_manifest"] = _write_figure_manifest(figures_dir, outputs)
    return outputs
