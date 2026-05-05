import math
import os

import numpy as np
import pandas as pd
import plotly.express as px
from scipy.stats import linregress, norm, spearmanr

from panr2.associations import build_unified_feature_matrix


def _fdr_bh(p_values):
    indexed = [(i, p) for i, p in enumerate(p_values) if pd.notna(p)]
    q_values = [np.nan] * len(p_values)
    if not indexed:
        return q_values
    indexed.sort(key=lambda item: item[1])
    m = len(indexed)
    running = 1.0
    for rank, (original_index, p_value) in reversed(list(enumerate(indexed, start=1))):
        running = min(running, (float(p_value) * m) / rank)
        q_values[original_index] = min(running, 1.0)
    return q_values


def _extract_years(metadata):
    if metadata is None or metadata.empty:
        return pd.Series(dtype="float64")
    if "Collection Year" in metadata.columns:
        years = pd.to_numeric(metadata["Collection Year"], errors="coerce")
    elif "Collection Date" in metadata.columns:
        years = pd.to_numeric(metadata["Collection Date"].astype(str).str.extract(r"(\d{4})", expand=False), errors="coerce")
    else:
        years = pd.Series(np.nan, index=metadata.index)
    return years


def _mann_kendall(values):
    values = [float(value) for value in values if pd.notna(value)]
    n = len(values)
    if n < 3:
        return 0.0, 0.0, 1.0, "insufficient_years"
    s_value = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            s_value += np.sign(values[j] - values[i])
    var_s = n * (n - 1) * (2 * n + 5) / 18
    if var_s == 0:
        return float(s_value), 0.0, 1.0, "no_trend"
    if s_value > 0:
        z_value = (s_value - 1) / math.sqrt(var_s)
    elif s_value < 0:
        z_value = (s_value + 1) / math.sqrt(var_s)
    else:
        z_value = 0.0
    p_value = 2 * norm.sf(abs(z_value))
    trend = "increasing" if z_value > 0 and p_value < 0.05 else "decreasing" if z_value < 0 and p_value < 0.05 else "not_significant"
    return float(s_value), float(z_value), float(p_value), trend


def _logistic_trend(years, y_values):
    years = np.asarray(years, dtype=float)
    y_values = np.asarray(y_values, dtype=float)
    if len(np.unique(y_values)) < 2 or len(years) < 4:
        return np.nan, np.nan, np.nan
    x = years - np.nanmean(years)
    design = np.column_stack([np.ones(len(x)), x])
    beta = np.zeros(2)
    for _ in range(50):
        eta = design @ beta
        prob = 1 / (1 + np.exp(-np.clip(eta, -30, 30)))
        weights = np.clip(prob * (1 - prob), 1e-6, None)
        z = eta + (y_values - prob) / weights
        xtw = design.T * weights
        try:
            beta_new = np.linalg.solve(xtw @ design, xtw @ z)
        except np.linalg.LinAlgError:
            return np.nan, np.nan, np.nan
        if np.max(np.abs(beta_new - beta)) < 1e-7:
            beta = beta_new
            break
        beta = beta_new
    try:
        cov = np.linalg.inv((design.T * weights) @ design)
        se = math.sqrt(cov[1, 1])
        z_score = beta[1] / se if se else np.nan
        p_value = 2 * norm.sf(abs(z_score)) if pd.notna(z_score) else np.nan
    except Exception:
        p_value = np.nan
    return float(beta[1]), float(np.exp(beta[1])), float(p_value)


def write_temporal_trends(output_dir, base_name, amr_tidy_df, feature_outputs=None, fig_format="png", min_years=3):
    """Write advanced temporal feature and burden trend summaries."""
    temporal_dir = os.path.join(output_dir, "temporal")
    analysis_dir = os.path.join(temporal_dir, "analysis")
    figures_dir = os.path.join(temporal_dir, "figures")
    html_dir = os.path.join(figures_dir, "html_files")
    for path in [analysis_dir, figures_dir, html_dir]:
        os.makedirs(path, exist_ok=True)

    feature_path = os.path.join(analysis_dir, "temporal_feature_trends.csv")
    burden_path = os.path.join(analysis_dir, "temporal_burden_trends.csv")
    matrix, feature_map, metadata = build_unified_feature_matrix(amr_tidy_df, feature_outputs or {})
    if matrix.empty or metadata.empty:
        pd.DataFrame().to_csv(feature_path, index=False)
        pd.DataFrame().to_csv(burden_path, index=False)
        return {"temporal_feature_trends": feature_path, "temporal_burden_trends": burden_path}

    metadata = metadata.drop_duplicates(subset=["Assembly Accession"]).set_index("Assembly Accession").reindex(matrix.index)
    years = _extract_years(metadata)
    valid = years.notna()
    matrix = matrix.loc[valid]
    years = years.loc[valid].astype(int)
    if matrix.empty or years.nunique() < min_years:
        pd.DataFrame().to_csv(feature_path, index=False)
        pd.DataFrame().to_csv(burden_path, index=False)
        return {"temporal_feature_trends": feature_path, "temporal_burden_trends": burden_path}

    feature_rows = []
    yearly = matrix.copy()
    yearly["collection_year"] = years.values
    for feature in matrix.columns:
        prevalence = yearly.groupby("collection_year")[feature].mean().reset_index()
        if len(prevalence) < min_years:
            continue
        s_value, mk_z, mk_p, mk_trend = _mann_kendall(prevalence[feature].tolist())
        try:
            spearman_r, spearman_p = spearmanr(prevalence["collection_year"], prevalence[feature])
        except Exception:
            spearman_r, spearman_p = np.nan, np.nan
        slope, odds_ratio, logistic_p = _logistic_trend(years.values, matrix[feature].values)
        feature_rows.append({
            "feature": feature,
            "database": feature_map.get(feature, feature.split(":", 1)[0]),
            "years_observed": int(len(prevalence)),
            "first_year": int(prevalence["collection_year"].min()),
            "last_year": int(prevalence["collection_year"].max()),
            "mean_prevalence": float(prevalence[feature].mean()),
            "mann_kendall_s": s_value,
            "mann_kendall_z": mk_z,
            "mann_kendall_p_value": mk_p,
            "mann_kendall_trend": mk_trend,
            "spearman_r": spearman_r,
            "spearman_p_value": spearman_p,
            "logistic_slope_per_year": slope,
            "logistic_odds_ratio_per_year": odds_ratio,
            "logistic_p_value": logistic_p,
        })
    feature_df = pd.DataFrame(feature_rows)
    if not feature_df.empty:
        feature_df["mann_kendall_q_value"] = _fdr_bh(feature_df["mann_kendall_p_value"].tolist())
        feature_df["spearman_q_value"] = _fdr_bh(feature_df["spearman_p_value"].tolist())
        feature_df["logistic_q_value"] = _fdr_bh(feature_df["logistic_p_value"].tolist())
        feature_df = feature_df.sort_values(["mann_kendall_p_value", "mean_prevalence", "feature"], ascending=[True, False, True])
    feature_df.to_csv(feature_path, index=False)

    burden = pd.DataFrame({"collection_year": years.values})
    for database in sorted(set(feature_map.values())):
        cols = [feature for feature, db in feature_map.items() if db == database and feature in matrix.columns]
        if cols:
            burden[f"{database.lower()}_feature_count"] = matrix[cols].sum(axis=1).values
    burden_rows = []
    for count_col in [col for col in burden.columns if col.endswith("_feature_count")]:
        per_sample = burden[["collection_year", count_col]].dropna()
        if per_sample["collection_year"].nunique() < min_years:
            continue
        try:
            lr = linregress(per_sample["collection_year"], per_sample[count_col])
            slope = lr.slope
            p_value = lr.pvalue
            r_value = lr.rvalue
        except Exception:
            slope, p_value, r_value = np.nan, np.nan, np.nan
        yearly_mean = per_sample.groupby("collection_year")[count_col].mean().tolist()
        s_value, mk_z, mk_p, mk_trend = _mann_kendall(yearly_mean)
        burden_rows.append({
            "burden_metric": count_col,
            "years_observed": int(per_sample["collection_year"].nunique()),
            "linear_slope_per_year": slope,
            "linear_r_value": r_value,
            "linear_p_value": p_value,
            "mann_kendall_s": s_value,
            "mann_kendall_z": mk_z,
            "mann_kendall_p_value": mk_p,
            "mann_kendall_trend": mk_trend,
        })
    burden_df = pd.DataFrame(burden_rows)
    if not burden_df.empty:
        burden_df["linear_burden_q_value"] = _fdr_bh(burden_df["linear_p_value"].tolist())
        burden_df["mann_kendall_q_value"] = _fdr_bh(burden_df["mann_kendall_p_value"].tolist())
    burden_df.to_csv(burden_path, index=False)

    if not feature_df.empty:
        top_features = feature_df.head(20)["feature"].tolist()
        trend_df = yearly.groupby("collection_year")[top_features].mean().reset_index().melt(
            id_vars="collection_year",
            var_name="feature",
            value_name="prevalence",
        )
        fig = px.line(trend_df, x="collection_year", y="prevalence", color="feature",
                      markers=True, title="Top feature prevalence trends")
        html_path = os.path.join(html_dir, "temporal_top_feature_trends.html")
        fig.write_html(html_path)
        return {
            "temporal_feature_trends": feature_path,
            "temporal_burden_trends": burden_path,
            "temporal_top_feature_trends_html": html_path,
        }
    return {"temporal_feature_trends": feature_path, "temporal_burden_trends": burden_path}
