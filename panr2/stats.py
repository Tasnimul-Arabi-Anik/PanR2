import logging
import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.stats import pearsonr, spearmanr


def correlation_scatterplot_analysis(tidy_file, output_dir, group_col, min_samples_per_group=5):
    """
    Generate scatterplots and correlation analysis for NUM_FOUND vs Collection Date (Year),
    grouped by the specified column (e.g., 'Continent', 'Geographic Location').
    Uses NumPy polyfit for regression line. Saves interactive Plotly HTML and correlation summary CSV.
    """
    try:
        df = pd.read_csv(tidy_file)
        df = df[["Assembly BioSample Accession", "Collection Date", group_col, "NUM_FOUND"]].dropna().drop_duplicates()

        # Clean collection date
        df = df[~df["Collection Date"].isin(["absent", "none", "", None])]

        def parse_year(x):
            try:
                return int(str(x)[:4])
            except Exception:
                return np.nan

        df["Year"] = df["Collection Date"].apply(parse_year)
        df = df.dropna(subset=["Year"])
        df["Year"] = df["Year"].astype(int)

        unique_values = []
        valid_groups = {}

        for val in sorted(df[group_col].dropna().unique()):
            sub_df = df[df[group_col] == val]
            if len(sub_df) >= min_samples_per_group:
                unique_values.append(val)
                valid_groups[val] = sub_df

        if not unique_values:
            corr_csv = os.path.join(output_dir, f"{group_col.replace(' ', '_')}_correlation_summary.csv")
            pd.DataFrame(columns=[
                group_col, "n_samples", "pearson_r", "pearson_p", "spearman_r", "spearman_p"
            ]).to_csv(corr_csv, index=False)
            logging.warning(
                f"No {group_col} groups have sufficient data for correlation analysis "
                f"(minimum {min_samples_per_group} samples required). Skipping."
            )
            return None

        records = []
        fig = go.Figure()

        for i, val in enumerate(unique_values):
            sub_df = valid_groups[val]
            x = sub_df["Year"].values
            y = sub_df["NUM_FOUND"].values

            pearson_corr, pearson_p = pearsonr(x, y)
            spearman_corr, spearman_p = spearmanr(x, y)

            records.append({
                group_col: val,
                "n_samples": len(sub_df),
                "pearson_r": pearson_corr,
                "pearson_p": pearson_p,
                "spearman_r": spearman_corr,
                "spearman_p": spearman_p
            })

            fig.add_trace(go.Scatter(
                x=x,
                y=y,
                mode="markers",
                name=val,
                visible=(i == 0),
                marker=dict(size=8, opacity=0.6),
                text=sub_df["Assembly BioSample Accession"],
                hovertemplate=(
                    f"<b>{group_col}</b>: {val}<br>"
                    "Year: %{x}<br>"
                    "NUM_FOUND: %{y}<br>"
                    "Sample: %{text}<extra></extra>"
                )
            ))

            # Add regression line using np.polyfit
            slope, intercept = np.polyfit(x, y, 1)
            y_pred = slope * x + intercept
            fig.add_trace(go.Scatter(
                x=x,
                y=y_pred,
                mode="lines",
                name=f"{val} Trend",
                line=dict(dash="dash", color="black"),
                visible=(i == 0),
                showlegend=False
            ))

        corr_df = pd.DataFrame(records)
        corr_csv = os.path.join(output_dir, f"{group_col.replace(' ', '_')}_correlation_summary.csv")
        corr_df.to_csv(corr_csv, index=False)
        logging.info(f"Saved correlation summary to {corr_csv}")

        dropdown_buttons = []
        for i, val in enumerate(unique_values):
            visibility = [False] * len(unique_values) * 2  # 2 traces per group (scatter + line)
            visibility[i * 2] = True       # scatter
            visibility[i * 2 + 1] = True   # line

            stats = records[i]
            dropdown_buttons.append(dict(
                label=f"{val} (n={stats['n_samples']})",
                method="update",
                args=[
                    {"visible": visibility},
                    {
                        "title": f"NUM_FOUND vs Year - {group_col}: {val}<br>" +
                                 f"<sub>Pearson r={stats['pearson_r']:.3f} (p={stats['pearson_p']:.3f}), " +
                                 f"Spearman ρ={stats['spearman_r']:.3f} (p={stats['spearman_p']:.3f})</sub>"
                    }
                ]
            ))

        initial_stats = records[0]
        initial_title = (f"NUM_FOUND vs Year - {group_col}: {unique_values[0]}<br>" +
                         f"<sub>Pearson r={initial_stats['pearson_r']:.3f} (p={initial_stats['pearson_p']:.3f}), " +
                         f"Spearman ρ={initial_stats['spearman_r']:.3f} (p={initial_stats['spearman_p']:.3f})</sub>")

        fig.update_layout(
            title=initial_title,
            xaxis_title="Collection Year",
            yaxis_title="Number of Resistance Genes (NUM_FOUND)",
            updatemenus=[dict(
                buttons=dropdown_buttons,
                direction="down",
                showactive=True,
                x=0.02,
                y=0.98,
                xanchor="left",
                yanchor="top",
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="rgba(0,0,0,0.2)",
                borderwidth=1
            )],
            margin=dict(t=120, l=60, r=60, b=60),
            template="plotly_white"
        )

        fig.add_annotation(
            text=f"Select {group_col}:",
            x=0.02,
            y=1.02,
            xref="paper",
            yref="paper",
            showarrow=False,
            font=dict(size=12, color="black")
        )

        html_out = os.path.join(output_dir, f"{group_col.replace(' ', '_')}_correlation_plot.html")
        fig.write_html(html_out)
        logging.info(f"Saved interactive plot to {html_out}")

        print(f"\nCorrelation Analysis Summary for {group_col}:")
        print(f"Total groups processed: {len(unique_values)}")
        for record in records:
            print(f"  {record[group_col]}: {record['n_samples']} samples, "
                  f"Pearson r={record['pearson_r']:.3f}")

        return fig

    except Exception as e:
        logging.error(f"Error in correlation analysis: {e}")
        raise

def combined_correlation_analysis(output_dir):
    """
    Generate combined correlation analysis for all grouping variables in the dataset.
    Saves individual correlation scatterplots and a summary CSV file.
    """
    try:
        # Read the three generated CSV files
        geo_location_df = pd.read_csv(os.path.join(output_dir, "Geographic_Location_correlation_summary.csv"))
        continent_df = pd.read_csv(os.path.join(output_dir, "Continent_correlation_summary.csv"))
        subcontinent_df = pd.read_csv(os.path.join(output_dir, "Subcontinent_correlation_summary.csv"))
        
        # Add Geographic_Level column to identify the type of geographic grouping
        geo_location_df['Geographic_Level'] = 'Geographic Location'
        continent_df['Geographic_Level'] = 'Continent'
        subcontinent_df['Geographic_Level'] = 'Subcontinent'
        
        # Rename the geographic columns to a consistent name
        geo_location_df = geo_location_df.rename(columns={'Geographic Location': 'Geographic_Region'})
        continent_df = continent_df.rename(columns={'Continent': 'Geographic_Region'})
        subcontinent_df = subcontinent_df.rename(columns={'Subcontinent': 'Geographic_Region'})
        
        # Combine all dataframes
        combined_df = pd.concat([geo_location_df, continent_df, subcontinent_df], ignore_index=True)
        
        # Reorder columns for better readability
        combined_df = combined_df[['Geographic_Level', 'Geographic_Region', 'n_samples', 
                                 'pearson_r', 'pearson_p', 'spearman_r', 'spearman_p']]
        
        # Save combined CSV
        combined_csv = os.path.join(output_dir, "combined_geographic_correlation_summary.csv")
        combined_df.to_csv(combined_csv, index=False)

        # Remove the three separate CSV files
        separate_files = [
            os.path.join(output_dir, "Geographic_Location_correlation_summary.csv"),
            os.path.join(output_dir, "Continent_correlation_summary.csv"),
            os.path.join(output_dir, "Subcontinent_correlation_summary.csv")
        ]
        
        print("Removing separate correlation...CSV files...")
        for file_path in separate_files:
            try:
                if os.path.exists(file_path):
                    os.remove(file_path)
                    print(f"  Removed: {os.path.basename(file_path)}")
                else:
                    print(f"  File not found: {os.path.basename(file_path)}")
            except Exception as e:
                print(f"  Error removing {os.path.basename(file_path)}: {e}")
        
        print("Cleanup completed. Only combined correlation...CSV file remains.")
    except Exception as e:
        print(f"Error in combined_correlation_analysis: {e}")

