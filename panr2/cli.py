import argparse
import glob
import logging
import os
import shutil

import pandas as pd

from panr2.analysis import generate_comprehensive_analysis_outputs
from panr2.features import analyze_abricate_feature_database
from panr2.filters import apply_analysis_filters
from panr2.io import (
    convert_tab_to_csv,
    convert_to_tidy_format,
    load_and_merge_data,
    save_merged_data,
    unique_input_files,
)
from panr2.plots import (
    analyze_gene_presence,
    generate_comparison_heatmap,
    generate_comparison_heatmap_plotly,
    generate_gene_identity_boxplot,
    generate_gene_identity_boxplot_plotly,
    generate_geographic_resistance_map_plotly,
    generate_index_html,
    generate_mean_arg_lollipop,
    generate_mean_arg_lollipop_plotly,
    generate_resistance_barplot,
    mean_Arg_resistance_analysis_plotly,
)
from panr2.qc import write_input_qc_report
from panr2.report import write_report
from panr2.stats import combined_correlation_analysis, correlation_scatterplot_analysis


PANR2_VERSION = "0.1.3-dev"


# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def main(ncbi_dir, abricate_dir, output_dir, fig_format, nseq, genep, min_identity=0.0, drop_unmatched_accessions=False, min_samples_per_group=5, core_threshold=95.0, rare_threshold=5.0, top_n=25, cooccurrence_min_prevalence=0.0, cooccurrence_top_n=25, vfdb_dir=None, plasmidfinder_dir=None):
    """Main function to process data and generate outputs."""
    logging.info("Starting the script.")

    if min_identity < 0 or min_identity > 100:
        raise ValueError("--min-identity must be between 0 and 100.")
    if min_samples_per_group < 1:
        raise ValueError("--min-samples-per-group must be at least 1.")
    if core_threshold < 0 or core_threshold > 100:
        raise ValueError("--core-threshold must be between 0 and 100.")
    if rare_threshold < 0 or rare_threshold > 100:
        raise ValueError("--rare-threshold must be between 0 and 100.")
    if rare_threshold > core_threshold:
        raise ValueError("--rare-threshold cannot be greater than --core-threshold.")
    if top_n < 1:
        raise ValueError("--top-n must be at least 1.")
    if cooccurrence_min_prevalence < 0 or cooccurrence_min_prevalence > 100:
        raise ValueError("--cooccurrence-min-prevalence must be between 0 and 100.")
    if cooccurrence_top_n < 1:
        raise ValueError("--cooccurrence-top-n must be at least 1.")
    
    # Define paths
    ncbi_clean_path = os.path.join(ncbi_dir, "ncbi_clean.csv")
    abricate_summary_files = unique_input_files(
        glob.glob(os.path.join(abricate_dir, "*summary.[ct][sa][bv]"))
    )  # Match .csv and .tab files
    abricate_results_files = unique_input_files(
        glob.glob(os.path.join(abricate_dir, "*results.[ct][sa][bv]"))
    )
    
    if not os.path.exists(ncbi_clean_path):
        raise FileNotFoundError(f"ncbi_clean.csv not found in {ncbi_dir}.")
    
    if not abricate_summary_files:
        raise FileNotFoundError(f"No CSV or TAB summary files (abricate) found in {abricate_dir}.")
        
    if not abricate_results_files:
        raise FileNotFoundError(f"No CSV or TAB results files (abricate) found in {abricate_dir}.")
    
    # Create output subdirectories
    merged_output_dir = os.path.join(output_dir, "merged_output")
    figures_dir = os.path.join(output_dir, "figures")
    os.makedirs(merged_output_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    first_summary_file = sorted(abricate_summary_files)[0]
    first_results_file = sorted(abricate_results_files)[0]
    write_input_qc_report(ncbi_clean_path, first_summary_file, first_results_file, output_dir)

    optional_feature_outputs = {}
    if vfdb_dir:
        optional_feature_outputs["vfdb"] = analyze_abricate_feature_database(
            ncbi_clean_path, vfdb_dir, output_dir, "vfdb", "virulence", min_identity=min_identity, fig_format=fig_format
        )
    if plasmidfinder_dir:
        optional_feature_outputs["plasmidfinder"] = analyze_abricate_feature_database(
            ncbi_clean_path, plasmidfinder_dir, output_dir, "plasmidfinder", "plasmid", min_identity=min_identity, fig_format=fig_format
        )
    
    for abricate_summary_file in abricate_summary_files:
        try:
            # Get base name without extension
            base = os.path.splitext(abricate_summary_file)[0]
            csv_file = base + ".csv"
            tab_file = base + ".tab"

            if os.path.exists(csv_file):
                abricate_summary_file = csv_file
            elif os.path.exists(tab_file):
                convert_tab_to_csv(tab_file, csv_file)
                abricate_summary_file = csv_file
            else:
                print(f"⚠️ No .csv or .tab found for {os.path.basename(base)}, skipping.")
                continue

            # Load, filter, and merge data
            merged_df = load_and_merge_data(ncbi_clean_path, abricate_summary_file)
            output_filename = f"ncbi_{os.path.basename(abricate_summary_file)}"
            base_name = os.path.basename(abricate_summary_file).replace("_summary.csv", "")
            merged_df = apply_analysis_filters(
                merged_df,
                min_identity=min_identity,
                drop_unmatched_accessions=drop_unmatched_accessions,
                output_dir=output_dir,
                base_name=base_name,
            )
            
            # Save merged data
            save_merged_data(merged_df, merged_output_dir, output_filename)
            
            # Convert to tidy format
            tidy_df = convert_to_tidy_format(merged_df)

            # Add RESISTANCE column from corresponding results file
            basename = os.path.basename(abricate_summary_file).replace("_summary.csv", "")
            expected_results_file = os.path.join(abricate_dir, f"{basename}_results.csv")

            # Convert .tab to .csv if needed
            if not os.path.exists(expected_results_file):
                tab_version = expected_results_file.replace(".csv", ".tab")
                if os.path.exists(tab_version):
                    convert_tab_to_csv(tab_version, expected_results_file)
                    logging.info(f"Converted {tab_version} to CSV.")
                else:
                    logging.warning(f"Results file not found for {basename}. Skipping RESISTANCE enrichment.")
                    expected_results_file = None
            # Merge RESISTANCE and GENE data
            if expected_results_file and os.path.exists(expected_results_file):
                try:
                    results_df = pd.read_csv(expected_results_file, dtype=str)
                    results_df = results_df[['GENE', 'RESISTANCE']].drop_duplicates()

                    # Perform the merge using 'Gene' from tidy_df and 'GENE' from results_df
                    tidy_df = tidy_df.merge(results_df, left_on='Gene', right_on='GENE', how='left')

                    # Drop redundant GENE column
                    tidy_df.drop(columns=['GENE'], inplace=True)

                    logging.info(f"Successfully added RESISTANCE info for {basename}.")
                except Exception as e:
                    logging.error(f"Error merging RESISTANCE data from {expected_results_file}: {e}")

            tidy_file = os.path.join(merged_output_dir, output_filename.replace("_summary.csv", "_tidy_summary.csv"))
            tidy_df.to_csv(tidy_file, index=False)
            logging.info(f"Tidied file saved to {tidy_file}")

            generate_comprehensive_analysis_outputs(
                tidy_df,
                output_dir,
                base_name,
                fig_format,
                core_threshold=core_threshold,
                rare_threshold=rare_threshold,
                top_n=top_n,
                cooccurrence_min_prevalence=cooccurrence_min_prevalence,
                cooccurrence_top_n=cooccurrence_top_n,
            )
            
            # Analyze gene prevalence and generate figures
            analyze_gene_presence(tidy_df, figures_dir, base_name, fig_format)
            

            # Generate boxplot for gene identity
            generate_gene_identity_boxplot(tidy_file, figures_dir, fig_format)
            # Generate interactive boxplot for gene identity
            generate_gene_identity_boxplot_plotly(tidy_file, figures_dir)

            
            # Create subdirectory 'mean_ARG' inside figures_dir
            mean_ARG_dir = os.path.join(figures_dir, "mean_ARG")
            os.makedirs(mean_ARG_dir, exist_ok=True)
                # Generate mean ARG based on different groups
            generate_mean_arg_lollipop(tidy_file, mean_ARG_dir, fig_format, group_by="Geographic Location")
            generate_mean_arg_lollipop(tidy_file, mean_ARG_dir, fig_format, group_by="Collection Date")
            generate_mean_arg_lollipop(tidy_file, mean_ARG_dir, fig_format, group_by="Continent")
            generate_mean_arg_lollipop(tidy_file, mean_ARG_dir, fig_format, group_by="Subcontinent")


            # Generate barplots for resistance compariosn
            generate_resistance_barplot(tidy_file, figures_dir, fig_format)

            # Create subdirectory 'heatmap' inside figures_dir
            heatmap_dir = os.path.join(figures_dir, "heatmap")
            os.makedirs(heatmap_dir, exist_ok=True)
                # Generate country comparison heatmap
            generate_comparison_heatmap(tidy_file, heatmap_dir, fig_format, group_col="Geographic Location", resistance_col="RESISTANCE", genep_threshold=genep, nseq_threshold=nseq)
            generate_comparison_heatmap(tidy_file, heatmap_dir, fig_format, group_col="Collection Date", resistance_col="RESISTANCE", genep_threshold=genep, nseq_threshold=nseq)
            generate_comparison_heatmap(tidy_file, heatmap_dir, fig_format, group_col="Continent", resistance_col="RESISTANCE", genep_threshold=genep, nseq_threshold=nseq)
            generate_comparison_heatmap(tidy_file, heatmap_dir, fig_format, group_col="Subcontinent", resistance_col="RESISTANCE", genep_threshold=genep, nseq_threshold=nseq)
            

            # Plotly plots
            # Generate interactive heatmap with multiple grouping options
            generate_comparison_heatmap_plotly(tidy_file, figures_dir)
            # Generate geographic resistance map
            generate_geographic_resistance_map_plotly(tidy_file, figures_dir)
            # Generate mean ARG resistance analysis box plot
            mean_Arg_resistance_analysis_plotly(tidy_file, figures_dir)
            # Generate interactive lollipop plot with dropdown for group selection
            generate_mean_arg_lollipop_plotly(tidy_file, figures_dir)


            # Generate correlation scatterplot analysis
            print("Generating Geographic Location analysis...")
            correlation_scatterplot_analysis(tidy_file, figures_dir, group_col="Geographic Location", min_samples_per_group=min_samples_per_group)
            print("Generating Continent analysis...")
            correlation_scatterplot_analysis(tidy_file, figures_dir, group_col="Continent", min_samples_per_group=min_samples_per_group)
            print("Generating Subcontinent analysis...")
            correlation_scatterplot_analysis(tidy_file, figures_dir, group_col="Subcontinent", min_samples_per_group=min_samples_per_group)
             # Combine the three CSV files
            print("Combining correlation summary CSV files...")
            combined_correlation_analysis(figures_dir)


            # Keeping all the html files in html_dir inside figures_dir
            html_dir = os.path.join(figures_dir, "html_files")
            os.makedirs(html_dir, exist_ok=True)
            # Move all HTML files to html_dir
            for file in os.listdir(figures_dir):
                if file.endswith(".html"):
                    src_path = os.path.join(figures_dir, file)
                    dest_path = os.path.join(html_dir, file)
                    os.rename(src_path, dest_path)
                    logging.info(f"Moved {file} to {html_dir}")
            # Generate index.html for easy navigation
            generate_index_html(html_dir, figures_dir)

            # Moves every .csv file in figures_dir into that subdirectory.
            stat_analysis_dir = os.path.join(figures_dir, "Stat_analysis")
            os.makedirs(stat_analysis_dir, exist_ok=True)

            # Move all .csv files from figures_dir to figures_dir/Stat_analysis
            for file in os.listdir(figures_dir):
                if file.endswith(".csv"):
                    src = os.path.join(figures_dir, file)
                    dst = os.path.join(stat_analysis_dir, file)
                    shutil.move(src, dst)

            write_report(
                output_dir,
                base_name,
                options={
                    "fig_format": fig_format,
                    "nseq": nseq,
                    "genep": genep,
                    "min_identity": min_identity,
                    "drop_unmatched_accessions": drop_unmatched_accessions,
                    "min_samples_per_group": min_samples_per_group,
                    "core_threshold": core_threshold,
                    "rare_threshold": rare_threshold,
                    "top_n": top_n,
                    "cooccurrence_min_prevalence": cooccurrence_min_prevalence,
                    "cooccurrence_top_n": cooccurrence_top_n,
                    "vfdb_dir": vfdb_dir or "not provided",
                    "plasmidfinder_dir": plasmidfinder_dir or "not provided",
                },
                panr2_version=PANR2_VERSION,
                feature_outputs=optional_feature_outputs,
                input_files={
                    "ncbi_clean": ncbi_clean_path,
                    "abricate_summary": abricate_summary_file,
                    "abricate_results": expected_results_file or "not available",
                },
            )

        except Exception as e:
            logging.error(f"Error processing {abricate_summary_file}: {e}")
    
    logging.info("panr run successfully.")


def run_cli():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process NCBI and Abricate data.")
    parser.add_argument("--ncbi-dir", required=True, help="Directory containing ncbi_clean.csv.")
    parser.add_argument("--abricate-dir", required=True, help="Directory containing Abricate summary CSV or TAB files.")
    parser.add_argument("--output-dir", required=True, help="Base output directory.")
    parser.add_argument("--genep", type=float, default=10.0, help="Minimum %% gene presence to include in heatmap.")
    parser.add_argument("--nseq", type=int, default=1, help="Minimum number of sequences required per group in heatmaps.")
    parser.add_argument("--format", default="tiff", choices=["tiff", "svg", "png", "pdf"], help="Output format for figures (tiff, svg, png, pdf).")
    parser.add_argument("--min-identity", type=float, default=0.0, help="Minimum ABRicate identity percentage to treat a gene call as present. Values below this are set to absent.")
    parser.add_argument("--drop-unmatched-accessions", action="store_true", help="Drop NCBI assemblies that do not have a matching ABRicate summary row after merge.")
    parser.add_argument("--min-samples-per-group", type=int, default=5, help="Minimum samples per group required for correlation analyses.")
    parser.add_argument("--core-threshold", type=float, default=95.0, help="Prevalence percentage used to classify core ARGs.")
    parser.add_argument("--rare-threshold", type=float, default=5.0, help="Prevalence percentage used to classify rare ARGs.")
    parser.add_argument("--top-n", type=int, default=25, help="Number of top genes/classes to include in compact summary plots.")
    parser.add_argument("--cooccurrence-min-prevalence", type=float, default=0.0, help="Minimum prevalence percentage for genes/classes included in co-occurrence matrices.")
    parser.add_argument("--cooccurrence-top-n", type=int, default=25, help="Number of top genes/classes or pairs to include in co-occurrence plots and pair tables.")
    parser.add_argument("--vfdb-dir", help="Optional directory containing ABRicate VFDB summary/results files.")
    parser.add_argument("--plasmidfinder-dir", help="Optional directory containing ABRicate PlasmidFinder summary/results files.")
    parser.add_argument('--version', action='version', version=f'PanR2 {PANR2_VERSION}')

    args = parser.parse_args()
    
    # Run the main function
    main(
        args.ncbi_dir,
        args.abricate_dir,
        args.output_dir,
        args.format,
        args.nseq,
        args.genep,
        min_identity=args.min_identity,
        drop_unmatched_accessions=args.drop_unmatched_accessions,
        min_samples_per_group=args.min_samples_per_group,
        core_threshold=args.core_threshold,
        rare_threshold=args.rare_threshold,
        top_n=args.top_n,
        cooccurrence_min_prevalence=args.cooccurrence_min_prevalence,
        cooccurrence_top_n=args.cooccurrence_top_n,
        vfdb_dir=args.vfdb_dir,
        plasmidfinder_dir=args.plasmidfinder_dir,
    )


if __name__ == "__main__":
    run_cli()
