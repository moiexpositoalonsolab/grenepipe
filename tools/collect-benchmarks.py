#!/usr/bin/env python3

"""
Combine snakemake per-sample benchmark files and create scatter plots
relating sample sizes to the benchmark columns.

This helps to fine-tune the runtime and memory requirements for each rule,
in order to scale up to very large dataset analyses.

Usage:
  python collect-benchmarks.py <analysis_dir>

for a given grenpipe analysis run.
"""

import os
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import yaml


def get_samples_table_path(config_file_path):
    """
    Load the config yaml file and retrieve the value of the 'samples-table' key
    under the 'data' entry.
    """
    with open(config_file_path, 'r', encoding='utf-8') as file:
        try:
            content = yaml.safe_load(file)
        except yaml.YAMLError as exc:
            raise yaml.YAMLError(f"Error parsing config yaml file: {exc}")

    if 'data' not in content or 'samples-table' not in content['data']:
        raise KeyError(
            "config yaml file does not contain the expected structure: 'data' -> 'samples-table'"
        )

    return content['data']['samples-table']


def parse_samples_table(samples_path):
    """
    Read the samples table and compute the total size of fq1 + fq2 (if present).
    Return two dictionaries:
      1) sample_unit_size: keyed by 'sample-unit'
      2) sample_size: keyed by 'sample'
    """
    sample_unit_size = {}
    sample_size = {}

    with open(samples_path, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        # We expect at least:
        #   sample, unit, platform, fq1, fq2
        # Let's just find their indices dynamically
        col_idx = {colname: i for i, colname in enumerate(header)}

        # We'll require at least 'sample' and 'fq1' columns
        for required_col in ["sample", "unit", "fq1", "fq2"]:
            if required_col not in col_idx:
                raise ValueError(f"Missing required column '{required_col}' in samples table!")

        # Process each row
        for line in f:
            line = line.strip()
            if not line:
                continue
            row = line.split("\t")
            sample_name = row[col_idx["sample"]]
            unit_name = row[col_idx["unit"]]
            fq1_path = row[col_idx["fq1"]]
            fq2_path = row[col_idx["fq2"]] if col_idx["fq2"] < len(row) else None

            # Compute file sizes (some might be empty)
            size_total = 0
            for fq_path in [fq1_path, fq2_path]:
                if fq_path:  # non-empty
                    if os.path.isfile(fq_path):
                        size_total += os.path.getsize(fq_path)
                    else:
                        # It's possible the path doesn't exist or is invalid.
                        # We'll just skip or treat it as size 0.
                        pass

            # Build "sample-unit" key if unit_name is present
            sample_unit_key = f"{sample_name}-{unit_name}"

            # Store in sample_unit_size
            if sample_unit_key not in sample_unit_size:
                sample_unit_size[sample_unit_key] = 0
            sample_unit_size[sample_unit_key] += size_total

            # Also store in sample_size (sum across all units)
            if sample_name not in sample_size:
                sample_size[sample_name] = 0
            sample_size[sample_name] += size_total

    return sample_unit_size, sample_size


def combine_benchmarks(top_level_dir, output_dir):
    """
    Traverse top_level_dir, and for each subdirectory (including nested),
    combine all files into a single TSV in output_dir. The name of the output
    file is derived from the subdirectory path with slashes replaced by underscores.
    The first column is 'filename', added for each data row.
    Return a list of paths to the newly created combined TSV files.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    combined_files = []  # to store the paths to the combined TSVs

    for root, dirs, files in os.walk(top_level_dir):
        print("Processing:", root)

        # If there are no files in this directory, skip it
        if not files:
            continue

        # Create a name for this subdirectory, relative to top_level_dir
        rel_dir = os.path.relpath(root, top_level_dir)

        # Convert slashes to underscores
        if rel_dir == ".":
            subdir_name = "top_level"
        else:
            subdir_name = rel_dir.replace(os.sep, "_")

        # Full path of the output file for this subdirectory
        output_file_name = subdir_name + ".tsv"
        output_file_path = os.path.join(output_dir, output_file_name)

        # We'll open the output file once and append lines as we go
        did_write_header = False
        with open(output_file_path, "w", encoding="utf-8") as outfile:
            # Process each regular file in the current directory
            for f in files:
                file_path = os.path.join(root, f)
                if not os.path.isfile(file_path):
                    continue  # skip directories, symlinks, etc.

                # Derive the base name by removing the extension
                base_name = os.path.splitext(f)[0]

                with open(file_path, "r", encoding="utf-8") as infile:
                    for i, line in enumerate(infile):
                        line = line.rstrip("\n")
                        if i == 0:
                            # Header line
                            if not did_write_header:
                                # Prepend "filename" to the header
                                final_header = "filename\t" + line
                                outfile.write(final_header + "\n")
                                did_write_header = True
                            # If we already wrote a header, ignore subsequent
                        else:
                            # Data line => prepend the filename
                            final_data = base_name + "\t" + line
                            outfile.write(final_data + "\n")

        # Check if we ended up writing anything
        if os.path.exists(output_file_path) and os.path.getsize(output_file_path) == 0:
            os.remove(output_file_path)
        else:
            print(f"Created: {output_file_path}")
            combined_files.append(output_file_path)

    return combined_files


def seconds_to_hms(value, pos):
    """Convert 'value' (in seconds) to a short human-readable form."""
    hours = int(value // 3600)
    minutes = int((value % 3600) // 60)
    seconds = int(value % 60)
    value = round(value, 1)
    if hours > 0:
        return f"{hours}h{minutes}m"
    elif minutes > 0:
        return f"{minutes}m{seconds}s"
    else:
        return f"{value}s"


def bytes_to_human_readable(num_bytes, pos=None):
    """
    Convert a number of bytes into a human-readable string.
    E.g., 1048576 -> '1.0 MB'.
    """
    # Use 1024-based units
    for unit in ["B", "KB", "MB", "GB", "TB", "PB"]:
        if abs(num_bytes) < 1024.0:
            return f"{num_bytes:3.1f} {unit}"
        num_bytes /= 1024.0
    return f"{num_bytes:.1f} PB"


def plot_rule_data(combined_file, sample_unit_size, sample_size):
    """
    Given a combined benchmark TSV (one rule's data),
    parse the columns for 'filename', 'max_rss' and 's'.
    Determine the combined FASTQ size by matching 'filename'
      - If 'filename' has a dash, interpret it as 'sample-unit'
      - Otherwise interpret it as 'sample'
    Then create two scatter subplots:
      1) x-axis = combined FASTQ size, y-axis = max_rss
      2) x-axis = combined FASTQ size, y-axis = s
    Returns a matplotlib figure object (and data arrays if needed).
    """

    # We'll parse the file in a straightforward manner:
    # Read the header, find the columns of interest, then parse each line.
    x_vals_maxrss = []
    y_vals_maxrss = []
    x_vals_s = []
    y_vals_s = []

    with open(combined_file, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        col_idx = {name: i for i, name in enumerate(header)}

        # We expect at least: filename, max_rss, s
        for required in ["filename", "max_rss", "s"]:
            if required not in col_idx:
                # If the file doesn't have these columns, skip plotting
                print("Warning: Wrong columns in", combined_file)
                return None

        for line in f:
            row = line.rstrip("\n").split("\t")
            filename_val = row[col_idx["filename"]]
            maxrss_val_str = row[col_idx["max_rss"]]
            s_val_str = row[col_idx["s"]]

            # If max_rss or s is blank (or not numeric), skip
            try:
                maxrss_val = float(maxrss_val_str)
                s_val = float(s_val_str)
            except ValueError:
                continue

            # Above, we created the "filename" column of the summary table by using the rule name,
            # and the sample (and unit if present). As this script here has no knowledge
            # of our internal snakemake rules (and should not have, to keep it independent),
            # we hence need to take those appart again.
            # Figure out if either "sample" or "sample-unit" appear in the filename.
            # We hence only assume that we always use a dash in our rule wildcards.
            fastq_size = [ sample_unit_size[key] for key in sample_unit_size if key in filename_val ]
            if not fastq_size:
                fastq_size = [ sample_size[key] for key in sample_size if key in filename_val ]
            if not fastq_size:
                continue  # no match found
            if len(fastq_size) > 1:
                print("Warning: Multiple entries for", filename_val, "in", combined_file)
            fastq_size = fastq_size[0]

            x_vals_maxrss.append(fastq_size)
            y_vals_maxrss.append(maxrss_val * 1024 * 1024)
            x_vals_s.append(fastq_size)
            y_vals_s.append(s_val)

    # If no data was collected, return None
    if not x_vals_maxrss:
        return None

    # Make a figure with 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Plot 1: s vs combined FASTQ size
    axes[0].scatter(x_vals_s, y_vals_s, alpha=0.7)
    axes[0].set_xlabel("Sample size")
    axes[0].set_ylabel("Runtime")
    axes[0].set_title("Sample FASTQ size vs runtime")
    axes[0].xaxis.set_major_formatter(ticker.FuncFormatter(bytes_to_human_readable))
    axes[0].yaxis.set_major_formatter(ticker.FuncFormatter(seconds_to_hms))

    # Plot 2: max_rss vs combined FASTQ size
    axes[1].scatter(x_vals_maxrss, y_vals_maxrss, alpha=0.7)
    axes[1].set_xlabel("Sample size")
    axes[1].set_ylabel("Memory")
    axes[1].set_title("Sample FASTQ size vs memory")
    axes[1].xaxis.set_major_formatter(ticker.FuncFormatter(bytes_to_human_readable))
    axes[1].yaxis.set_major_formatter(ticker.FuncFormatter(bytes_to_human_readable))

    fig.tight_layout()
    return fig


def main():
    if len(sys.argv) != 2:
        print("Usage: python collect-benchmarks.py <analysis_dir>")
        sys.exit(1)

    # Get the paths as needed
    analysis_dir = sys.argv[1]
    benchmarks_dir = os.path.join(analysis_dir, "benchmarks")
    output_dir = os.path.join(analysis_dir, "benchmarks-summary")

    # Check that they are in order
    if not os.path.isdir(analysis_dir):
        print(f"Error: '{analysis_dir}' is not a valid directory.")
        sys.exit(1)
    if not os.path.isdir(benchmarks_dir):
        print(
            f"Error: Benchmarks directory '{benchmarks_dir}' is not a valid directory.",
            "Did you run grenepipe here before?"
        )
        sys.exit(1)

    # 0) Read the config file
    samples_path = get_samples_table_path(os.path.join( analysis_dir, "config.yaml" ))

    # 1) Read & parse the samples table to get FASTQ sizes
    sample_unit_size, sample_size = parse_samples_table(samples_path)
    print("Parsed sample FASTQ sizes:")
    print("  sample_unit_size:", dict(list(sample_unit_size.items())[:5]), "...")
    print("  sample_size:", dict(list(sample_size.items())[:5]), "...")

    # 2) Combine benchmarks
    combined_files = combine_benchmarks(benchmarks_dir, output_dir)

    # 3) For each combined file, create scatter plots
    for cf in combined_files:
        fig = plot_rule_data(cf, sample_unit_size, sample_size)
        if fig is not None:
            # Name the plot based on the TSV filename,
            # e.g. "ruleA_subdir.tsv" -> "ruleA_subdir_plots.png"
            base_name = os.path.splitext(os.path.basename(cf))[0]
            plot_path = os.path.join(output_dir, base_name + "_plots.png")
            fig.savefig(plot_path, dpi=150)
            plt.close(fig)  # free memory
            print(f"Created scatter plot: {plot_path}")
        else:
            print(f"Skipping plot for {cf} (missing columns or no data).")


if __name__ == "__main__":
    main()
