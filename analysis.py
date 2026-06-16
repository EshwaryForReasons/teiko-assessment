from pathlib import Path
import sqlite3

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
from paths_store import paths

CELL_POPULATIONS = [
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]

ALPHA = 0.05

BASELINE_FILTER_CTE = """
    WITH baseline_miraclib_melanoma_pbmc AS (
        SELECT
            s.project,
            s.subject_id,
            s.indication,
            s.age,
            s.gender,
            s.treatment,
            COALESCE(NULLIF(LOWER(TRIM(s.response)), ''), 'unknown') AS response,
            sm.sample_id AS sample,
            sm.sample_type,
            sm.time_from_treatment_start
        FROM samples sm
        JOIN subjects s
            ON sm.subject_id = s.subject_id
        WHERE
            LOWER(s.indication) = 'melanoma'
            AND LOWER(s.treatment) = 'miraclib'
            AND LOWER(sm.sample_type) = 'pbmc'
            AND sm.time_from_treatment_start = 0
    )
"""

# =====================
# SQL query functions
# =====================

def get_cell_frequency_table(conn: sqlite3.Connection) -> pd.DataFrame:
    query = """
        WITH sample_totals AS (
            SELECT
                sample_id,
                SUM(cell_count) AS total_count
            FROM cell_counts
            GROUP BY sample_id
        )

        SELECT
            cc.sample_id AS sample,
            st.total_count AS total_count,
            cp.population_name AS population,
            cc.cell_count AS count,
            ROUND(100.0 * cc.cell_count / st.total_count, 4) AS percentage
        FROM cell_counts cc
        JOIN sample_totals st
            ON cc.sample_id = st.sample_id
        JOIN cell_populations cp
            ON cc.population_id = cp.population_id
        ORDER BY
            cc.sample_id,
            cp.population_id;
    """
    return pd.read_sql_query(query, conn)

def get_miraclib_melanoma_pbmc_frequencies(conn: sqlite3.Connection) -> pd.DataFrame:
    query = """
        WITH sample_totals AS (
            SELECT
                sample_id,
                SUM(cell_count) AS total_count
            FROM cell_counts
            GROUP BY sample_id
        )

        SELECT
            s.subject_id,
            sm.sample_id AS sample,
            s.indication,
            s.treatment,
            LOWER(s.response) AS response,
            sm.sample_type,
            sm.time_from_treatment_start,
            st.total_count,
            cp.population_name AS population,
            cc.cell_count AS count,
            ROUND(100.0 * cc.cell_count / st.total_count, 4) AS percentage
        FROM cell_counts cc
        JOIN sample_totals st
            ON cc.sample_id = st.sample_id
        JOIN cell_populations cp
            ON cc.population_id = cp.population_id
        JOIN samples sm
            ON cc.sample_id = sm.sample_id
        JOIN subjects s
            ON sm.subject_id = s.subject_id
        WHERE
            LOWER(s.indication) = 'melanoma'
            AND LOWER(s.treatment) = 'miraclib'
            AND LOWER(sm.sample_type) = 'pbmc'
            AND LOWER(COALESCE(s.response, '')) IN ('yes', 'no')
        ORDER BY
            cp.population_id,
            s.response,
            sm.sample_id;
    """
    return pd.read_sql_query(query, conn)

def get_baseline_miraclib_melanoma_pbmc_subset(conn: sqlite3.Connection) -> pd.DataFrame:
    query = (
        BASELINE_FILTER_CTE
        + """
        SELECT
            project,
            subject_id,
            indication,
            age,
            gender,
            treatment,
            response,
            sample,
            sample_type,
            time_from_treatment_start
        FROM baseline_miraclib_melanoma_pbmc
        ORDER BY
            project,
            subject_id,
            sample;
        """
    )
    return pd.read_sql_query(query, conn)

def get_baseline_miraclib_melanoma_pbmc_summaries(conn: sqlite3.Connection) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    project_counts_query = (
        BASELINE_FILTER_CTE
        + """
        SELECT
            project,
            COUNT(DISTINCT sample) AS n_samples
        FROM baseline_miraclib_melanoma_pbmc
        GROUP BY project
        ORDER BY project;
        """
    )

    response_counts_query = (
        BASELINE_FILTER_CTE
        + """
        SELECT
            response,
            COUNT(DISTINCT subject_id) AS n_subjects
        FROM baseline_miraclib_melanoma_pbmc
        GROUP BY response
        ORDER BY response;
        """
    )

    gender_counts_query = (
        BASELINE_FILTER_CTE
        + """
        SELECT
            gender,
            COUNT(DISTINCT subject_id) AS n_subjects
        FROM baseline_miraclib_melanoma_pbmc
        GROUP BY gender
        ORDER BY gender;
        """
    )

    project_counts = pd.read_sql_query(project_counts_query, conn)
    response_counts = pd.read_sql_query(response_counts_query, conn)
    gender_counts = pd.read_sql_query(gender_counts_query, conn)
    return project_counts, response_counts, gender_counts

def get_melanoma_male_responder_baseline_b_cell_average(conn: sqlite3.Connection) -> pd.DataFrame:
    query = """
        SELECT
            ROUND(AVG(cc.cell_count), 2) AS average_b_cell_count,
            COUNT(DISTINCT sm.sample_id) AS n_samples
        FROM cell_counts cc
        JOIN cell_populations cp
            ON cc.population_id = cp.population_id
        JOIN samples sm
            ON cc.sample_id = sm.sample_id
        JOIN subjects s
            ON sm.subject_id = s.subject_id
        WHERE
            LOWER(s.indication) = 'melanoma'
            AND LOWER(TRIM(s.response)) = 'yes'
            AND LOWER(TRIM(s.gender)) IN ('m', 'male')
            AND sm.time_from_treatment_start = 0
            AND cp.population_name = 'b_cell';
    """
    return pd.read_sql_query(query, conn)

def get_sample_metadata_table(conn: sqlite3.Connection) -> pd.DataFrame:
    query = """
        SELECT
            sm.sample_id AS sample,
            s.project,
            s.subject_id,
            s.indication,
            s.treatment,
            COALESCE(NULLIF(LOWER(TRIM(s.response)), ''), 'unknown') AS response,
            s.gender,
            sm.sample_type,
            sm.time_from_treatment_start
        FROM samples sm
        JOIN subjects s
            ON sm.subject_id = s.subject_id
        ORDER BY
            s.project,
            s.subject_id,
            sm.sample_id;
    """
    return pd.read_sql_query(query, conn)

# =====================
# Compute statistics
# =====================

def compute_response_statistics(frequencies: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for population in CELL_POPULATIONS:
        population_data = frequencies[frequencies["population"] == population]

        responder_values = population_data[
            population_data["response"] == "yes"
        ]["percentage"].dropna()

        non_responder_values = population_data[
            population_data["response"] == "no"
        ]["percentage"].dropna()

        if len(responder_values) == 0 or len(non_responder_values) == 0:
            p_value = None
            u_statistic = None
            rank_biserial_correlation = None
        else:
            test_result = mannwhitneyu(
                responder_values,
                non_responder_values,
                alternative="two-sided",
            )

            p_value = test_result.pvalue
            u_statistic = test_result.statistic
            rank_biserial_correlation = (
                (2 * u_statistic) / (len(responder_values) * len(non_responder_values))
            ) - 1

        rows.append(
            {
                "population": population,
                "n_responder_samples": len(responder_values),
                "n_non_responder_samples": len(non_responder_values),
                "responder_mean_percentage": responder_values.mean(),
                "non_responder_mean_percentage": non_responder_values.mean(),
                "mean_difference_percentage_points": (
                    responder_values.mean() - non_responder_values.mean()
                ),
                "responder_median_percentage": responder_values.median(),
                "non_responder_median_percentage": non_responder_values.median(),
                "median_difference_percentage_points": (
                    responder_values.median() - non_responder_values.median()
                ),
                "mann_whitney_u_statistic": u_statistic,
                "rank_biserial_correlation": rank_biserial_correlation,
                "p_value": p_value,
            }
        )

    stats = pd.DataFrame(rows)
    stats["p_value_adj_bh"] = pd.NA

    valid_p_values = stats["p_value"].notna()
    if valid_p_values.any():
        sorted_indices = stats.loc[valid_p_values, "p_value"].sort_values().index
        num_tests = len(sorted_indices)
        previous_adjusted_p_value = 1.0

        for rank, index in reversed(list(enumerate(sorted_indices, start=1))):
            adjusted_p_value = stats.loc[index, "p_value"] * num_tests / rank
            adjusted_p_value = min(adjusted_p_value, previous_adjusted_p_value, 1.0)

            stats.loc[index, "p_value_adj_bh"] = adjusted_p_value
            previous_adjusted_p_value = adjusted_p_value

    stats["significant_at_fdr_0_05"] = (
        stats["p_value_adj_bh"].fillna(1.0).astype(float) < ALPHA
    )

    return stats.sort_values("p_value_adj_bh")

# =====================
# Offline (no dashboard) visualization
# =====================

def plot_response_boxplots(frequencies: pd.DataFrame, output_path: Path) -> None:
    boxplot_data = []
    boxplot_labels = []

    for population in CELL_POPULATIONS:
        population_data = frequencies[frequencies["population"] == population]

        responder_values = population_data[
            population_data["response"] == "yes"
        ]["percentage"].dropna()

        non_responder_values = population_data[
            population_data["response"] == "no"
        ]["percentage"].dropna()

        boxplot_data.append(responder_values)
        boxplot_labels.append(f"{population}\nresponder")

        boxplot_data.append(non_responder_values)
        boxplot_labels.append(f"{population}\nnon-responder")

    plt.figure(figsize=(14, 7))
    plt.boxplot(boxplot_data, tick_labels=boxplot_labels, showmeans=True)
    plt.ylabel("Relative frequency (%)")
    plt.xlabel("Immune cell population and response group")
    plt.title(
        "Relative Immune Cell Frequencies in Melanoma PBMC Samples Treated with Miraclib"
    )
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def main() -> None:
    if not paths.db_filepath.exists():
        raise FileNotFoundError(
            "Could not find cell_counts.db. Run `python load_data.py` first."
        )

    # Perform analysis for all parts at once to ensure only database connection is only required once.
    with sqlite3.connect(paths.db_filepath) as conn:
        # part 2--get the frequency table and save it to CSV with the schema
        # requested in the technical assessment instructions.
        frequency_table = get_cell_frequency_table(conn)
        response_frequencies = get_miraclib_melanoma_pbmc_frequencies(conn)
        baseline_subset = get_baseline_miraclib_melanoma_pbmc_subset(conn)
        baseline_project_counts, baseline_response_counts, baseline_gender_counts = get_baseline_miraclib_melanoma_pbmc_summaries(conn)
        melanoma_male_responder_baseline_b_cell_average = get_melanoma_male_responder_baseline_b_cell_average(conn)

    # Ensure we have results. Ideally, this should have some fallback mechanism or allow for graceful failure.
    # For the sake of this technical assessment, I have decided to simply raise an error.
    if response_frequencies.empty:
        raise ValueError(
            "No melanoma PBMC samples treated with miraclib and response yes/no were found."
        )
    if baseline_subset.empty:
        raise ValueError(
            "No baseline melanoma PBMC samples treated with miraclib were found."
        )

    response_stats = compute_response_statistics(response_frequencies)

    # Save everything as .csv files for the interactive dashboard and potentially for separate offline analysis.
    
    # Part 2
    frequency_table.to_csv(paths.cell_frequencies_output_filepath, index=False)
    # Part 3
    response_frequencies.to_csv(paths.response_analysis_input_output_filepath, index=False)
    response_stats.to_csv(paths.response_stats_output_filepath, index=False)
    # Part 4
    baseline_subset.to_csv(paths.baseline_subset_output_filepath, index=False)
    baseline_project_counts.to_csv(paths.baseline_project_counts_output_filepath, index=False)
    baseline_response_counts.to_csv(paths.baseline_response_counts_output_filepath, index=False)
    baseline_gender_counts.to_csv(paths.baseline_gender_counts_output_filepath, index=False)
    # Special question
    melanoma_male_responder_baseline_b_cell_average.to_csv(
        paths.melanoma_male_responder_baseline_b_cell_average_output_filepath,
        index=False,
    )
    
    # Running logs
    print("\n========== Running logs ==========")
    print(f"Saved Part 2 table to {paths.cell_frequencies_output_filepath.name}")
    print(f"Saved Part 3 analysis table to {paths.response_analysis_input_output_filepath.name}")
    print(f"Saved Part 3 statistics to {paths.response_stats_output_filepath.name}")
    print(f"Saved Part 4 baseline subset to {paths.baseline_subset_output_filepath.name}")
    print(f"Saved Part 4 project counts to {paths.baseline_project_counts_output_filepath.name}")
    print(f"Saved Part 4 response counts to {paths.baseline_response_counts_output_filepath.name}")
    print(f"Saved Part 4 gender counts to {paths.baseline_gender_counts_output_filepath.name}")
    print(
        "Saved special question output to "
        f"{paths.melanoma_male_responder_baseline_b_cell_average_output_filepath.name}"
    )


    print("\n========== Part 2 ==========")
    print("\nCell frequency table:")
    print(frequency_table[:5].to_string(index=False))
    print(
        "only showing the first 5 rows here for brevity, please use the interactive dashboard or check "
        f"the saved .csv file at {paths.cell_frequencies_output_filepath.name} for the full table."
    )


    print("\n========== Part 3 ==========")
    print(
        "My understanding is the interactive dashboard should provide all the requested abilities for Part 3.\n"
        "Nonetheless, have decided to provide some brief statistics and a boxplot for visualized to allow\n"
        "offline comparisons and immediate visualization of results."
    )
    print("\nMiraclib melanoma PBMC responder versus non-responder statistics:\n")
    print(response_stats.to_string(index=False))

    significant_populations = response_stats[
        response_stats["significant_at_fdr_0_05"]
    ]["population"].tolist()

    if significant_populations:
        print(
            "\nSignificant populations after Benjamini-Hochberg correction: "
            f"{', '.join(significant_populations)}"
        )
    else:
        print("\nNo populations were significant after Benjamini-Hochberg correction.")
        
    # Visualize in boxplot as requested in the technical assessment instructions.
    plot_response_boxplots(response_frequencies, paths.response_boxplot_output_filepath)
    print(f"Saved Part 3 boxplot to {paths.response_boxplot_output_filepath.name}")
    
    
    print("\n========== Part 4 ==========")
    print(
        "As before, I only print (at most) the first 5 rows of each subset here for brevity, but the interactive "
        "dashboard and saved .csv files contain the full results. Please use those for further analysis."
    )

    print("\nBaseline melanoma PBMC samples treated with miraclib:")
    print(baseline_subset[:5].to_string(index=False))

    print("\nBaseline sample counts by project:")
    print(baseline_project_counts[:5].to_string(index=False))

    print("\nBaseline subject counts by response:")
    print(baseline_response_counts[:5].to_string(index=False))

    print("\nBaseline subject counts by gender:")
    print(baseline_gender_counts[:5].to_string(index=False))
    
    
    print("\n========== Special Question ==========")
    average_b_cell_count = melanoma_male_responder_baseline_b_cell_average.loc[0, "average_b_cell_count"]
    n_samples = melanoma_male_responder_baseline_b_cell_average.loc[0, "n_samples"]

    print(
        "Considering melanoma males, the average number of B cells for responders\n"
        f"at time=0 is {average_b_cell_count:.2f} across {n_samples} samples."
        "\n\nNOTE: this does not carry forward the previous constraints of only using miraclib and PBMC samples.\n"
        "Since this is a separate question, I have performed the analysis on only the specified constraints of\n"
        "males with Melanoma responsive at time=0."
    )

if __name__ == "__main__":
    main()