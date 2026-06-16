import sqlite3
import pandas as pd
import plotly.express as px
import streamlit as st
from dataclasses import dataclass
from analysis import (
    paths,
    get_baseline_miraclib_melanoma_pbmc_subset,
    get_baseline_miraclib_melanoma_pbmc_summaries,
    get_cell_frequency_table,
    get_miraclib_melanoma_pbmc_frequencies,
    compute_response_statistics,
    get_melanoma_male_responder_baseline_b_cell_average,
    get_sample_metadata_table
)

@dataclass
class DashboardTables:
    frequency_table: pd.DataFrame
    sample_metadata: pd.DataFrame
    response_frequencies: pd.DataFrame
    response_stats: pd.DataFrame
    baseline_subset: pd.DataFrame
    baseline_project_counts: pd.DataFrame
    baseline_response_counts: pd.DataFrame
    baseline_gender_counts: pd.DataFrame
    melanoma_male_responder_baseline_b_cell_average: pd.DataFrame

st.set_page_config( page_title="Loblaw Bio Immune Cell Analysis", layout="wide")
st.title("Loblaw Bio Immune Cell Analysis")
st.caption(
    "Interactive dashboard for immune cell population frequencies, "
    "miraclib response analysis, and baseline melanoma PBMC subset summaries."
)

if not paths.db_filepath.exists():
    st.error(f"Could not find {paths.db_filepath}. Run `make pipeline` before launching the dashboard.")
    st.stop()

@st.cache_data
def load_dashboard_tables(db_modified_time: float) -> DashboardTables:
    with sqlite3.connect(paths.db_filepath) as conn:
        frequency_table = get_cell_frequency_table(conn)
        response_frequencies = get_miraclib_melanoma_pbmc_frequencies(conn)
        baseline_subset = get_baseline_miraclib_melanoma_pbmc_subset(conn)
        baseline_project_counts, baseline_response_counts, baseline_gender_counts = get_baseline_miraclib_melanoma_pbmc_summaries(conn)
        melanoma_male_responder_baseline_b_cell_average = get_melanoma_male_responder_baseline_b_cell_average(conn)
        sample_metadata = get_sample_metadata_table(conn)

    response_stats = compute_response_statistics(response_frequencies)
    return DashboardTables(
        frequency_table=frequency_table,
        sample_metadata=sample_metadata,
        response_frequencies=response_frequencies,
        response_stats=response_stats,
        baseline_subset=baseline_subset,
        baseline_project_counts=baseline_project_counts,
        baseline_response_counts=baseline_response_counts,
        baseline_gender_counts=baseline_gender_counts,
        melanoma_male_responder_baseline_b_cell_average=melanoma_male_responder_baseline_b_cell_average,
    )

dashboard_tables = load_dashboard_tables(paths.db_filepath.stat().st_mtime)

# =====================
# summary of everything
# =====================

st.subheader("Executive Summary")

summary_col_1, summary_col_2, summary_col_3, summary_col_4 = st.columns(4)
with summary_col_1:
    st.metric("Total samples", dashboard_tables.frequency_table["sample"].nunique())
with summary_col_2:
    st.metric(
        "Miraclib melanoma PBMC samples",
        dashboard_tables.response_frequencies["sample"].nunique(),
    )
with summary_col_3:
    st.metric(
        "Baseline miraclib melanoma PBMC samples",
        dashboard_tables.baseline_subset["sample"].nunique(),
    )
with summary_col_4:
    significant_count = dashboard_tables.response_stats[
        dashboard_tables.response_stats["significant_at_fdr_0_05"]
    ].shape[0]
    st.metric("FDR-significant populations", significant_count)
    
b_cell_answer = (
    dashboard_tables
    .melanoma_male_responder_baseline_b_cell_average
    .loc[0, "average_b_cell_count"]
)
st.metric(
    "Avg B cells: melanoma male responders at time=0",
    f"{b_cell_answer:.2f}",
)

# =====================
# Tabs for each part
# =====================

part_2_tab, part_3_tab, part_4_tab = st.tabs(
    [
        "Part 2: Cell Frequencies",
        "Part 3: Response Analysis",
        "Part 4: Baseline Subset",
    ]
)

with part_2_tab:
    st.header("Relative frequency of each cell population in each sample")

    populations = sorted(dashboard_tables.frequency_table["population"].unique())
    selected_populations = st.multiselect(
        "Cell populations",
        options=populations,
        default=populations,
    )

    sample_metadata = dashboard_tables.sample_metadata.copy()
    with st.expander("Filter samples", expanded=True):
        filter_col_1, filter_col_2, filter_col_3, filter_col_4 = st.columns(4)
        with filter_col_1:
            selected_projects = st.multiselect(
                "Project",
                options=sorted(sample_metadata["project"].dropna().unique()),
                default=sorted(sample_metadata["project"].dropna().unique()),
            )
            selected_indications = st.multiselect(
                "Indication",
                options=sorted(sample_metadata["indication"].dropna().unique()),
                default=sorted(sample_metadata["indication"].dropna().unique()),
            )

        with filter_col_2:
            selected_treatments = st.multiselect(
                "Treatment",
                options=sorted(sample_metadata["treatment"].dropna().unique()),
                default=sorted(sample_metadata["treatment"].dropna().unique()),
            )
            selected_responses = st.multiselect(
                "Response",
                options=sorted(sample_metadata["response"].dropna().unique()),
                default=sorted(sample_metadata["response"].dropna().unique()),
            )

        with filter_col_3:
            selected_genders = st.multiselect(
                "Gender",
                options=sorted(sample_metadata["gender"].dropna().unique()),
                default=sorted(sample_metadata["gender"].dropna().unique()),
            )
            selected_sample_types = st.multiselect(
                "Sample type",
                options=sorted(sample_metadata["sample_type"].dropna().unique()),
                default=sorted(sample_metadata["sample_type"].dropna().unique()),
            )

        with filter_col_4:
            selected_timepoints = st.multiselect(
                "Time from treatment start",
                options=sorted(sample_metadata["time_from_treatment_start"].dropna().unique()),
                default=sorted(sample_metadata["time_from_treatment_start"].dropna().unique()),
            )

    filtered_sample_metadata = sample_metadata[
        sample_metadata["project"].isin(selected_projects)
        & sample_metadata["indication"].isin(selected_indications)
        & sample_metadata["treatment"].isin(selected_treatments)
        & sample_metadata["response"].isin(selected_responses)
        & sample_metadata["gender"].isin(selected_genders)
        & sample_metadata["sample_type"].isin(selected_sample_types)
        & sample_metadata["time_from_treatment_start"].isin(selected_timepoints)
    ]

    st.caption(f"{len(filtered_sample_metadata):,} samples match the selected metadata filters.")
    selection_mode = st.radio(
        "Choose samples to display",
        options=[
            "First N matching samples",
            "Random N matching samples",
            "Paste sample IDs",
            "All matching samples",
        ],
        horizontal=True,
    )

    if filtered_sample_metadata.empty:
        selected_samples = []
        st.info("No samples match the selected filters.")
    else:
        if selection_mode == "First N matching samples":
            n_samples = st.number_input(
                "Number of matching samples to include",
                min_value=1,
                max_value=len(filtered_sample_metadata),
                value=min(50, len(filtered_sample_metadata)),
                step=10,
            )
            selected_samples = (
                filtered_sample_metadata["sample"]
                .head(n_samples)
                .tolist()
            )
        elif selection_mode == "Random N matching samples":
            n_samples = st.number_input(
                "Number of random matching samples to include",
                min_value=1,
                max_value=len(filtered_sample_metadata),
                value=min(50, len(filtered_sample_metadata)),
                step=10,
            )
            selected_samples = (
                filtered_sample_metadata["sample"]
                .sample(n=n_samples, random_state=42)
                .tolist()
            )
        elif selection_mode == "Paste sample IDs":
            pasted_sample_ids = st.text_area(
                "Paste sample IDs, separated by commas, spaces, or new lines",
                height=100,
            )
            raw_sample_ids = (
                pasted_sample_ids.replace(",", "\n")
                .replace(" ", "\n")
                .splitlines()
            )
            selected_samples = [
                sample_id.strip()
                for sample_id in raw_sample_ids
                if sample_id.strip()
            ]
            valid_sample_ids = set(filtered_sample_metadata["sample"])
            missing_sample_ids = [
                sample_id
                for sample_id in selected_samples
                if sample_id not in valid_sample_ids
            ]
            if missing_sample_ids:
                st.warning(
                    "Some pasted sample IDs do not match the current filters: "
                    + ", ".join(missing_sample_ids[:10])
                )
            selected_samples = [
                sample_id
                for sample_id in selected_samples
                if sample_id in valid_sample_ids
            ]
        else:
            selected_samples = filtered_sample_metadata["sample"].tolist()

        filtered_frequency_table = dashboard_tables.frequency_table[
            dashboard_tables.frequency_table["population"].isin(selected_populations)
            & dashboard_tables.frequency_table["sample"].isin(selected_samples)
        ]

        st.dataframe(filtered_frequency_table, use_container_width=True)
        st.download_button(
            label="Download selected frequency table",
            data=filtered_frequency_table.to_csv(index=False),
            file_name="selected_cell_frequencies.csv",
            mime="text/csv",
        )
        
        # If the user has selected a large number of samples, plotting all of them will probably crash the dashboard
        # and be unreadable even if it doesn't. The following code avoid's that issue.
        MAX_PLOTTED_SAMPLES = 200
        plotted_samples = selected_samples[:MAX_PLOTTED_SAMPLES]
        plot_frequency_table = filtered_frequency_table[
            filtered_frequency_table["sample"].isin(plotted_samples)
        ]

        if len(selected_samples) > MAX_PLOTTED_SAMPLES:
            st.warning(
                f"{len(selected_samples):,} samples are selected. "
                f"Only the first {MAX_PLOTTED_SAMPLES:,} are shown in the plot, "
                "but the table and download include the full selected set."
            )

        if not plot_frequency_table.empty:
            figure = px.bar(
                plot_frequency_table,
                x="sample",
                y="percentage",
                color="population",
                barmode="stack",
                hover_data=["total_count", "count"],
                title="Cell population relative frequencies by sample",
            )

            figure.update_layout(
                yaxis_title="Relative frequency (%)",
                xaxis_title="Sample",
                yaxis=dict(range=[0, 100]),
            )

            st.plotly_chart(figure, use_container_width=True)
            
with part_3_tab:
    st.header("Miraclib response analysis in melanoma PBMC samples")
    st.write(
        "This analysis compares relative cell population frequencies between "
        "miraclib-treated melanoma PBMC responder and non-responder samples."
    )


    st.subheader("Responder versus non-responder boxplot")
    response_boxplot = px.box(
        dashboard_tables.response_frequencies,
        x="population",
        y="percentage",
        color="response",
        points="all",
        hover_data=["subject_id", "sample", "time_from_treatment_start"],
        title="Relative frequencies by response group",
    )
    response_boxplot.update_layout(
        yaxis_title="Relative frequency (%)",
        xaxis_title="Cell population",
    )
    st.plotly_chart(response_boxplot, use_container_width=True)


    st.subheader("Statistical test results")
    st.dataframe(dashboard_tables.response_stats, use_container_width=True)
    significant_populations = dashboard_tables.response_stats[
        dashboard_tables.response_stats["significant_at_fdr_0_05"]
    ]["population"].tolist()
    if significant_populations:
        st.success(
            "Significant populations after Benjamini-Hochberg correction: "
            + ", ".join(significant_populations)
        )
    else:
        st.info("No populations were significant after Benjamini-Hochberg correction.")
    st.caption(
        "The statistical comparison uses a two-sided Mann-Whitney U test for each "
        "cell population, followed by Benjamini-Hochberg false-discovery-rate correction."
    )
    
    st.subheader("Key finding")
    significant_stats = dashboard_tables.response_stats[
        dashboard_tables.response_stats["significant_at_fdr_0_05"]
    ]
    if significant_stats.empty:
        st.info(
            "No immune cell populations showed a statistically significant difference "
            "between responders and non-responders after Benjamini-Hochberg correction."
        )
    else:
        for _, row in significant_stats.iterrows():
            direction = (
                "higher in responders"
                if row["mean_difference_percentage_points"] > 0
                else "lower in responders"
            )
            st.success(
                f"{row['population']} was significantly different between groups "
                f"and was {direction}. "
                f"Adjusted p-value: {row['p_value_adj_bh']:.4g}."
            )

with part_4_tab:
    st.header("Baseline melanoma PBMC samples treated with miraclib")
    st.write(
        "This subset includes melanoma PBMC samples where time_from_treatment_start is 0 "
        "and the patient treatment is miraclib."
    )

    st.subheader("Matching baseline samples")
    st.dataframe(dashboard_tables.baseline_subset, use_container_width=True)
    part_4_col_1, part_4_col_2, part_4_col_3 = st.columns(3)
    with part_4_col_1:
        st.metric("Baseline samples", dashboard_tables.baseline_subset["sample"].nunique())
    with part_4_col_2:
        st.metric("Baseline subjects", dashboard_tables.baseline_subset["subject_id"].nunique())
    with part_4_col_3:
        st.metric("Projects represented", dashboard_tables.baseline_subset["project"].nunique())

    project_column, response_column, gender_column = st.columns(3)
    with project_column:
        st.subheader("Samples by project")
        st.dataframe(dashboard_tables.baseline_project_counts, use_container_width=True)
        if not dashboard_tables.baseline_project_counts.empty:
            project_figure = px.bar(
                dashboard_tables.baseline_project_counts,
                x="project",
                y="n_samples",
                title="Baseline sample count by project",
            )
            st.plotly_chart(project_figure, use_container_width=True)

    with response_column:
        st.subheader("Subjects by response")
        st.dataframe(dashboard_tables.baseline_response_counts, use_container_width=True)
        if not dashboard_tables.baseline_response_counts.empty:
            response_figure = px.bar(
                dashboard_tables.baseline_response_counts,
                x="response",
                y="n_subjects",
                title="Baseline subject count by response",
            )
            st.plotly_chart(response_figure, use_container_width=True)

    with gender_column:
        st.subheader("Subjects by gender")
        st.dataframe(dashboard_tables.baseline_gender_counts, use_container_width=True)
        if not dashboard_tables.baseline_gender_counts.empty:
            gender_figure = px.bar(
                dashboard_tables.baseline_gender_counts,
                x="gender",
                y="n_subjects",
                title="Baseline subject count by gender",
            )
            st.plotly_chart(gender_figure, use_container_width=True)