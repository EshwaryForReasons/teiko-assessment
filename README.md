# Loblaw Bio Immune Cell Analysis

This repository analyzes immune cell count data from `cell-count.csv` for Loblaw Bio's clinical trial. The project loads the CSV into a normalized SQLite database, computes relative immune cell frequencies, compares miraclib responders versus non-responders in melanoma PBMC samples, summarizes a baseline treatment subset, and serves the results in an interactive Streamlit dashboard.

## Dashboard

After starting the dashboard with `make dashboard`, open:

[http://localhost:8501](http://localhost:8501)

When running in GitHub Codespaces, open the forwarded port `8501`.

## Repository Structure

```text
.
├── cell-count.csv                                              # Input data file
├── load_data.py                                                # Part 1: database schema creation and CSV loading
├── analysis.py                                                 # Parts 2-4: analysis tables, statistics, and static plot generation
├── dashboard.py                                                # Interactive Streamlit dashboard
├── requirements.txt                                            # Python dependencies
├── Makefile                                                    # setup, pipeline, and dashboard commands
├── cell_counts.db                                              # Generated SQLite database
└── out/                                                        # Generated analysis outputs
    ├── cell_frequencies.csv                                    # Generated Part 2 output
    ├── miraclib_melanoma_pbmc_frequencies.csv                  # Generated Part 3 analysis input table
    ├── miraclib_melanoma_pbmc_response_stats.csv               # Generated Part 3 statistical results
    ├── miraclib_melanoma_pbmc_response_boxplot.png             # Generated Part 3 boxplot
    ├── baseline_miraclib_melanoma_pbmc_samples.csv             # Generated Part 4 baseline subset
    ├── baseline_miraclib_melanoma_pbmc_project_counts.csv
    ├── baseline_miraclib_melanoma_pbmc_response_counts.csv
    ├── baseline_miraclib_melanoma_pbmc_gender_counts.csv
    └── melanoma_male_responder_baseline_b_cell_average.csv
```

Generated files are created by `make pipeline`.

## How to Run in GitHub Codespaces

Place `cell-count.csv` in the repository root, then run:

```bash
make setup
make pipeline
make dashboard
```

The commands do the following:

```bash
make setup
```

Installs all Python dependencies from `requirements.txt`.

```bash
make pipeline
```

Runs the complete data pipeline from start to finish. This creates `cell_counts.db`, loads all rows from `cell-count.csv`, generates the Part 2 cell frequency table, runs the Part 3 responder versus non-responder statistical analysis, creates the Part 3 boxplot, and generates the Part 4 baseline subset summaries.

```bash
make dashboard
```

Starts the Streamlit dashboard on port `8501`.

## Relational Database Schema

The SQLite database uses a normalized schema with four core tables and one convenience view.

### `subjects`

One row per patient subject.

Columns:
- `subject_id`: subject identifier from the CSV `subject` column
- `project`: project identifier
- `indication`: disease indication, such as melanoma, carcinoma, or healthy
- `age`: subject age
- `gender`: subject sex from the CSV `sex` column
- `treatment`: treatment group, such as miraclib, phauximab, or none
- `response`: treatment response, stored as text and allowed to be `NULL` for healthy controls or missing values

### `samples`

One row per biological sample.

Columns:
- `sample_id`: sample identifier from the CSV `sample` column
- `subject_id`: foreign key to `subjects.subject_id`
- `sample_type`: sample type, such as PBMC
- `time_from_treatment_start`: timepoint relative to treatment start

### `cell_populations`

One row per immune cell population.

Columns:
- `population_id`: integer primary key
- `population_name`: immune cell population name, such as `b_cell`, `cd8_t_cell`, `cd4_t_cell`, `nk_cell`, or `monocyte`

### `cell_counts`

Long-form fact table containing one row per sample and cell population.

Columns:
- `sample_id`: foreign key to `samples.sample_id`
- `population_id`: foreign key to `cell_populations.population_id`
- `cell_count`: observed cell count for that population in that sample

The primary key is `(sample_id, population_id)`, ensuring each sample has at most one count for each population.

### `cell_counts_long`

A view joining the subject metadata, sample metadata, cell population names, and cell counts. This view is useful for downstream analytics and dashboards because it presents the normalized database in an analysis-ready long format.

## Schema Rationale and Scalability

The data starts as a wide CSV with one column per cell population. The database stores cell counts in long format instead. This design avoids hardcoding the five current immune cell populations into the database schema. If future datasets add more populations, such as additional T cell subsets or myeloid populations, the database can add rows to `cell_populations` without adding new columns to the schema.

This normalized structure also scales better for analytics. With hundreds of projects and thousands of samples, queries can filter by project, indication, treatment, response, sample type, timepoint, and population using joins instead of repeatedly reshaping wide tables. For larger deployments, indexes could be added on frequently filtered columns such as `subjects.project`, `subjects.indication`, `subjects.treatment`, `subjects.response`, `samples.sample_type`, `samples.time_from_treatment_start`, and `cell_counts.population_id`. The same schema can also support additional analytics such as longitudinal modeling, response prediction, baseline-only comparisons, project-level summaries, and dashboard filters.

## Analysis Overview (brief code structure rationale)

### Part 1: Data Management

`load_data.py` creates the SQLite database `cell_counts.db`, initializes the schema, and loads all rows from `cell-count.csv`. It stores blank response values as `NULL`, which is important for healthy subjects where response is not applicable.

### Part 2: Initial Analysis

`analysis.py` computes the relative frequency of each immune cell population in each sample. For each sample, it sums the five cell population counts to get `total_count`, then calculates each population's percentage of that total. The output is stored as a .csv file matching the columns requested in the assessment instructions.

### Part 3: Statistical Analysis

The program filters to melanoma PBMC samples from subjects treated with miraclib and compares responders (`response = yes`) against non-responders (`response = no`). It visualizes the cell population relative frequencies using boxplots and reports statistical evidence for differences between groups.

After some research, I found a combination of Mann-Whitney and Benjamini-Hochberg tests is good practice here. Thus, the statistical test is a two-sided Mann-Whitney U test for each immune cell population, followed by Benjamini-Hochberg correction across the five populations to control the false discovery rate.

### Part 4: Data Subset Analysis

The program identifies baseline melanoma PBMC samples where `time_from_treatment_start = 0` and treatment is miraclib. It then reports:
- number of samples from each project
- number of responder and non-responder subjects
- number of male and female subjects