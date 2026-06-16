"""
NOTE: In a real project it is, in my experience, common to have a separate file for paths.
I have attempted to emulate that behavior here. This is not necessary for a project
of a size this small, but it is a good practice nonetheless.
"""

from pathlib import Path
from dataclasses import dataclass

@dataclass(frozen=True)
class Paths:
    """
    I find a singleton design is common for this type of class. I have decided against that to avoid
    unnecessary complexity for such a small project.
    """
    # Standard dirs
    root_dir: Path
    out_dir: Path

    # Input / database filepaths
    db_filepath: Path
    csv_filepath: Path

    # Part 2 outputs
    cell_frequencies_output_filepath: Path

    # Part 3 outputs
    response_analysis_input_output_filepath: Path
    response_stats_output_filepath: Path
    response_boxplot_output_filepath: Path

    # Part 4 outputs
    baseline_subset_output_filepath: Path
    baseline_project_counts_output_filepath: Path
    baseline_response_counts_output_filepath: Path
    baseline_gender_counts_output_filepath: Path
    
    # Special question output
    melanoma_male_responder_baseline_b_cell_average_output_filepath: Path

    @staticmethod
    def create() -> "Paths":
        root_dir = Path(__file__).resolve().parent
        out_dir = root_dir / "out"
        out_dir.mkdir(exist_ok=True)

        db_filepath = root_dir / "cell_counts.db"
        csv_filepath = root_dir / "cell-count.csv"

        return Paths(
            root_dir=root_dir,
            out_dir=out_dir,
            db_filepath=db_filepath,
            csv_filepath=csv_filepath,
            cell_frequencies_output_filepath=out_dir / "cell_frequencies.csv",
            response_analysis_input_output_filepath=out_dir / "miraclib_melanoma_pbmc_frequencies.csv",
            response_stats_output_filepath=out_dir / "miraclib_melanoma_pbmc_response_stats.csv",
            response_boxplot_output_filepath=out_dir / "miraclib_melanoma_pbmc_response_boxplot.png",
            baseline_subset_output_filepath=out_dir / "baseline_miraclib_melanoma_pbmc_samples.csv",
            baseline_project_counts_output_filepath=out_dir / "baseline_miraclib_melanoma_pbmc_project_counts.csv",
            baseline_response_counts_output_filepath=out_dir / "baseline_miraclib_melanoma_pbmc_response_counts.csv",
            baseline_gender_counts_output_filepath=out_dir / "baseline_miraclib_melanoma_pbmc_gender_counts.csv",
            melanoma_male_responder_baseline_b_cell_average_output_filepath=(
                out_dir / "melanoma_male_responder_baseline_b_cell_average.csv"
            )
        )
        
paths = Paths.create()