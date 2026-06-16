from pathlib import Path
import csv
import sqlite3
from paths_store import paths

CELL_POPULATIONS = [
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]

REQUIRED_COLUMNS = {
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
    *CELL_POPULATIONS,
}

def create_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        PRAGMA foreign_keys = ON;

        DROP VIEW IF EXISTS cell_counts_long;
        DROP TABLE IF EXISTS cell_counts;
        DROP TABLE IF EXISTS cell_populations;
        DROP TABLE IF EXISTS samples;
        DROP TABLE IF EXISTS subjects;

        CREATE TABLE subjects (
            subject_id TEXT PRIMARY KEY,
            project TEXT NOT NULL,
            indication TEXT NOT NULL,
            age INTEGER NOT NULL CHECK (age >= 0),
            gender TEXT NOT NULL,
            treatment TEXT NOT NULL,
            response TEXT
        );

        CREATE TABLE samples (
            sample_id TEXT PRIMARY KEY,
            subject_id TEXT NOT NULL,
            sample_type TEXT NOT NULL,
            time_from_treatment_start INTEGER NOT NULL CHECK (time_from_treatment_start >= 0),

            FOREIGN KEY (subject_id)
                REFERENCES subjects(subject_id)
                ON DELETE CASCADE
        );

        CREATE TABLE cell_populations (
            population_id INTEGER PRIMARY KEY AUTOINCREMENT,
            population_name TEXT NOT NULL UNIQUE
        );

        CREATE TABLE cell_counts (
            sample_id TEXT NOT NULL,
            population_id INTEGER NOT NULL,
            cell_count INTEGER NOT NULL CHECK (cell_count >= 0),

            PRIMARY KEY (sample_id, population_id),

            FOREIGN KEY (sample_id)
                REFERENCES samples(sample_id)
                ON DELETE CASCADE,

            FOREIGN KEY (population_id)
                REFERENCES cell_populations(population_id)
                ON DELETE CASCADE
        );

        CREATE VIEW cell_counts_long AS
        SELECT
            s.project,
            s.subject_id,
            s.indication,
            s.age,
            s.gender,
            s.treatment,
            s.response,
            sm.sample_id,
            sm.sample_type,
            sm.time_from_treatment_start,
            cp.population_name,
            cc.cell_count
        FROM cell_counts cc
        JOIN samples sm
            ON cc.sample_id = sm.sample_id
        JOIN subjects s
            ON sm.subject_id = s.subject_id
        JOIN cell_populations cp
            ON cc.population_id = cp.population_id;
        """
    )

def read_rows(data_path: Path) -> list[dict[str, str]]:
    with data_path.open("r", newline="", encoding="utf-8-sig") as file:
        sample = file.read(4096)
        file.seek(0)

        try:
            dialect = csv.Sniffer().sniff(sample)
        except csv.Error:
            dialect = csv.excel

        reader = csv.DictReader(file, dialect=dialect)
        rows = list(reader)

    if not rows:
        raise ValueError(f"No rows found in {data_path.name}.")

    missing_columns = REQUIRED_COLUMNS - set(rows[0].keys())
    if missing_columns:
        raise ValueError(
            f"Missing required columns in {data_path.name}: "
            f"{sorted(missing_columns)}"
        )

    return rows

def insert_data(conn: sqlite3.Connection, rows: list[dict[str, str]]) -> None:
    population_ids: dict[str, int] = {}

    for population in CELL_POPULATIONS:
        conn.execute(
            """
            INSERT INTO cell_populations (population_name)
            VALUES (?);
            """,
            (population,),
        )

        cursor = conn.execute(
            """
            SELECT population_id
            FROM cell_populations
            WHERE population_name = ?;
            """,
            (population,),
        )

        population_ids[population] = cursor.fetchone()[0]

    for row in rows:
        subject_id = row["subject"].strip()
        sample_id = row["sample"].strip()
        response = row["response"].strip().lower() or None

        if not subject_id:
            raise ValueError(f"Found row with empty subject ID: {row}")
        if not sample_id:
            raise ValueError(f"Found row with empty sample ID: {row}")

        conn.execute(
            """
            INSERT INTO subjects (
                subject_id,
                project,
                indication,
                age,
                gender,
                treatment,
                response
            )
            VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(subject_id) DO UPDATE SET
                project = excluded.project,
                indication = excluded.indication,
                age = excluded.age,
                gender = excluded.gender,
                treatment = excluded.treatment,
                response = excluded.response;
            """,
            (
                subject_id,
                row["project"].strip(),
                row["condition"].strip(),
                int(row["age"]),
                row["sex"].strip(),
                row["treatment"].strip(),
                response,
            ),
        )

        conn.execute(
            """
            INSERT INTO samples (
                sample_id,
                subject_id,
                sample_type,
                time_from_treatment_start
            )
            VALUES (?, ?, ?, ?)
            ON CONFLICT(sample_id) DO UPDATE SET
                subject_id = excluded.subject_id,
                sample_type = excluded.sample_type,
                time_from_treatment_start = excluded.time_from_treatment_start;
            """,
            (
                sample_id,
                subject_id,
                row["sample_type"].strip(),
                int(row["time_from_treatment_start"]),
            ),
        )

        for population in CELL_POPULATIONS:
            conn.execute(
                """
                INSERT INTO cell_counts (
                    sample_id,
                    population_id,
                    cell_count
                )
                VALUES (?, ?, ?)
                ON CONFLICT(sample_id, population_id) DO UPDATE SET
                    cell_count = excluded.cell_count;
                """,
                (
                    sample_id,
                    population_ids[population],
                    int(row[population]),
                ),
            )

def main() -> None:
    if not paths.csv_filepath.exists():
        raise FileNotFoundError("Could not find cell-count.csv in the repository root.")

    if paths.db_filepath.exists():
        paths.db_filepath.unlink()

    with sqlite3.connect(paths.db_filepath) as conn:
        conn.execute("PRAGMA foreign_keys = ON;")
        create_schema(conn)

        rows = read_rows(paths.csv_filepath)
        insert_data(conn, rows)

        conn.commit()

    print(f"Loaded {len(rows)} samples into {paths.db_filepath.name}")

if __name__ == "__main__":
    main()