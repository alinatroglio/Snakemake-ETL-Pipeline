# Snakemake-ETL-Pipeline

This repository contains a reproducible, modular ETL (Extract–Transform–Load) pipeline built using [Snakemake](https://snakemake.readthedocs.io/en/stable/) as the workflow manager, [PostgreSQL](https://www.postgresql.org/) for storage, and [Streamlit](https://streamlit.io/) for interactive data visualization. The project is fully containerized with Docker Compose, ensuring that the entire setup runs seamlessly across systems.

## Project Structure
```bash
├── Dockerfile                     # Defines Python image and installs dependencies
├── docker-compose.yml             # ETL, DB, and Dashboard containers
├── requirements.txt               # Python dependencies 
├── Snakefile                      # Snakemake workflow with defined rules
│
├── config/
│   ├── config.yaml                # Main configuration (env variables, mappings)
│   └── mappings/
│       ├── encounter_type_map.yaml
│       ├── gender_map.yaml
│       └── id_map.yaml
│
├── data/                          # Input datasets
│
├── output/                        # Generated outputs (transformed & merged data)
│   ├── transformed/
│   ├── merged/
│   └── ok_files/                  # Completion marker files
│
├── scripts/
│   ├── load_data.py               # Reads and pre-processes raw data
│   ├── transform_data.py          # Helper functions for transformation logic
│   ├── merge_data.py              # Combines transformed datasets
│   ├── load_to_database.py        # Loads merged data into PostgreSQL
│   ├── etl_logging.py             # Logging for data quality report
│   ├── config_utils.py            # Helpers to read YAML configuration
│   └── dashboard.py               # Streamlit dashboard for data overview
│
└── shared/                        # Shared folder for logs and intermediate files
```

## How it works
### 1. Workflow Management
The Snakefile defines a set of ETL rules, each specifying:
- Input: source files or outputs from previous rules
- Output: processed files or database-ready datasets
- Script: the Python script that performs the steps
  
This modular structure ensures each step can be independently tested, cached, and re-run only when inputs change.
### 2. Container Setup
The system runs three containers:
| Service   | Description            | Port        |
|------------|------------------------|--------------|
| `db`       | PostgreSQL database    | 5430 → 5432  |
| `etl`      | Snakemake ETL pipeline |              |
| `streamlit`| Streamlit dashboard    | 8501         |

All containers share volumes for data and logs, ensuring the results are accessible for the dashboard.
### 3. Build and Run
#### Prerequisites
Before running the project, make sure you have the following installed:

- [Docker Desktop](https://docs.docker.com/desktop/setup/install/windows-install/) (or Docker Engine + Docker Compose if you’re on Linux)
- [Git](https://git-scm.com/downloads/win) – to clone the repository

Once Docker Desktop is running, you can proceed in the repository folder with:
```bash
docker compose build
docker compose up
```
When running for the first time, Snakemake executes the ETL workflow automatically.

### 4. Access the Dashboard
Once containers are up, to see the interactive dashboard, open your browser and navigate to
[http://0.0.0.0:8501](http://127.0.0.1:8501/)

## Design Choices
- Snakemake for reproducibility, each step explicitly declares its inputs/outputs, reruns only what changed
- Docker for consistent environment on every machine
- Modular scripts with clear separation of extraction, transformation, and loading
- Easily adaptable to new datasets or database connections

## Known Limitations / Future Improvements
- **Data quality rules**: The logging framework captures basic inconsistencies but not semantic issues (e.g. discharge_dt before admit_dt).
- **Performance**: Steps use `DataFrame.apply` row-by-row. Vectorized would be more efficient.
- **Diagnoses ID/Encounter ID**: Currently, encounter IDs are assumed unique per record. In practice, one patient can have multiple diagnoses for the same encounter, so the merge logic may duplicate rows. This should be re-designed with a separate `diagnoses` table.
- **Data transformation**: Some transformation rules were developed quickly to illustrate the pipeline flow rather than full robustness. In future solutions, they should be re-implemented with consistent schema validation, reusable mapping logic, and externalized reference data.
