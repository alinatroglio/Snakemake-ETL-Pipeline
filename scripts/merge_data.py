import pandas as pd
from snakemake.script import snakemake
from config_utils import load_mapping, map_values

# Data paths for merge
encounters_file = snakemake.input.transformed_encounters
diagnoses_file = snakemake.input.transformed_diagnoses
patients_file = snakemake.input.transformed_patients

# Load input data
df_encounters = pd.read_csv(encounters_file)
df_diagnoses = pd.read_csv(diagnoses_file)
df_patients = pd.read_csv(patients_file)

# Apply ID column name mapping 
id_map = load_mapping(snakemake.config["mappings"]["id"])
df_diagnoses = df_diagnoses.rename(columns=id_map)

# Merge encounter - patients - diagnoses
df_merged = (df_encounters
    .merge(df_patients, on="patient_id", how="outer")
    .merge(df_diagnoses, on="encounter_id", how="outer")
)

# Normalize identifier columns
ID_COLS = ["row_id", "encounter_id", "patient_id"]
for c in ID_COLS:
    if c in df_merged.columns:
        # Normalize whitespace and ensure pandas nullable string dtype
        df_merged[c] = df_merged[c].astype(str).str.strip().astype("string")

# Save merged dataset
df_merged.to_csv(snakemake.output.merged_data, index=False)