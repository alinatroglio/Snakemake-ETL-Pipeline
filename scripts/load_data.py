import pandas as pd
from snakemake.script import snakemake
import csv
import re 
import pint 
from config_utils import load_mapping, map_values
from dateutil import parser
from transform_data import *
from etl_logging import reset_loggin

# Include maps for manual FHIR mapping of sex and encounter_type
FHIR_encounter_map = load_mapping(snakemake.config["mappings"]["encounter_type"])
FHIR_gender_map = load_mapping(snakemake.config["mappings"]["gender"])

# reset logging file
reset_loggin()

# Define standardized column groups used for consistent data cleaning:
# dates: columns expected to contain date/time values
# strip_safe: string columns that can be safely stripped of surrounding whitespace
# cap_safe: name-like columns suitable for capitalization normalization
dates = ["admit_dt", "discharge_dt", "dob", "recordedAt"]
strip_safe = [ "encounter_type", "source_file", "patient_id", "encounter_id", "given_name", "family_name", "sex", "weight", "height"]
cap_safe = ["given_name", "family_name"]

# Clean and standardize encounters.csv:
# 1. Detect and normalize the delimiter, load into a DataFrame
# 2. Unify and standardize datetime columns
# 3. Clean string fields 
# 4. Remove or fix duplicate encounter IDs
# 5. Map encounter types to FHIR codes
df_encounters = (
    clean_csv_to_df(snakemake.input.encounters)
    .pipe(unify_datetime, dates, file_name= snakemake.input.encounters)
    .pipe(clean_str_data, strip_safe)
    .pipe(fix_ids, "encounter_id", file_name= snakemake.input.encounters)
    .pipe(map_values, "encounter_type", "encounter_type_FHIR", FHIR_encounter_map))

# Clean and standardize the patients.csv:
# 1. Detect and normalize the delimiter, then load into a DataFrame
# 2. Unify and standardize datetime columns
# 3. Clean string fields 
# 4. Capitalize name fields consistently
# 5. Remove or fix duplicate patient IDs
# 6. Map sex values to FHIR codes
# 7. Normalize height and weight units to standard (cm, kg)
df_patients = (
    clean_csv_to_df(snakemake.input.patients)
    .pipe(unify_datetime, dates, file_name= snakemake.input.patients)
    .pipe(clean_str_data, strip_safe)
    .pipe(capitalize, cap_safe)
    .pipe(fix_ids, "patient_id", file_name= snakemake.input.patients)
    .pipe(map_values, "sex", "gender_FHIR", FHIR_gender_map))
df_patients["height_cm"] = df_patients.apply(lambda row: unify_height(row["height"], row.name, file_name= snakemake.input.encounters), axis=1)
df_patients["weight_kg"] = df_patients.apply(lambda row: unify_weight(row["weight"], row.name, file_name= snakemake.input.encounters ), axis=1)

# Clean and validate the diagnoses.xml:
# 1. Load diagnoses data from XML into a DataFrame
# 2. Standardize datetime columns
# 3. Remove or fix duplicate encounter IDs
# 4. Validate ICD-10 diagnosis codes and flag invalid entries
df_diagnoses = pd.read_xml(snakemake.input.diagnoses)
df_diagnoses = (df_diagnoses.pipe(unify_datetime, dates, file_name = snakemake.input.diagnoses)
                .pipe(fix_ids, "encounterId", file_name = snakemake.input.diagnoses))
df_diagnoses["valid_idc10_code"] = df_diagnoses.apply(lambda row: validate_icd10_code(row["code"], row.name, file_name=snakemake.input.diagnoses), axis=1)

# Save transformed DataFrames to csv
df_encounters.to_csv(snakemake.output.transformed, index=False)
df_patients.to_csv(snakemake.output.transformed_patients, index=False)
df_diagnoses.to_csv(snakemake.output.transformed_diagnoses, index=False)