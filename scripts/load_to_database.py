from snakemake.script import snakemake
from sqlalchemy import create_engine, select, Table, MetaData
import psycopg2
import pandas as pd
import streamlit as st
from pathlib import Path

# Get DB URL and map ID types to string
db_url = snakemake.config["env"]["db_url"]
id_type = {"encounter_id": "string", "patient_id": "string", "row_id": "string"}

# Load merged data
df_merged_data = pd.read_csv(snakemake.input.merged_data, dtype=id_type,na_values="")

# Load data to SQL database
if db_url:
    engine = create_engine(db_url)
    df_merged_data.to_sql("patient_data", engine, if_exists="replace", index=False)

# Create OK file for Snakemake
ok_path = Path(snakemake.output.ok)
ok_path.parent.mkdir(parents=True, exist_ok=True)
ok_path.touch()
