import streamlit as st
import pandas as pd
from etl_logging import get_log_collection
from pathlib import Path
import yaml 


# Path to your config file (adjust if needed)
config_path = Path("/app/config/config.yaml")

# Load config setting
with open(config_path, "r") as f:
    config = yaml.safe_load(f)

# Establish SQL connection 
db_url = config["env"]["db_url"]
conn_st = st.connection(name="postgres", type='sql', url = db_url)

# Select data from database
patient_data = conn_st.query("SELECT * FROM patient_data")

# Dashboad setup
st.title("Patient Data")

# Add dataframe overview
st.header("Data Preview")
st.dataframe(patient_data, width= "stretch")

# Display missing values per column
st.header("Data Quality")
nans = patient_data.isna().sum().to_frame("missing_values")
st.subheader("Missing values per column")
st.dataframe(nans)

# Load and display inconsistencies
logs = get_log_collection()

if logs: 
    df_logs = pd.DataFrame(logs)
    st.subheader("Data Quality Issues")
    st.dataframe(df_logs)
else:
    st.info("No data inconsistency logged yet.")

# Show data distributions
st.header("Data Distributions")
numeric_cols = patient_data.select_dtypes(include="number").columns.tolist()
str_cols = patient_data.select_dtypes(exclude="number").columns.tolist()
if numeric_cols:
    col = st.selectbox("Numeric column", numeric_cols) 
    if col:
        st.bar_chart(patient_data[col].dropna())
if str_cols:
    col = st.selectbox("Non-numeric column", str_cols)
    if col:
        st.bar_chart(patient_data[col].value_counts(dropna=False))
