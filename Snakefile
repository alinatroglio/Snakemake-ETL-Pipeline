configfile: "config/config.yaml"

rule all:
    input:
        "output/ok_files/load_to_database.ok"

rule load_to_database:
    input:
        merged_data = "output/merged/merged_data.csv"
    output:
        ok = "output/ok_files/load_to_database.ok"
    script:
        "scripts/load_to_database.py"

rule merge_data:
    input:
        transformed_encounters = "output/transformed/transformed_encounters.csv",
        transformed_patients = "output/transformed/transformed_patients.csv",
        transformed_diagnoses= "output/transformed/transformed_diagnoses.csv"
    output:
        merged_data = "output/merged/merged_data.csv"
    script:
        "scripts/merge_data.py"

rule load_data:
    input:
        encounters = "data/encounters.csv",
        patients= "data/patients.csv",
        diagnoses="data/diagnoses.xml"
    output:
        transformed = "output/transformed/transformed_encounters.csv",
        transformed_patients = "output/transformed/transformed_patients.csv",
        transformed_diagnoses= "output/transformed/transformed_diagnoses.csv"
    script:
        "scripts/load_data.py"