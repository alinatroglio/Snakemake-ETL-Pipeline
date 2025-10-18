import csv
import re
import pandas as pd
import pint 
from dateutil import parser
from etl_logging import log_inconsistency

# Initialize Pint unit registry for handling physical units and conversions
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

def clean_dup_header_and_spaces(content, file_name):
    """
    Detect and remove duplicate headers and extra spaces from a file's content.
    
    Parameters
    ----------
    content 
        Lines of file content.
    file_name 
        Name of source file for logging. 
    Returns
    -------
    Cleaned list of lines with duplicates and extra spaces removed.
    """
    header = content[0]
    indices = [i+1 for i, s in enumerate(content[1:]) if header in s]
    if indices:
        for idx in indices:
            log_inconsistency(idx, file_name, "Duplicate Header", f"Duplicated header in row {idx}")
    return [item.strip() for i, item in enumerate(content) if i not in indices]
     
def clean_delimiter_errors(cleaned_content, file_name):
    """
    Normalize delimiter in csv-like text and clean header. 

    Notes
    -----
    - Rows with mismatched column counts are logged.
    - Extra columns are truncated to the header length.

    Parameters
    ----------
    cleanded_content 
        Cleaned lines of file content.
    file_name 
        Name of source file for logging. 
    Returns
    -------
    Cleaned list of lines with replaced delimiter, detected delimtier, 
    and header with extra spaces removed.
    """
    if not cleaned_content:
        return cleaned_content, ",", []
    sample = "".join(cleaned_content[:2048]) if len(cleaned_content) > 0 else ""
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=[',',';','\t'])
        delimiter = dialect.delimiter
    except Exception:
        delimiter = ","
        log_inconsistency(0, file_name, "Delimiter sniff failed", "Use , as default")
    header = cleaned_content[0].split(delimiter)
    header_cleaned = [h.strip() for h in header]
    column_num = len(header)

    for i, row in enumerate(cleaned_content):
        if not len(row.split(delimiter)) == column_num:
            log_inconsistency(i, file_name, "Delimiter Error", f"Row has not same column count with {delimiter}")
            try:
                dialect = csv.Sniffer().sniff(row, delimiters=[',',';','\t'])
                delimiter_new = dialect.delimiter
                cleaned_content[i] = row.replace(delimiter_new, delimiter)
                log_inconsistency(i, file_name, "Delimiter Fixed", "Replaced sniffed delimiter")
                if len(cleaned_content[i].split(delimiter)) > column_num:
                    split_row = cleaned_content[i].split(delimiter)
                    cleaned_content[i] = delimiter.join(split_row[:column_num])
                    log_inconsistency(i, file_name, "Remove extra column", "Row had extra column, mismatched number was removed")
            except Exception:
                continue
    return cleaned_content, delimiter, header_cleaned

def clean_csv_to_df(input_file):
    """
    Read a csv file, clean formatting issues and load clean data as DataFrame.
    
    Parameters
    ----------
    input_file 
        Path to raw input file.
    Returns
    -------
    DataFrame containing the cleade data.
    """
    try:
        with open(input_file, encoding='utf-8-sig') as f_in:
            f_in_content = [line for line in f_in.readlines() if line.strip()]
            cleaned_content = clean_dup_header_and_spaces(f_in_content, input_file)
            cleaned_content, delimiter, header = clean_delimiter_errors(cleaned_content, input_file)
            rows = [row.split(delimiter) for row in cleaned_content[1:]]
            df = pd.DataFrame(rows, columns=header)
            return df
    except Exception as e:
        log_inconsistency(0, input_file, "Could not read file", f"Error while reading file {input_file} with {e}")

def clean_str_data(df, str_safe):
    """
    Strip whitspaces of selected columns, where it is safe.

    Parameters
    ----------
    df 
        Input DataFrame.
    str_safe
        Columns safe to clean.
    Returns
    -------
    DataFrame with trimmed string values.
    """
    for col in df.columns.intersection(str_safe):
        df[col] = df[col].astype(str).str.strip()
    return df

def capitalize(df, cap_safe):
    """
    Capitalize selected columns, that are name-like.

    Parameters
    ----------
    df 
        Input DataFrame.
    str_safe
        Columns safe to cap.
    Returns
    -------
    DataFrame with capitalized values.
    """
    for col in df.columns.intersection(cap_safe):
        df[col] = df[col].astype(str).str.lower().str.title()
    return df

def fix_ids(df, id_col, file_name):
    """
    Remove duplicated and missing IDs.

    Parameters
    ----------
    df 
        Input DataFrame.
    id_col
        Column name of ID.
    file_name
        File name for logging.
    Returns
    -------
    DataFrame with cleaned IDs only.
    """
    dup_mask = df[id_col].duplicated(keep="first")
    for idx in df[dup_mask].index:
        log_inconsistency(idx, file_name, "Duplicate ID", f"Duplicate {id_col}: row {idx}, dupliactes are removed")
    nan_mask = df[id_col].isna() | (df[id_col].astype(str).str.strip() == "")
    for idx in df[nan_mask].index:
        log_inconsistency(idx, file_name, "Missing ID", f"Missing {id_col} in row {idx}")
    df = df[~dup_mask & ~nan_mask]
    return df


def parse_date(row, row_id, file_name):
    """
    Parse date for unification. 

    Parameters
    ----------
    row 
        Date string to parse. 
    row_id
        Row ID for logging.
    file_name
        File name for logging.
    Returns
    -------
    Parsed datetime or NaT if invalid.
    """
    try:
        return parser.parse(row)
    except Exception:
        log_inconsistency(row_id, file_name, "Invalid date", f"Not parseable {row}")
        return pd.NaT
    
def unify_datetime(df, dates, file_name):
    """
    Normalize datetime or NaT if invalid.

    Parameters
    ----------
    df 
        Input DataFrame.
    dates
        Date columns to unify.
    file_name
        File name for logging.
    Returns
    -------
    DataFrame with standardized datetime columns.    
    """
    for col in df.columns.intersection(dates):
            df[col] = df[col].str.replace('-', '/')
            df[col] = df.apply(lambda row: parse_date(row[col], row.name, file_name), axis=1)
    return df


def unify_height(value, row_id, file_name):
    """
    Normalize height to cm or None if invalid.

    Parameters
    ----------
    value 
        Height to normalize.
    row_id
        Row ID for logging.
    file_name
        File name for logging.
    Returns
    -------
    Height in centimeters, or None if invalid.
    """
    # lower case text
    value = str(value).lower()
    numbers = re.findall(r"\d+(?:\.\d+)?", value)
    if not numbers:
        log_inconsistency(row_id, file_name, "Invalid heights", "No numeric value, return None")
        return None
    # case cm:
    if re.search(r"cm\b", value):
        return round(float(numbers[0]), 1)
    # case meter:
    if re.search(r"\bm\b", value):
        return round(Q_(float(numbers[0]),"meter").to("cm").magnitude, 1)
        # case feet and inch:
    if re.search(r"ft\b", value):
        if re.search(r"in\b",value):
            print(value)
            return round((Q_(float(numbers[0]), "feet") + Q_(float(numbers[1]), "inch")).to("cm").magnitude, 1)
    # case inch:
    if re.search(r"in\b", value):
        return round(Q_(float(numbers[0]),"inch").to("cm").magnitude, 1)
    # default cm if no unit
    try:
        return round(float(numbers[0]), 1)
    # otherwise return None
    except Exception as e: 
        log_inconsistency(row_id, file_name, "Height cannot be unified", f"Error {e}")
        return None
    
def unify_weight(value, row_id, file_name):
    """
    Normalize weight to kg or None if invalid.

    Parameters
    ----------
    value 
        Weight to normalize.
    row_id
        Row ID for logging.
    file_name
        File name for logging.
    Returns
    -------
    Weight in kg, or None if invalid.
    """
    # lower case text
    value = str(value).lower()
    numbers = re.findall(r"\d+(?:\.\d+)?", value)
    if not numbers:
        log_inconsistency(row_id, file_name, "Invalid weight", "No numeric value, return None")
        return None
    # case kg:
    if re.search(r"kg\b", value):
        return round(float(numbers[0]), 1)
    # case meter:
    if re.search(r"lbs?\b", value):
        return round(Q_(float(numbers[0]),"lb").to("kg").magnitude, 1)
    # default kg if no unit
    try:
        return round(float(numbers[0]), 1)
    # otherwise return None
    except Exception as e: 
        log_inconsistency(row_id, file_name, "Height cannot be unified", f"Error {e}")
        return None

def validate_icd10_code(value, row_id, file_name):
    """
    Validate ICD-10 code format.

    Parameters
    ----------
    value
        ICD-10 code to validate.
    row_id
        Row index for logging.
    file_name
        File name for logging.
    Returns
    -------
    True if valid ICD-10 code, False otherwise.
    """
    if pd.isna(value) or str(value).strip() == "":
        log_inconsistency(row_id, file_name, "Empty ICD-10 code", "Empty or nan")
        return False
    match = re.search(r"^[A-Za-z][0-9]{2}(\.[A-Z0-9]{1,2})?", value)
    if not match:
        log_inconsistency(row_id, file_name, "Invalid ICD-10 code", f"Code {value} is not in the right format")
        return False
    return True 