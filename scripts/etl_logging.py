import logging
from datetime import datetime
import json, os

LOG_FILE = "/app/shared/logs.json"

def reset_loggin():
    """
    Reset the shared log file at workflow start. 
    """
    os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
    with open(LOG_FILE, "w") as f:
        json.dump([], f)

def log_inconsistency(row_id, file_name, issue_type, details):
    """
    Append a structured data inconsistency record to the log file.  

    Parameters
    ----------
    row_id 
        Row ID that contains inconsistency.
    file_name 
        Name of source file where inconsistency occured.
    issue_type
        Short label describing the type of inconsistency
    details
        Additional context
    """
    log_element = {"timestamp": datetime.now().isoformat(), 
                   "row_id": row_id,
                   "file_name": file_name,
                   "issue_type": issue_type,
                   "details": details}
    log_collection = []
    if os.path.exists(LOG_FILE):
        with open(LOG_FILE, "r") as f:
            try:
                log_collection = json.load(f)
            except:
                pass
    log_collection.append(log_element)
    with open(LOG_FILE, "w") as f:
        json.dump(log_collection, f)

def get_log_collection():
    """
    Load and return all logged inconsistency records.

    Returns
    -------
        List of logged inconsistency entries, or an empty list if none exist.
    """
    if os.path.exists(LOG_FILE):
        with open(LOG_FILE, "r") as f:
            return json.load(f)
    return []