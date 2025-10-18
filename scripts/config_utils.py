import yaml

def load_yaml(path):
    with open(path, "r") as f_yaml:
        return yaml.safe_load(f_yaml)
        
def load_mapping(path):
    yaml_data = load_yaml(path)
    return yaml_data

def map_values(df, col,col_mapped, map):
    df[col_mapped] = df[col].str.upper().map(map)
    return df