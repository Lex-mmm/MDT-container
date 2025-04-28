import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import json

def flatten_json(d, parent_key='', sep='.'):
    """
    Recursively flattens a nested JSON/dictionary into a single level dict.
    
    Args:
        d (dict): The nested dictionary.
        parent_key (str): The base key string (for recursion).
        sep (str): Separator between nested keys.
    
    Returns:
        dict: A flattened dictionary.
    """
    items = {}
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.update(flatten_json(v, new_key, sep=sep))
        elif isinstance(v, list):
            for idx, item in enumerate(v):
                list_key = f"{new_key}[{idx}]"
                if isinstance(item, dict):
                    items.update(flatten_json(item, list_key, sep=sep))
                else:
                    items[list_key] = item
        else:
            items[new_key] = v
    return items
with open("healthyFlat.json", "r") as f:
    healthy = json.load(f)

print(healthy["params.M_CO2"]["value"])