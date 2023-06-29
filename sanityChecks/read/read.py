"""
This file will contain functions required to read constraint-based model files in MatLab format
"""
import cobra
from cobra import Model
from pathlib import Path


def get_model(model_path: Path) -> Model:
    model: Model
    if model_path.suffix == ".json":
        model = cobra.io.load_json_model(filename=model_path)
    elif model_path.suffix == ".mat":
        model = cobra.io.load_matlab_model(infile_path=model_path)
    elif model_path.suffix in [".yml", ".yaml"]:
        model = cobra.io.load_yaml_model(filename=model_path)
    
    # For each reaction in the model, if it ends with "_c", "_m", or "_e", change it to "[c]", "[m]", or "[e]" respectively
    for reaction in model.reactions:
        if reaction.id.endswith("_c"):
            reaction.id = reaction.id.replace("_c", "[c]")
        elif reaction.id.endswith("_m"):
            reaction.id = reaction.id.replace("_m", "[m]")
        elif reaction.id.endswith("_e"):
            reaction.id = reaction.id.replace("_e", "[e]")
    
    if model.name == "":
        model.name = model_path.stem
    
    return model
