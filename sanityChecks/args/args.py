import argparse
from cobra import Model
from pathlib import Path
from dataclasses import dataclass
import read


@dataclass
class Arguments:
    model_path: Path
    solver: str
    is_human: bool
    __model: Model = Model()
    
    def __post_init__(self):
        self.model_path = Path(self.model_path)
        if not self.model_path.exists():
            raise FileNotFoundError(f"Model file {self.model_path} not found")
        
        self.__model = read.get_model(self.model_path)
    
    @property
    def model(self) -> Model:
        return self.__model


def parse_args() -> Arguments:
    parser = argparse.ArgumentParser(description="Sanity Check for Model")
    parser.add_argument(
        "-m", "--model",
        type=str,
        dest="model_path",
        required=True,
        help="Path to model",
    )
    parser.add_argument(
        "-s", "--solver",
        type=str,
        dest="solver",
        required=True,
        help="Solver to use",
        choices=["gurobi", "glpk"],
    )
    
    # Create a new argument group called "organism" with two options: "human" and "non-human"
    # Only one if them is allowed, but one is required
    organism = parser.add_mutually_exclusive_group(required=True)
    organism.add_argument(
        "-H", "--human",
        action="store_true",
        dest="is_human",
        required=False,
        help="Model is human",
    )
    organism.add_argument(
        "-N", "--non-human",
        action="store_false",
        dest="is_human",
        required=False,
        help="Model is NOT human",
    )
    
    args = Arguments(**vars(parser.parse_args()))
    return args
