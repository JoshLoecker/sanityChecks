import cobra
from cobra import Model
from pathlib import Path

# Custom functions
from read import read
from args.args import parse_args
from logger.logger import get_logger

# Sanity Checks
from checks import atp_yield

get_logger()


def matlab(model_path: Path) -> Model:
    return cobra.io.load_model(str(model_path))


def main():
    args_ = parse_args()
    model: Model = read.get_model(args_.model_path)
    atp_yield.ATPYield(arguments=args_)


if __name__ == "__main__":
    import sys
    
    sys.argv.extend([
        "--model", "/Users/joshl/PycharmProjects/sanityChecks/data/naiveB_SpecificModel_IMAT.mat",
        "--solver", "gurobi",
        "--non-human"
    
    ])
    main()
