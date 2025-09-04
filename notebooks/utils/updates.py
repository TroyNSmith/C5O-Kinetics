import yaml

from pathlib import Path


def update_yaml(click):
    summary = {}
    calc_dir = Path.home() / "C5O-Kinetics/calc"
    for smiles in filter(lambda p: p.is_dir(), calc_dir.iterdir()):
        smile_dict = {}
        for procedure in filter(lambda p: p.is_dir(), smiles.iterdir()):
            procedure_dict = {}
            run_dir = procedure / "run"
            if run_dir.exists():
                energies = {"SPC": [], "ZPV": []}
                method_dict = {}
                for method in filter(lambda p: p.is_dir(), run_dir.iterdir()):
                    xyzs = [str(x) for x in method.rglob("*xyz")]
                    log_files = list(method.rglob("*.log"))
                    for log_file in log_files:
                        if log_file.name == "REVDSD.log":
                            with open(log_file) as lf:
                                energies["ZPV"] += [
                                    line.strip().split()[-2]
                                    for line in lf
                                    if "Zero point energy" in line
                                    and "kcal/mol" in line
                                ]
                        elif log_file.name == "CCSDT.log":
                            with open(log_file) as lf:
                                energies["SPC"] += [
                                    str(float(line.strip().split()[-1]) * 627.509)
                                    for line in lf
                                    if "FINAL SINGLE POINT ENERGY" in line
                                ]
                    if xyzs:
                        method_dict[method.name] = {"xyzs": xyzs}
                if method_dict:
                    procedure_dict["methods"] = method_dict
                if energies["SPC"] or energies["ZPV"]:
                    procedure_dict["energies"] = energies
            if procedure_dict:
                smile_dict[procedure.name] = procedure_dict
        if smile_dict:
            summary[smiles.name] = smile_dict

    with open("results.yaml", "w") as f:
        yaml.safe_dump(summary, f)
