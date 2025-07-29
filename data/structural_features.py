import os
from multiprocessing import Pool
import pandas as pd
from Bio.PDB import PDBParser, DSSP


def extract_accessibilities(pdb_filepath: str) -> dict[int, float]:
    """Takes a PDB filepath and returns a mapping of site to surface accessibility."""

    name = os.path.basename(pdb_filepath).split(".")[0]
    temp_filepath = f"{name}_temp.pdb"
    wrote_header = False
    with open(temp_filepath, "w") as temp_pdb:
        temp_pdb.write("HEADER\n")
        with open(pdb_filepath) as original_pdb:
            for line in original_pdb:
                if not (line.startswith("HEADER") or line.startswith("CRYST1")):
                    if not wrote_header and (
                        line.startswith("ORIGX1") or line.startswith("ATOM")
                    ):
                        temp_pdb.write(
                            "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"
                        )
                        wrote_header = True
                    temp_pdb.write(line)
    assert wrote_header

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PDB_structure", temp_filepath)
    features = {res[0]: res[3] for res in DSSP(structure[0], temp_filepath)}
    os.remove(temp_filepath)

    return features


if __name__ == "__main__":
    # Extract structural features for all the PDBs in the structures directory

    pdb_dir = "structures/"
    pdb_fns = [fn for fn in os.listdir(pdb_dir) if fn[-4:] == ".pdb"]
    pdb_paths = [pdb_dir + pdb_fn for pdb_fn in pdb_fns]
    structural_features = Pool().map(extract_accessibilities, pdb_paths)
    structural_features_list = []
    for fn, protein_structural_features in zip(pdb_fns, structural_features):
        name = os.path.basename(fn).split(".")[0]
        for pos_features in protein_structural_features.items():
            structural_features_list.append([name, *pos_features])

    out_df = pd.DataFrame(structural_features_list)
    out_df.columns = ["name", "position", "surface_accessibility"]
    out_df.to_csv("structural_features.csv", index=False)
