import os

from dotenv import load_dotenv
from mp_api.client import MPRester
from pymatgen.core import Structure

import moyopy


def main():
    material_id = "mp-550745"
    print(f"{material_id=}")
    with MPRester(api_key=os.environ.get("MP_API_KEY")) as mpr:
        doc = mpr.materials.summary.search(material_ids=[material_id])[0]
        structure: Structure = doc.structure

        basis = structure.lattice.matrix
        positions = structure.frac_coords
        numbers = [site.specie.Z for site in structure]

        moyopy_cell = moyopy.Cell(basis.tolist(), positions.tolist(), numbers)
        with open(f"{material_id}.json", "w") as f:
            f.write(moyopy_cell.serialize_json())

        for symprec in [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1]:
            moyopy_dataset = moyopy.MoyoDataset(moyopy_cell, symprec=symprec)
            print(symprec, moyopy_dataset.number)


if __name__ == "__main__":
    assert load_dotenv()
    main()
