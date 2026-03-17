# Auxiliary functions to support combining data files

from lmp_define import LammpsData
import re
import warnings

def combine_lammps_system(file_specs, spacing=5.0):
    """
    file_specs = [
        ("Li.data", 100),
        ("FSI.data", 100),
        ("PEO.data", 100),
        ("EC.data", 100),
    ]
    """

    data_objects = [LammpsData(fname) for fname, _ in file_specs]

    global_masses, atom_type_maps = build_global_type_map_dedup(data_objects, "masses")
    global_pair_coeffs, pair_type_maps = build_global_type_map_dedup(data_objects, "pair_coeffs")
    global_bond_coeffs, bond_type_maps = build_global_type_map_dedup(data_objects, "bond_coeffs")
    global_angle_coeffs, angle_type_maps = build_global_type_map_dedup(data_objects, "angle_coeffs")
    global_dihedral_coeffs, dihedral_type_maps = build_global_type_map_dedup(data_objects, "dihedral_coeffs")
    global_improper_coeffs, improper_type_maps = build_global_type_map_dedup(data_objects, "improper_coeffs")

    system = {
        "masses": global_masses,
        "pair_coeffs": global_pair_coeffs,
        "bond_coeffs": global_bond_coeffs,
        "angle_coeffs": global_angle_coeffs,
        "dihedral_coeffs": global_dihedral_coeffs,
        "improper_coeffs": global_improper_coeffs,
        "atoms": [],
        "bonds": [],
        "angles": [],
        "dihedrals": [],
        "impropers": [],
    }

    atom_id_offset = 0
    bond_id_offset = 0
    angle_id_offset = 0
    dihedral_id_offset = 0
    improper_id_offset = 0
    mol_id_offset = 0

    xshift = 0.0

    for ifile, ((fname, ncopy), data) in enumerate(zip(file_specs, data_objects)):
        atom_type_map = atom_type_maps[ifile]
        bond_type_map = bond_type_maps[ifile]
        angle_type_map = angle_type_maps[ifile]
        dihedral_type_map = dihedral_type_maps[ifile]
        improper_type_map = improper_type_maps[ifile]

        print(f'Processing {fname}')
        
        if data.atoms:
            xmin = min(a["x"] for a in data.atoms)
            xmax = max(a["x"] for a in data.atoms)
            dx = (xmax - xmin) + spacing
        else:
            dx = spacing

        local_mol_ids = sorted(set(a["mol"] for a in data.atoms))
        mol_reindex = {old: i + 1 for i, old in enumerate(local_mol_ids)}
        nmol_template = len(local_mol_ids)

        for icopy in range(ncopy):
            atom_id_map = {}

            for a in data.atoms:
                new_id = atom_id_offset + a["id"]
                atom_id_map[a["id"]] = new_id

                system["atoms"].append({
                    "id": new_id,
                    "mol": mol_id_offset + mol_reindex[a["mol"]],
                    "type": atom_type_map[a["type"]],
                    "charge": a["charge"],
                    "x": a["x"] + xshift,
                    "y": a["y"],
                    "z": a["z"],
                    "extra": a.get("extra",[]),
                    "comments": a.get("comments","")
                })

            for b in data.bonds:

                if b["type"] not in bond_type_map:
                    raise ValueError(
                        f"{fname}: bond type {b['type']} appears in Bonds section but was not found in Bond Coeffs. "
                        f"Parsed bond coeff types = {sorted(bond_type_map.keys())}"
                    )

                system["bonds"].append({
                    "id": bond_id_offset + b["id"],
                    "type": bond_type_map[b["type"]],
                    "a1": atom_id_map[b["a1"]],
                    "a2": atom_id_map[b["a2"]],
                })
                

            for ang in data.angles:
                system["angles"].append({
                    "id": angle_id_offset + ang["id"],
                    "type": angle_type_map[ang["type"]],
                    "a1": atom_id_map[ang["a1"]],
                    "a2": atom_id_map[ang["a2"]],
                    "a3": atom_id_map[ang["a3"]],
                })

            for dih in data.dihedrals:
                system["dihedrals"].append({
                    "id": dihedral_id_offset + dih["id"],
                    "type": dihedral_type_map[dih["type"]],
                    "a1": atom_id_map[dih["a1"]],
                    "a2": atom_id_map[dih["a2"]],
                    "a3": atom_id_map[dih["a3"]],
                    "a4": atom_id_map[dih["a4"]],
                })

            for imp in data.impropers:
                system["impropers"].append({
                    "id": improper_id_offset + imp["id"],
                    "type": improper_type_map[imp["type"]],
                    "a1": atom_id_map[imp["a1"]],
                    "a2": atom_id_map[imp["a2"]],
                    "a3": atom_id_map[imp["a3"]],
                    "a4": atom_id_map[imp["a4"]],
                })

            atom_id_offset += len(data.atoms)
            bond_id_offset += len(data.bonds)
            angle_id_offset += len(data.angles)
            dihedral_id_offset += len(data.dihedrals)
            improper_id_offset += len(data.impropers)
            mol_id_offset += nmol_template
            xshift += dx

    return system


def normalize_label(label):
    return label.strip().replace(" ", "").upper()


def build_global_type_map_dedup(data_objects, attr_name, require_coeff_match=True):
    """
    Deduplicate by label.
    Reuse same global type if label repeats.

    Example attr_name:
        masses
        pair_coeffs
        bond_coeffs
        angle_coeffs
        dihedral_coeffs
        improper_coeffs
    """

    global_defs = []
    local_to_global = []

    label_to_global = {}
    label_to_coeffs = {}

    next_id = 1

    for data in data_objects:
        local_dict = getattr(data, attr_name)
        this_map = {}

        for local_id in sorted(local_dict.keys()):
            entry = local_dict[local_id]

            raw_label = entry["label"]
            coeffs = tuple(entry["coeffs"])

            if raw_label.strip() == "":
                norm_label = f"__NO_LABEL__{attr_name}__{data.filename}__{local_id}"
                out_label = raw_label
            else:
                norm_label = normalize_label(raw_label)
                out_label = raw_label.strip()

            if norm_label in label_to_global:
                if require_coeff_match:
                    if label_to_coeffs[norm_label] != coeffs:
                        raise ValueError(
                            f"Mismatch for label '{raw_label}' in file {data.filename}. "
                            f"Old coeffs = {label_to_coeffs[norm_label]}, "
                            f"new coeffs = {coeffs}"
                        )

                this_map[local_id] = label_to_global[norm_label]

            else:
                gid = next_id
                next_id += 1

                label_to_global[norm_label] = gid
                label_to_coeffs[norm_label] = coeffs
                this_map[local_id] = gid

                global_defs.append({
                    "id": gid,
                    "coeffs": list(coeffs),
                    "label": out_label,
                })

        local_to_global.append(this_map)

    return global_defs, local_to_global

    
def write_coeff_section(f, title, entries):
    if not entries:
        return

    f.write(f"\n{title}\n\n")
    for entry in entries:
        if entry["coeffs"]:
            line = f"{entry['id']} " + " ".join(entry["coeffs"])
            if entry["label"]:
                line += f" # {entry['label']}"
        else:
            if entry["label"]:
                line = f"# {entry['id']} {entry['label']}"
            else:
                line = f"{entry['id']}"
        f.write(line.rstrip() + "\n")

    
def write_lammps_data(filename, system, box=None, write_atoms=None):

    if write_atoms is not None:
        with open(write_atoms,'w') as fat:

            fat.write("Just the atom data for verification\n\n")

            fat.write(f"{len(system['atoms'])} atoms\n")
            fat.write(f"{len(system['masses'])} atom types\n")

            fat.write("\nAtoms\n\n")

            for a in system["atoms"]:
                extra = ""
                comments = a.get("comments","")
                if a.get("extra"):
                    extra = " " + " ".join(a["extra"])
                if comments:
                    extra += " " + comments

                fat.write(
                    f"{a['id']} {a['mol']} {a['type']} {a['charge']:.8f} "
                    f"{a['x']:.8f} {a['y']:.8f} {a['z']:.8f}{extra}\n"
                )

    
    if box is None:
        if system["atoms"]:
            xs = [a["x"] for a in system["atoms"]]
            ys = [a["y"] for a in system["atoms"]]
            zs = [a["z"] for a in system["atoms"]]

            pad = 10.0
            box = (
                (min(xs) - pad, max(xs) + pad),
                (min(ys) - pad, max(ys) + pad),
                (min(zs) - pad, max(zs) + pad),
            )
        else:
            box = ((0.0, 100.0), (0.0, 100.0), (0.0, 100.0))

    with open(filename, "w") as f:
        f.write("Combined LAMMPS data file\n\n")

        f.write(f"{len(system['atoms'])} atoms\n")
        f.write(f"{len(system['bonds'])} bonds\n")
        f.write(f"{len(system['angles'])} angles\n")
        f.write(f"{len(system['dihedrals'])} dihedrals\n")
        f.write(f"{len(system['impropers'])} impropers\n\n")

        f.write(f"{len(system['masses'])} atom types\n")
        f.write(f"{len(system['bond_coeffs'])} bond types\n")
        f.write(f"{len(system['angle_coeffs'])} angle types\n")
        f.write(f"{len(system['dihedral_coeffs'])} dihedral types\n")
        f.write(f"{len(system['improper_coeffs'])} improper types\n\n")

        (xlo, xhi), (ylo, yhi), (zlo, zhi) = box
        f.write(f"{xlo:.6f} {xhi:.6f} xlo xhi\n")
        f.write(f"{ylo:.6f} {yhi:.6f} ylo yhi\n")
        f.write(f"{zlo:.6f} {zhi:.6f} zlo zhi\n")

        write_coeff_section(f, "Masses", system["masses"])
        write_coeff_section(f, "# Pair Coeffs", system["pair_coeffs"])
        write_coeff_section(f, "# Bond Coeffs", system["bond_coeffs"])
        write_coeff_section(f, "# Angle Coeffs", system["angle_coeffs"])
        write_coeff_section(f, "# Dihedral Coeffs", system["dihedral_coeffs"])
        write_coeff_section(f, "# Improper Coeffs", system["improper_coeffs"])

        f.write("\nAtoms\n\n")
        for a in system["atoms"]:

            extra = ""
            comments = a.get("comments","")
            if a.get("extra"):
                extra = " " + " ".join(a["extra"])
            if comments:
                extra += " " + comments

            f.write(
                f"{a['id']} {a['mol']} {a['type']} {a['charge']:.8f} "
                f"{a['x']:.8f} {a['y']:.8f} {a['z']:.8f}{extra}\n"
            )

        if system["bonds"]:
            f.write("\nBonds\n\n")
            for b in system["bonds"]:
                f.write(f"{b['id']} {b['type']} {b['a1']} {b['a2']}\n")

        if system["angles"]:
            f.write("\nAngles\n\n")
            for ang in system["angles"]:
                f.write(f"{ang['id']} {ang['type']} {ang['a1']} {ang['a2']} {ang['a3']}\n")

        if system["dihedrals"]:
            f.write("\nDihedrals\n\n")
            for dih in system["dihedrals"]:
                f.write(
                    f"{dih['id']} {dih['type']} {dih['a1']} {dih['a2']} "
                    f"{dih['a3']} {dih['a4']}\n"
                )

        if system["impropers"]:
            f.write("\nImpropers\n\n")
            for imp in system["impropers"]:
                f.write(
                    f"{imp['id']} {imp['type']} {imp['a1']} {imp['a2']} "
                    f"{imp['a3']} {imp['a4']}\n"
                )


