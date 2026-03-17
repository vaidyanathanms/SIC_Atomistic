# Class to define the LAMMPS classes for combining data files
# Currently this works only for atomstyle - full
import re
from collections import defaultdict

class LammpsData:
    def __init__(self, filename):
        self.filename = filename

        self.header_counts = {
            "atoms": 0,
            "bonds": 0,
            "angles": 0,
            "dihedrals": 0,
            "impropers": 0,
            "atom types": 0,
            "bond types": 0,
            "angle types": 0,
            "dihedral types": 0,
            "improper types": 0,
        }

        self.masses = {}
        self.pair_coeffs = {}
        self.bond_coeffs = {}
        self.angle_coeffs = {}
        self.dihedral_coeffs = {}
        self.improper_coeffs = {}

        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.impropers = []

        self.box = {
            "xlo xhi": None,
            "ylo yhi": None,
            "zlo zhi": None,
            "xy xz yz": None,
        }

        self.read_data()

    def read_data(self):
        with open(self.filename, "r") as f:
            lines = f.readlines()

        section = None

        known_sections = {
            "Masses",
            "Pair Coeffs",
            "Bond Coeffs",
            "Angle Coeffs",
            "Dihedral Coeffs",
            "Improper Coeffs",
            "Atoms",
            "Bonds",
            "Angles",
            "Dihedrals",
            "Impropers",
        }

        i = 0
        while i < len(lines):
            raw = lines[i].rstrip("\n")
            raw_stripped = raw.strip()
            
            if raw_stripped == "":
                i += 1
                continue

            left, comment = self.split_comment(raw)
            stripped = left.strip()

            # section header like "Atoms # full"
            if stripped in known_sections:
                section = stripped
                i += 1
                continue

            # commented header like: # Bond Coeffs
            if raw_stripped.startswith("#"):
                comment_text = raw_stripped[1:].strip()
                if comment_text in known_sections:
                    section = comment_text
                    i += 1
                    continue

            # header counts
            m = re.match(
                r"^\s*(\d+)\s+(atoms|bonds|angles|dihedrals|impropers|atom types|bond types|angle types|dihedral types|improper types)\s*$",
                stripped
            )
            if m:
                self.header_counts[m.group(2)] = int(m.group(1))
                i += 1
                continue

            # box lines
            m = re.match(r"^\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+(xlo xhi|ylo yhi|zlo zhi)\s*$", stripped)
            if m:
                self.box[m.group(3)] = (float(m.group(1)), float(m.group(2)))
                i += 1
                continue

            m = re.match(r"^\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+(xy xz yz)\s*$", stripped)
            if m:
                self.box[m.group(4)] = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
                i += 1
                continue

            if section == "Masses":
                self.parse_coeff_line(raw, self.masses)

            elif section == "Pair Coeffs":
                self.parse_coeff_line(raw, self.pair_coeffs)

            elif section == "Bond Coeffs":
                self.parse_coeff_line(raw, self.bond_coeffs)

            elif section == "Angle Coeffs":
                self.parse_coeff_line(raw, self.angle_coeffs)

            elif section == "Dihedral Coeffs":
                self.parse_coeff_line(raw, self.dihedral_coeffs)

            elif section == "Improper Coeffs":
                self.parse_coeff_line(raw, self.improper_coeffs)

            elif section == "Atoms":
                self.parse_atom_line(raw)

            elif section == "Bonds":
                self.parse_bond_line(raw)

            elif section == "Angles":
                self.parse_angle_line(raw)

            elif section == "Dihedrals":
                self.parse_dihedral_line(raw)

            elif section == "Impropers":
                self.parse_improper_line(raw)

            i += 1


    def split_comment(self, line):
        if "#" in line:
            left, right = line.split("#", 1)
            return left.rstrip(), right.strip()
        return line.rstrip(), ""

            
    def parse_coeff_line(self, raw, store):
        left, comment = self.split_comment(raw)
        left = left.strip()

        if left:
            parts = left.split()
            try:
                type_id = int(parts[0])
                coeffs = parts[1:]
                label = comment
                store[type_id] = {"coeffs": coeffs, "label": label}
            except ValueError:
                pass
            return

        if comment:
            m = re.match(r"^\s*(\d+)\s+(.+)$", comment)
            if m:
                type_id = int(m.group(1))
                label = m.group(2).strip()
                store[type_id] = {"coeffs": [], "label": label}

    def parse_atom_line(self, raw):
        left, right = self.split_comment(raw)
        parts = left.split()

        if len(parts) < 7:
            return

        comment = ""
        if right.strip():
            comment = "#" + right.strip()
        
        self.atoms.append({
            "id": int(parts[0]),
            "mol": int(parts[1]),
            "type": int(parts[2]),
            "charge": float(parts[3]),
            "x": float(parts[4]),
            "y": float(parts[5]),
            "z": float(parts[6]),
            "extra": parts[7:],
            "comments": comment,
        })

    def parse_bond_line(self, raw):
        left, _ = self.split_comment(raw)
        parts = left.split()
        if len(parts) < 4:
            return

        self.bonds.append({
            "id": int(parts[0]),
            "type": int(parts[1]),
            "a1": int(parts[2]),
            "a2": int(parts[3]),
        })

    def parse_angle_line(self, raw):
        left, _ = self.split_comment(raw)
        parts = left.split()
        if len(parts) < 5:
            return

        self.angles.append({
            "id": int(parts[0]),
            "type": int(parts[1]),
            "a1": int(parts[2]),
            "a2": int(parts[3]),
            "a3": int(parts[4]),
        })

    def parse_dihedral_line(self, raw):
        left, _ = self.split_comment(raw)
        parts = left.split()
        if len(parts) < 6:
            return

        self.dihedrals.append({
            "id": int(parts[0]),
            "type": int(parts[1]),
            "a1": int(parts[2]),
            "a2": int(parts[3]),
            "a3": int(parts[4]),
            "a4": int(parts[5]),
        })

    def parse_improper_line(self, raw):
        left, _ = self.split_comment(raw)
        parts = left.split()
        if len(parts) < 6:
            return

        self.impropers.append({
            "id": int(parts[0]),
            "type": int(parts[1]),
            "a1": int(parts[2]),
            "a2": int(parts[3]),
            "a3": int(parts[4]),
            "a4": int(parts[5]),
        })

