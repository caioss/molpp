import xml.etree.ElementTree as ET
import requests


COUNTS_URL = "http://ligand-expo.rcsb.org/dictionaries/cc-counts.tdd"
RESIDUES_URL = "http://ligand-expo.rcsb.org/reports/{0}/{1}/{1}.xml"
MIN_COUNT = 500


PDB_ORDERS = {"sing": 1, "doub": 2, "trip": 3, "quad": 4}

response = requests.get(COUNTS_URL)
assert response.status_code == requests.status_codes.codes.OK

print("extern ResiduesTable const mol::internal::RESIDUES_TABLE({")
for line in response.text.splitlines()[1:]:
    resname, count = line.split()
    if int(count) < MIN_COUNT:
        break

    if resname == "UNX":
        continue

    response = requests.get(RESIDUES_URL.format(resname[0], resname))
    assert response.status_code == requests.status_codes.codes.OK
    xml_root = ET.fromstring(response.text)

    name = xml_root.findall("{*}chem_compCategory/{*}chem_comp/{*}name")[0].text

    atoms = []
    atom_ids = {}
    for i, atom_tag in enumerate(xml_root.iterfind("{*}chem_comp_atomCategory/*")):
        atom_ids[atom_tag.attrib["atom_id"]] = i
        atoms.append('{{"{}", {:d}}}'.format(atom_tag.attrib["atom_id"], i))

    bonds = []
    for bond_tag in xml_root.iterfind("{*}chem_comp_bondCategory/*"):
        atom1 = atom_ids[bond_tag.attrib["atom_id_1"]]
        atom2 = atom_ids[bond_tag.attrib["atom_id_2"]]

        bond_order = PDB_ORDERS.get(bond_tag.find("{*}value_order").text, 1)

        if bond_tag.find("{*}pdbx_aromatic_flag").text == "Y":
            is_aromatic = "true"
        else:
            is_aromatic = "false"

        bonds.append('{{{}, {:d}, {:d}, {:d}}}'.format(is_aromatic, bond_order, atom1, atom2))

    print("    // {}".format(name))
    print('    {{"{}",'.format(resname))
    print("         // Atoms")
    print("        {{", end="")
    print(",\n          ".join(atoms), end="")
    print("},")
    print("         // Bonds")
    print("         {", end="")
    print(",\n          ".join(bonds), end="")
    print("}}},")

print("});")
