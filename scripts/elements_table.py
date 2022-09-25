import xml.etree.ElementTree as ET
import requests


ELEMENTS_URL = "https://svn.code.sf.net/p/bodr/code/trunk/bodr/elements/elements.xml"


response = requests.get(ELEMENTS_URL)
assert response.status_code == requests.status_codes.codes.OK
xml_root = ET.fromstring(response.text)

print("ElementsTable const& mol::ELEMENTS_TABLE()\n{\n    static ElementsTable table{")

for atom_tag in xml_root.iterfind("{*}atom"):
    vdw_tags = atom_tag.findall(".//*[@dictRef='bo:radiusVDW']")
    vdw_radius = vdw_tags[0].text if vdw_tags else "0.0"

    covalent_tags = atom_tag.findall(".//*[@dictRef='bo:radiusCovalent']")
    covalent_radius = covalent_tags[0].text if covalent_tags else "0.0"

    data = [
        atom_tag.findall(".//*[@dictRef='bo:atomicNumber']")[0].text,
        atom_tag.findall(".//*[@dictRef='bo:mass']")[0].text,
        covalent_radius,
        vdw_radius,
        atom_tag.findall(".//*[@dictRef='bo:symbol']")[0].attrib["value"],
        atom_tag.findall(".//*[@dictRef='bo:name']")[0].attrib["value"],
    ]

    print('        {{{}, {}, {}, {}, "{}", "{}"}},'.format(*data))

print("    };\n\n    return table;\n}")
