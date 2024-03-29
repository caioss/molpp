#include <molpp/ElementsTable.hpp>
#include <algorithm>

using namespace mol;

ElementsTable::ElementsTable(std::initializer_list<Element> data)
: m_max_vdw{0},
  m_max_covalent{0}
{
    m_atomic.reserve(data.size());
    m_mass.reserve(data.size());
    m_covalent.reserve(data.size());
    m_vdw.reserve(data.size());
    m_symbol.reserve(data.size());
    m_name.reserve(data.size());

    for (Element const &element : data)
    {
        m_atomic.push_back(element.atomic_number);
        m_mass.push_back(element.mass);
        m_covalent.push_back(element.covalent_radius);
        m_vdw.push_back(element.VDW_radius);
        m_symbol.push_back(element.symbol);
        m_name.push_back(element.name);
        m_max_covalent = std::max(m_max_covalent, element.covalent_radius);
        m_max_vdw = std::max(m_max_vdw, element.VDW_radius);
    }
}

// The following lines are generated automatically.
// Modifications done here are not guarenteed to be preserved.
ElementsTable const& mol::ELEMENTS_TABLE()
{
    static ElementsTable table{
        {0, 0.0000, 0.0, 0, "Xx", "Dummy"},
        {1, 1.008, 0.37, 1.2, "H", "Hydrogen"},
        {2, 4.002602, 0.32, 1.4, "He", "Helium"},
        {3, 6.94, 1.34, 2.2, "Li", "Lithium"},
        {4, 9.012182, 0.90, 1.9, "Be", "Beryllium"},
        {5, 10.81, 0.82, 1.8, "B", "Boron"},
        {6, 12.011, 0.77, 1.7, "C", "Carbon"},
        {7, 14.007, 0.75, 1.6, "N", "Nitrogen"},
        {8, 15.999, 0.73, 1.55, "O", "Oxygen"},
        {9, 18.9984032, 0.71, 1.5, "F", "Fluorine"},
        {10, 20.1797, 0.69, 1.54, "Ne", "Neon"},
        {11, 22.98976928, 1.54, 2.4, "Na", "Sodium"},
        {12, 24.305, 1.30, 2.2, "Mg", "Magnesium"},
        {13, 26.9815386, 1.18, 2.1, "Al", "Aluminium"},
        {14, 28.085, 1.11, 2.1, "Si", "Silicon"},
        {15, 30.973762, 1.06, 1.95, "P", "Phosphorus"},
        {16, 32.06, 1.02, 1.8, "S", "Sulfur"},
        {17, 35.45, 0.99, 1.8, "Cl", "Chlorine"},
        {18, 39.948, 0.97, 1.88, "Ar", "Argon"},
        {19, 39.0983, 1.96, 2.8, "K", "Potassium"},
        {20, 40.078, 1.74, 2.4, "Ca", "Calcium"},
        {21, 44.955912, 1.44, 2.3, "Sc", "Scandium"},
        {22, 47.867, 1.36, 2.15, "Ti", "Titanium"},
        {23, 50.9415, 1.25, 2.05, "V", "Vanadium"},
        {24, 51.9961, 1.27, 2.05, "Cr", "Chromium"},
        {25, 54.938045, 1.39, 2.05, "Mn", "Manganese"},
        {26, 55.845, 1.25, 2.05, "Fe", "Iron"},
        {27, 58.933195, 1.26, 2, "Co", "Cobalt"},
        {28, 58.6934, 1.21, 2, "Ni", "Nickel"},
        {29, 63.546, 1.38, 2, "Cu", "Copper"},
        {30, 65.38, 1.31, 2.1, "Zn", "Zinc"},
        {31, 69.723, 1.26, 2.1, "Ga", "Gallium"},
        {32, 72.630, 1.22, 2.1, "Ge", "Germanium"},
        {33, 74.92160, 1.19, 2.05, "As", "Arsenic"},
        {34, 78.96, 1.16, 1.9, "Se", "Selenium"},
        {35, 79.904, 1.14, 1.9, "Br", "Bromine"},
        {36, 83.798, 1.10, 2.02, "Kr", "Krypton"},
        {37, 85.4678, 2.11, 2.9, "Rb", "Rubidium"},
        {38, 87.62, 1.92, 2.55, "Sr", "Strontium"},
        {39, 88.90585, 1.62, 2.4, "Y", "Yttrium"},
        {40, 91.224, 1.48, 2.3, "Zr", "Zirconium"},
        {41, 92.90638, 1.37, 2.15, "Nb", "Niobium"},
        {42, 95.96, 1.45, 2.1, "Mo", "Molybdenum"},
        {43, 97, 1.56, 2.05, "Tc", "Technetium"},
        {44, 101.07, 1.26, 2.05, "Ru", "Ruthenium"},
        {45, 102.90550, 1.35, 2, "Rh", "Rhodium"},
        {46, 106.42, 1.31, 2.05, "Pd", "Palladium"},
        {47, 107.8682, 1.53, 2.1, "Ag", "Silver"},
        {48, 112.411, 1.48, 2.2, "Cd", "Cadmium"},
        {49, 114.818, 1.44, 2.2, "In", "Indium"},
        {50, 118.710, 1.41, 2.25, "Sn", "Tin"},
        {51, 121.760, 1.38, 2.2, "Sb", "Antimony"},
        {52, 127.60, 1.35, 2.1, "Te", "Tellurium"},
        {53, 126.90447, 1.33, 2.1, "I", "Iodine"},
        {54, 131.293, 1.30, 2.16, "Xe", "Xenon"},
        {55, 132.9054519, 2.25, 3, "Cs", "Caesium"},
        {56, 137.327, 1.98, 2.7, "Ba", "Barium"},
        {57, 138.90547, 1.69, 2.5, "La", "Lanthanum"},
        {58, 140.116, 0.0, 2.48, "Ce", "Cerium"},
        {59, 140.90765, 0.0, 2.47, "Pr", "Praseodymium"},
        {60, 144.242, 0.0, 2.45, "Nd", "Neodymium"},
        {61, 145, 0.0, 2.43, "Pm", "Promethium"},
        {62, 150.36, 0.0, 2.42, "Sm", "Samarium"},
        {63, 151.964, 0.0, 2.4, "Eu", "Europium"},
        {64, 157.25, 0.0, 2.38, "Gd", "Gadolinium"},
        {65, 158.92535, 0.0, 2.37, "Tb", "Terbium"},
        {66, 162.500, 0.0, 2.35, "Dy", "Dysprosium"},
        {67, 164.93032, 0.0, 2.33, "Ho", "Holmium"},
        {68, 167.259, 0.0, 2.32, "Er", "Erbium"},
        {69, 168.93421, 0.0, 2.3, "Tm", "Thulium"},
        {70, 173.054, 0.0, 2.28, "Yb", "Ytterbium"},
        {71, 174.9668, 1.60, 2.27, "Lu", "Lutetium"},
        {72, 178.49, 1.50, 2.25, "Hf", "Hafnium"},
        {73, 180.94788, 1.38, 2.2, "Ta", "Tantalum"},
        {74, 183.84, 1.46, 2.1, "W", "Tungsten"},
        {75, 186.207, 1.59, 2.05, "Re", "Rhenium"},
        {76, 190.23, 1.28, 2, "Os", "Osmium"},
        {77, 192.217, 1.37, 2, "Ir", "Iridium"},
        {78, 195.084, 1.28, 2.05, "Pt", "Platinum"},
        {79, 196.966569, 1.44, 2.1, "Au", "Gold"},
        {80, 200.592, 1.49, 2.05, "Hg", "Mercury"},
        {81, 204.38, 1.48, 2.2, "Tl", "Thallium"},
        {82, 207.2, 1.47, 2.3, "Pb", "Lead"},
        {83, 208.98040, 1.46, 2.3, "Bi", "Bismuth"},
        {84, 209, 0.0, 2, "Po", "Polonium"},
        {85, 210, 0.0, 2, "At", "Astatine"},
        {86, 222, 1.45, 2, "Rn", "Radon"},
        {87, 223, 0.0, 2, "Fr", "Francium"},
        {88, 226, 0.0, 2, "Ra", "Radium"},
        {89, 227, 0.0, 2, "Ac", "Actinium"},
        {90, 232.03806, 0.0, 2.4, "Th", "Thorium"},
        {91, 231.03588, 0.0, 2, "Pa", "Protactinium"},
        {92, 238.02891, 0.0, 2.3, "U", "Uranium"},
        {93, 237, 0.0, 2, "Np", "Neptunium"},
        {94, 244, 0.0, 2, "Pu", "Plutonium"},
        {95, 243, 0.0, 2, "Am", "Americium"},
        {96, 247, 0.0, 2, "Cm", "Curium"},
        {97, 247, 0.0, 2, "Bk", "Berkelium"},
        {98, 251, 0.0, 2, "Cf", "Californium"},
        {99, 252, 0.0, 2, "Es", "Einsteinium"},
        {100, 257, 0.0, 2, "Fm", "Fermium"},
        {101, 258, 0.0, 2, "Md", "Mendelevium"},
        {102, 259, 0.0, 2, "No", "Nobelium"},
        {103, 262, 0.0, 2, "Lr", "Lawrencium"},
        {104, 267, 0.0, 2, "Rf", "Rutherfordium"},
        {105, 270, 0.0, 2, "Db", "Dubnium"},
        {106, 271, 0.0, 2, "Sg", "Seaborgium"},
        {107, 270, 0.0, 2, "Bh", "Bohrium"},
        {108, 277, 0.0, 2, "Hs", "Hassium"},
        {109, 276, 0.0, 2, "Mt", "Meitnerium"},
        {110, 281, 0.0, 0.0, "Ds", "Darmstadtium"},
        {111, 282, 0.0, 0.0, "Rg", "Roentgenium"},
        {112, 285, 0.0, 0.0, "Cn", "Copernicium"},
        {113, 285, 0.0, 0.0, "Uut", "Ununtrium"},
        {114, 289, 0.0, 0.0, "Fl", "Flerovium"},
        {115, 289, 0.0, 0.0, "Uup", "Ununpentium"},
        {116, 293, 0.0, 0.0, "Lv", "Livermorium"},
        {117, 294, 0.0, 0.0, "Uus", "Ununseptium"},
        {118, 294, 0.0, 0.0, "Uuo", "Ununoctium"},
    };

    return table;
}
