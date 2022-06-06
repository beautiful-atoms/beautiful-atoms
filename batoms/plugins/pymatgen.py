"""
Search public database

Pymatgen
Pubchem
Nomad
"""

from batoms import Batoms
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


def pymatgen_search(key, id):
    from pymatgen.ext.matproj import MPRester
    with MPRester(key) as m:
        # Structure for material id
        structure = m.get_structure_by_material_id(id)
        label = id.replace('-', '_')
        batoms = Batoms(label, from_pymatgen=structure)
        # Dos for material id
        return batoms
