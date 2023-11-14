"""Search  the Materials Project"""

from batoms import Batoms
import logging

logger = logging.getLogger(__name__)


def mp_search(key, mpid):
    from mp_api.client import MPRester

    with MPRester(api_key=key) as mpr:
        structure = mpr.materials.search(material_ids=[mpid])
        label = mpid.replace("-", "_")
        return Batoms(label, from_pymatgen=structure)
