"""
"""
from batoms import Batoms
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


def pubchem_search(cid):
    from ase.data.pubchem import pubchem_search
    # import ssl
    # ssl._create_default_https_context = ssl._create_unverified_context
    data = pubchem_search(cid=cid)
    batoms = Batoms(label='cid_%s' % cid, from_ase=data.atoms)
    return batoms
