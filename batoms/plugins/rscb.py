import urllib.request
from batoms import Batoms
from batoms.pdbparser import read_pdb
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


def rscb_import(name):
    if not name.endswith('.pdb'):
        name += '.pdb'
    urllib.request.urlretrieve(
        'http://files.rcsb.org/download/%s' % name, name)
    atoms = read_pdb(name)
    batoms = Batoms('pdb_%s' % name[0:-4], from_ase=atoms)
    batoms.ribbon.draw()
    return batoms


if __name__ == '__main__':
    rscb_import('1ema')
