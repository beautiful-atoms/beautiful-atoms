def test_species():
    """
    """
from batoms.butils import removeAll
from batoms.bspecies import Bspecies
import numpy as np
removeAll()
h2o = Bspecies('h2o', species = {'O':{'elements':{'O':0.8, 'N': 0.2}}, 
                                 'H':{'elements': {'H': 0.8}}})
h2o['O'] = {'elements': {'Al':0.8, 'Si':0.2}}