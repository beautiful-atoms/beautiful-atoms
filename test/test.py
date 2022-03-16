import bpy
import numpy as np
from batoms import Batoms

bpy.ops.batoms.delete()
bpy.ops.batoms.bulk_add(cubic = True)
au = Batoms('Au')
au.rotate(45, 'Z')
print(au.positions)
assert np.allclose(au.positions[1], np.array([1.44249, 1.44249, 2.04]))
#
au.repeat([2, 1, 1])