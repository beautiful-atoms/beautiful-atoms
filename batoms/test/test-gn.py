
def test_gn():
    import bpy
    from batoms.butils import removeAll
    from time import time
    from batoms.butils import get_nodes_by_name
    removeAll()
    bpy.ops.mesh.primitive_plane_add(size=2)
    obj = bpy.context.object
    obj.modifiers.new(name = 'test', type = 'NODES')
    tstart = time()
    gn = obj.modifiers[0]
    for i in range(10):
        tstart = time()
        for j in range(100):
            CompareSelect = get_nodes_by_name(gn.node_group.nodes, 
                        'select_%s_%s'%(i, j),
                        'GeometryNodeObjectInfo')
        print('time: %s'%(time() - tstart))


if __name__ == '__main__':
    test_gn()
    print('\n GN: All pass! \n')

"""
Blender 3.0.0 (hash f1cca3055776 built 2021-12-03 00:34:54)
Read prefs: /home/xing/.config/blender/3.0/config/userpref.blend
time: 0.013685464859008789
time: 0.041017770767211914
time: 0.07170844078063965
time: 0.10091757774353027
time: 0.1293931007385254
time: 0.15906882286071777
time: 0.18669986724853516
time: 0.2155919075012207
time: 0.23975181579589844
time: 0.2664918899536133

 GN: All pass!
 
 
 
Blender 3.1.0 Beta (hash 07514def194a built 2022-01-31 00:34:46)
Read prefs: /home/xing/.config/blender/3.1/config/userpref.blend
time: 0.02619338035583496
time: 0.07931733131408691
time: 0.13824701309204102
time: 0.2018442153930664
time: 0.2677316665649414
time: 0.34605860710144043
time: 0.4333493709564209
time: 0.5064260959625244
time: 0.5932877063751221
time: 0.6814370155334473

 GN: All pass! 


Blender 2.93.5 (hash a791bdabd0b2 built 2021-10-06 06:26:10)
Read prefs: /home/xing/.config/blender/2.93/config/userpref.blend
time: 0.008690834045410156
time: 0.023457050323486328
time: 0.034780263900756836
time: 0.0353388786315918
time: 0.04577922821044922
time: 0.05544900894165039
time: 0.06362771987915039
time: 0.0733647346496582
time: 0.08342456817626953
time: 0.0934457778930664

"""