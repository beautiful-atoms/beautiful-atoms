
# This allows you to right click on a button and link to documentation
def batoms_manual_map():
    url_manual_prefix = "https://wiki.fysik.dtu.dk/ase/ase/"
    url_manual_mapping = (
        ("bpy.ops.batoms.atoms_add", "atoms.html?highlight=atoms#ase.Atoms"),
        ("bpy.ops.batoms.molecule_add",
         "build/build.html?highlight=molecule#ase.build.molecule"),
        ("bpy.ops.batoms.bulk_add", "build/build.html?highlight=bulk#ase.build.bulk"),
        ("bpy.ops.batoms.surface_add",
         "build/build.html?highlight=surface#ase.build.surface"),
        ("bpy.ops.batoms.root_surface_add",
         "build/build.html?highlight=surface#ase.build.root_surface"),
    )
    return url_manual_prefix, url_manual_mapping


def batoms_ase_manual_map():
    url_manual_prefix = "https://wiki.fysik.dtu.dk/ase/ase/"
    url_manual_mapping = (
        ("bpy.ops.surface.fcc100_add",
         "build/surface.html?highlight=fcc100#ase.build.fcc100"),
        ("bpy.ops.surface.fcc110_add",
         "build/surface.html?highlight=fcc110#ase.build.fcc110"),
        ("bpy.ops.surface.fcc111_add",
         "build/surface.html?highlight=fcc111#ase.build.fcc111"),
        ("bpy.ops.surface.fcc211_add",
         "build/surface.html?highlight=fcc211#ase.build.fcc211"),
        ("bpy.ops.surface.fcc111_root_add",
         "build/surface.html?highlight=fcc111_root#ase.build.fcc111_root"),
        ("bpy.ops.surface.bcc100_add",
         "build/surface.html?highlight=bcc100#ase.build.bcc100"),
        ("bpy.ops.surface.bcc110_add",
         "build/surface.html?highlight=bcc110#ase.build.bcc110"),
        ("bpy.ops.surface.bcc111_add",
         "build/surface.html?highlight=bcc111#ase.build.bcc111"),
        ("bpy.ops.surface.bcc111_root_add",
         "build/surface.html?highlight=bcc111_root#ase.build.bcc111_root"),
        ("bpy.ops.surface.hcp0001_add",
         "build/surface.html?highlight=hcp0001#ase.build.hcp0001"),
        ("bpy.ops.surface.hcp10m10_add",
         "build/surface.html?highlight=hcp10m10#ase.build.hcp10m10"),
        ("bpy.ops.surface.hcp0001_root_add",
         "build/surface.html?highlight=hcp0001_root#ase.build.hcp0001_root"),
        ("bpy.ops.surface.diamond100_add",
         "build/surface.html?highlight=diamond100#ase.build.diamond100"),
        ("bpy.ops.surface.diamond111_add",
         "build/surface.html?highlight=diamond111#ase.build.diamond111"),
        ("bpy.ops.nano.nanotube_add",
         "build/build.html?highlight=nanotube#ase.build.nanotube"),
        ("bpy.ops.nano.nanoribbon_add",
         "build/build.html?highlight=nanoribbon#ase.build.nanoribbon"),
        ("bpy.ops.nano.decahedron_add",
         "cluster/cluster.html?highlight=decahedron#ase.cluster.Decahedron"),
        ("bpy.ops.nano.icosahedron_add",
         "cluster/cluster.html?highlight=icosahedron#ase.cluster.Icosahedron"),
        ("bpy.ops.nano.octahedron_add",
         "cluster/cluster.html?highlight=octahedron#ase.cluster.Octahedron"),
    )
    return url_manual_prefix, url_manual_mapping
