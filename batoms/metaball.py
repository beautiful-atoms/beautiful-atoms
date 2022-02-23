
    def build_metaball(self, render_resolution = 0.2,
                resolution = 0.4,
                threshold = 1e-4,
                update_method = 'FAST'):
        if self.sas_name in bpy.data.metaballs:
            mb = bpy.data.metaballs.get(self.sas_name)
            bpy.data.metaballs.remove(mb, do_unlink = True)
        if self.sas_name in bpy.data.materials:
            mat = bpy.data.materials.get(self.sas_name)
            bpy.data.materials.remove(mat, do_unlink = True)
        mb = bpy.data.metaballs.new(self.sas_name)
        mb.render_resolution = render_resolution
        mb.resolution = resolution
        mb.update_method = update_method
        mb.threshold = threshold
        obj = bpy.data.objects.new(self.sas_name, mb)
        coll = self.batoms.coll.children['%s_surface'%self.label]
        coll.objects.link(obj)
        return obj
    
    def draw_SAS_mb(self, render_resolution = 0.2,
                resolution = 0.4,
                threshold = 1e-4,
                stiffness = 1,
                indices = None,
                update_method = 'FAST'
                ):
        """
        Algorithm: Metaball
        Computes density from given metaball at given position.
        Metaball equation is: 
        dens = `(1 - r^2 / R^2)^3 * s`
        r = distance from center
        R = metaball radius
        s - metaball stiffness
        field = threshold - dens;
        """
        from batoms.material import create_material
        stiffness = min(stiffness, 10)
        frames = self.batoms.frames
        n = len(frames[0])
        if indices is None:
            indices = range(n)
        radii = self.batoms.radii_vdw
        tstart = time()
        #
        self.build_metaball(render_resolution = render_resolution,
                resolution = resolution,
                threshold = threshold,
                update_method = update_method)
        color = default_colors[0]
        mat = create_material(self.sas_name,
                        # material_style = 'plastic',
                        color = color)
        obj = self.sas_obj
        mb = obj.data
        obj.data.materials.append(mat)
        atoms = frames[0]
        positions = atoms.positions
        scale = (1 - (threshold/stiffness)**(1.0/3))**(0.5)
        # print('scale: %s', scale)
        mb.elements.clear()
        for i in indices:
            mbele = mb.elements.new(type = 'BALL')
            mbele.co = positions[i]
            mbele.radius = (radii[i] + self.probe)/scale
            mbele.stiffness = stiffness
        # add frames
        nframe = len(frames)
        self.nframe = nframe
        for i in range(1, nframe):
            positions = frames[i].positions
            for j in indices:
                mb.elements[j].co = positions[j]
                mb.elements[j].keyframe_insert(data_path='co', frame = i)
        print('Build_SAS: %s'%(time() - tstart))
        return obj