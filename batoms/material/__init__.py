"""
module for materials
"""
import bpy

material_styles_dict = {
    'default': {'type': 'Principled BSDF',
                'inputs': {'Metallic': 0.10, 'Roughness': 0.20,
                           'Specular': 0.2}},
    'metallic': {'type': 'Principled BSDF',
                 'inputs': {'Metallic': 1.0, 'Roughness': 0.2,
                            'Specular': 0.2, 'Subsurface': 0.1}},
    'ceramic': {'type': 'Principled BSDF',
                'inputs': {'Metallic': 0.02, 'Roughness': 0.00,
                           'Specular': 0.5, 'Subsurface': 0.1}},
    'plastic': {'type': 'Principled BSDF',
                'inputs': {'Metallic': 0.00, 'Roughness': 1.00,
                           'Specular': 0.5}},
    'mirror': {'type': 'Principled BSDF',
               'inputs': {'Metallic': 0.99, 'Roughness': 0.01,
                          'Specular': 2.0}},
    'glass': {'type': 'Glass BSDF',
              'inputs': {'Roughness': 0.5, 'IOR': 1.45}},
    'emission': {'type': 'Emission',
                 'inputs': {'Strength': 5.0}},
}


def create_material(name,
                    color,
                    node_type='Principled BSDF',
                    node_inputs=None,
                    material_style='default',
                    backface_culling=True,
                    vertex_color=None,
                    color_by_attribute=None,
                    alpha = None,):
    """
    Creat material
    """
    if node_inputs is None:
        node_inputs = material_styles_dict[material_style]['inputs']
        node_type = material_styles_dict[material_style]['type']
    material = bpy.data.materials.new(name)
    material.diffuse_color = color
    for key, value in node_inputs.items():
        if hasattr(material, key):
            setattr(material, key, value)
    material.blend_method = 'BLEND'
    material.use_nodes = True
    #
    nodes = material.node_tree.nodes
    mat_links = material.node_tree.links
    if node_type == 'Glass BSDF':
        nodes.new('ShaderNodeBsdfGlass')
    elif node_type == 'Emission':
        nodes.new('ShaderNodeEmission')
    node = nodes[node_type]
    node.inputs[0].default_value = color
    # Alpha transpancy
    if 'Alpha' in node.inputs:
        if alpha is not None:
            node.inputs['Alpha'].default_value = alpha
        else:
            node.inputs['Alpha'].default_value = color[3]

    for key, value in node_inputs.items():
        node.inputs[key].default_value = value
    if backface_culling:
        material.use_backface_culling = True
    material.show_transparent_back = False

    # vertex color
    if vertex_color not in [None, 'None']:
        if bpy.data.version > (3, 0, 0):
            vcol = nodes.new(type="ShaderNodeVertexColor")
        else:
            vcol = nodes.new(type="ShaderNodeVertexColor")
        vcol.layer_name = vertex_color
        mat_links.new(vcol.outputs['Color'], node.inputs['Base Color'])
        mat_links.new(vcol.outputs['Alpha'], node.inputs['Alpha'])
    # color by attribute
    if color_by_attribute is not None:
        attr_node = nodes.new(type="ShaderNodeAttribute")
        attr_node.attribute_name = color_by_attribute['attribute_name']
        if 'ValToRGB' in color_by_attribute:
            color_ramp = nodes.new(type="ShaderNodeValToRGB")
            color_ramp.color_ramp.elements[0].color = color_by_attribute['ValToRGB'][0]
            color_ramp.color_ramp.elements[1].color = color_by_attribute['ValToRGB'][1]
            mat_links.new(attr_node.outputs['Fac'], color_ramp.inputs['Fac'])
            mat_links.new(color_ramp.outputs['Color'], node.inputs['Base Color'])
            mat_links.new(color_ramp.outputs['Alpha'], node.inputs['Alpha'])
        else:
            mat_links.new(attr_node.outputs['Color'], node.inputs['Base Color'])
        # mat_links.new(vcol.outputs['Alpha'], node.inputs['Alpha'])

    return material


def clean_node_tree(nodes):
    """
    Clear all nodes except the output.
    """
    for node in list(nodes):
        if not node.type == 'OUTPUT_MATERIAL':
            nodes.remove(node)
