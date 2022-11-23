import bpy
from batoms.batoms import Batoms


def get_active_bpy_data(btype = 'cell'):
    """Helper function.

    Args:
        btype (str, optional): _description_. Defaults to 'cell'.
    """
    def getter():
        """Get the collection of the active Batoms

        When get the attribute of Batoms object,
        if the attribute if saved in the Batoms.coll.batoms,
        we only need to read data form the colleciton,
        it is faster than get data from the Batoms itself.

        Returns:
            bpy.type.collection: _description_
        """
        context = bpy.context
        if context.object and context.object.batoms.type != 'OTHER':
            return getattr(bpy.data.collections[context.object.batoms.label], btype)
        return None
    return getter


def get_active_bpy_data_batoms(btype = 'cell'):
    """Helper function.

    Args:
        btype (str, optional): _description_. Defaults to 'cell'.
    """
    def getter():
        """

        Returns:
            bpy.type.collection: _description_
        """
        context = bpy.context
        if context.object and context.object.batoms.type != 'OTHER':
            return getattr(bpy.data.collections[context.object.batoms.label].batoms, btype)
        return None
    return getter


def get_active_module(module_name = 'cell'):
    """Helper function.

    Args:
        module_name (str, optional): _description_. Defaults to 'cell'.

    Returns:
        function: function to get the module
    """
    def getter():
        """Get the module of Batoms by name

        Returns:
            _type_: _description_
        """
        context = bpy.context
        if context.object and context.object.batoms.type != 'OTHER':
            mode = context.object.mode
            batoms = Batoms(label=context.object.batoms.label)
            bpy.context.view_layer.objects.active = batoms.obj
            bpy.ops.object.mode_set(mode=mode)
            return getattr(batoms, module_name)
        return None
    return getter

def set_module_attr(module_name):
    """Helper function.

    Args:
        module_name (_type_): _description_
    """
    def setter(key, value):
        module = get_active_module(module_name)()
        if module is not None:
            setattr(module, key, value)
            # bpy.context.view_layer.objects.active = module.obj
    return setter



def get_enum_attr(name, func):
    """Helper function to easily get enum property.

    Args:
        name (str): name of the attribute
    """

    def getter(self):
        batoms = func()
        if batoms is not None:
            return int(getattr(batoms, name))
        else:
            return 0

    return getter


def set_enum_attr(name, func):
    """Helper function to easily set enum property.

    Args:
        name (str): name of the attribute
    """

    def setter(self, value):
        items = self.bl_rna.properties[name].enum_items
        item = items[value]
        identifier = item.identifier
        self[name] = identifier
        func(name, value)

    return setter


def get_attr(name, func):
    """Helper function to easily get property.

    Args:
        name (str): name of the attribute
    """

    def getter(self):
        batoms = func()
        if batoms is not None:
            return getattr(batoms, name)
        else:
            prop = self.bl_rna.properties[name]
            return prop.default
    return getter


def set_attr(name, func):
    """Helper function to easily set property.

    Args:
        name (str): name of the attribute
    """

    def setter(self, value):
        # set value
        self[name] = value
        # set property of the module
        func(name, value)

    return setter
