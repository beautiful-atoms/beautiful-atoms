import bpy
import console_python
from batoms.batoms import Batoms

def console_hook():
    """add batoms to namespace of python console
    """
    from batoms.utils.butils import read_batoms_list
    import console_python
    items = read_batoms_list()
    for area in bpy.context.screen.areas:
        if area.type == 'CONSOLE':
            for region in area.regions:
                if region.type == 'WINDOW':
                    console, stdout, stderr = console_python.get_console(hash(region))
                    for item in items:
                        item = item.replace('-', '_')
                        item = item.replace('.', '')
                        if item[:1].isdigit():
                            item = 'b_' + item
                        if item not in console.locals:
                            console.locals[item] = Batoms(item)


def register_hook():
    console_python.execute.hooks.append((console_hook, ()))

def unregister_hook():
    console_python.execute.hooks.remove((console_hook, ()))
