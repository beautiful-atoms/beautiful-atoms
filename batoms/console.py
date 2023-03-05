import bpy
import console_python
from batoms.batoms import Batoms
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


def console_hook():
    """add batoms to namespace of python console
    """
    from batoms.utils.butils import read_batoms_list
    items = read_batoms_list()
    for area in bpy.context.screen.areas:
        if area.type == 'CONSOLE':
            for region in area.regions:
                if region.type == 'WINDOW':
                    console, stdout, stderr = console_python.get_console(hash(region))
                    if "Batoms" not in console.locals:
                        console.push("from batoms import Batoms")
                    for item in items:
                        item = item.replace('-', '_')
                        item = item.replace('.', '')
                        if item[:1].isdigit():
                            item = 'b_' + item
                        if item not in console.locals:
                            logger.info("Python console: Add Batoms {}".format(item))
                            console.push("{}=Batoms('{}')".format(item, item))
                            console.push("")
                            # console.locals[item] = Batoms(item)


def register_hook():
    console_python.execute.hooks.append((console_hook, ()))

def unregister_hook():
    console_python.execute.hooks.remove((console_hook, ()))
