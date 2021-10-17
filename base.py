import bpy

class Base():
    def __init__(self, name):
        self.name = name
    @property
    def coll(self):
        return self.get_coll()
    def get_coll(self):
        return bpy.data.collections[self.name]