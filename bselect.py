"""
select expression

all
none


Logical
not Sele
sele1 and sele2


Properties
element O


Coordinates
z<5



Chemcial classses

solvent: water
hydrogens: 
metals


"""


class Select():
    def __init__(self) -> None:
        pass
    @property    
    def data(self):
        return self.get_data()
    
    def __getitem__(self, index):
        name = tuple2string(index)
        item = self.find(name)
        if item is None:
            raise Exception('%s not in %s setting'%(name, self.name))
        return item
    
    def __setitem__(self, index, setdict):
        """
        Set properties
        """
        name = tuple2string(index)
        subset = self.find(name)
        if subset is None:
            subset = self.collection.add()
        subset.name = name
        for key, value in setdict.items():
            setattr(subset, key, value)
        subset.label = self.label
        subset.flag = True
    
    def add(self, label):
        pass
    
    def remove(self, index):
        """
        index: list of tuple
        """
        if isinstance(index, (str, int, float)):
            index = [(index)]
        if isinstance(index[0], (str, int, float)):
            index = [index]
        for key in index:
            name = tuple2string(key)
            i = self.collection.find(name)
            if i != -1:
                self.collection.remove(i)
            else:
                raise Exception('%s is not in %s'%(name, self.name))
    