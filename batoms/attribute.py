import numpy as np




class Attribute:
    def __init__(self, name, dtype = None, shape = (1),
                domain = "POINT",
                array = None,
                batoms = None) -> None:
        self.name = name
        self.shape = shape
        self.domain = domain
        self.batoms = batoms
        self.dtype = dtype
        self.array = np.array(array)
        self.shape = self.array.shape
        self.subatt = {}
        self.scatter()
    
    def scatter(self):
        self.subatt = {}
        if len(self.shape) == 1:
            if self.dtype is None:
                self.dtype = self.get_dtype(self.array[0])
            self.subatt["{}".format(self.name)] = self.array
        if len(self.shape) == 2:
            if self.dtype is None:
                self.dtype = self.get_dtype(self.array[0][0])
            for i in range(self.shape[1]):
                self.subatt["{}{}".format(self.name, i)] = self.array[:, i]
        if len(self.shape) == 3:
            if self.dtype is None:
                self.dtype = self.get_dtype(self.array[0][0][0])
            for i in range(self.shape[1]):
                for j in range(self.shape[2]):
                    self.subatt["{}{}{}".format(self.name, i, j)] = self.array[:, i, j]

    def gather(self):
        pass

    
    def get_dtype(self, data):
        dtype = type(data)
        if np.issubdtype(dtype, int):
            dtype = 'INT'
        elif np.issubdtype(dtype, float):
            dtype = 'FLOAT'
        elif np.issubdtype(dtype, str):
            dtype = 'STRING'
        elif np.issubdtype(dtype, bool):
            dtype = 'BOOLEAN'
        else:
            raise KeyError('Attribute: {}, {} is not supported.'.format(self.name, dtype))
        return dtype

    def from_list(self, data):
        pass

