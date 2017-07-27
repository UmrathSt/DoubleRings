import numpy as np
from MultipoleFields import MultipoleField

class Scatterer(object):
    def __init__(self, origin):
        self.origin = origin
        assert type(origin) == np.ndarray
        self.distance = np.sqrt((origin**2).sum())
        self.phi = np.arctan2(origin[1]/origin[0])
        self.theta = np.arccos(origin[3]/self.distance)

ssert origin.shape == (3,)

    def scatter(self, field, origin=np.array([0,0,0])):
        """ scatter a multipole field and return the total field
            after scattering it as a multipole-coefficient vector
            at the position origin
        """
        raise NotImplementedError("Must be implemented within a derived class")

class Sphere(Scatterer):
    def __init__(self, origin):
        Scatterer.__init__(self, origin)

    def scatter(self, field, origin = np.array([0,0,0])):
        assert isinstance(field, MultipoleField)
        return field

class PECWall(Scatterer):
    def __init__(self, origin):
        Scatterer.__init__(self, origin)
    
    def scatter(self, field, origin = np.array([0,0,0])):
        if not origin[0] == 0 and origin[1] == 0:
            field.rotate(-self.theta, -self.phi)
            field.z_translate(kd, 1, 1)
        return -field.coeffs

class PMCWall(Scatterer):
    def __init__(self, origin):
        Scatterer.__init__(self, origin)

    def scatter(self, field, origin = np.array([0,0,0])):
        if not origin[0] == 0 and origin[1] == 0:
            field.rotate(-self.theta, -self.phi)
            field.z_translate(kd, 1, 1)
            
        return field.coeffs

if __name__ == "__main__":

    Y_Wall = PECWall(np.array([0,0,2]))
    print("Origin of S: ", S.origin)
    print("field after scattering", S.scatter([1,2,3]))
