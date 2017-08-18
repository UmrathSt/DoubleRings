from MultipoleFields import HarmonicField
import numpy as np

class Scatterer(object):
    def __init__(self, origin):
        self.origin = origin
        assert type(origin) == np.ndarray
        self.distance = np.sqrt((origin**2).sum())
        self.phi = np.arctan2(origin[1], origin[0])
        self.theta = np.arccos(origin[2]/self.distance)
        assert origin.shape == (3,)

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
        assert isinstance(field, HarmonicField)
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
    print("Origin of S: ", Y_Wall.origin)
    field = HarmonicField(np.array([-1,0,1,-2,0,1]), 1, 1)
    print("field after scattering", Y_Wall.scatter(field))
