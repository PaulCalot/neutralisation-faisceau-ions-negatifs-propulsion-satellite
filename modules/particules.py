# import
import warnings

# "my classes" import
from  .vector import MyVector

# constant
NUCLEON_MASS = 1.672e-27 # kg
ELECTRON_MASS = 9.11e-31
IODINE_RADIUS = 1.98e-10 # 198 pm

def get_mass_part(electrons_nb, protons_number, neutrons_number):
    return (neutrons_number+protons_number)*NUCLEON_MASS+electrons_nb*ELECTRON_MASS

IODINE_MASS = get_mass_part(53, 53, 88) # Iodine : 53 protons, 53 electrons, 131 nucleons => 88 neutrons

class Particule(object):
    # TODO : consider adding an id to the particule ?

    e = 1.602e-19 # C
    types = ["I", "I2", "I-", "I+", "e"]
    def __init__(self, charge = 0, mass = IODINE_MASS, pos = MyVector(0,0,0), speed = MyVector(0,0,0), \
        part_type = "I", radius = IODINE_RADIUS, verbose = False, status = -1):
        """ Create a particule, by default an neutral atom of mass ... (Iode).

        Args:
            charge (int): the charge of the particules will then be charge * e 
            mass (float): mass of the particule
            pos (MyVector): position of the particules
            speed (MyVector): speed of the particules
            type (String): a name for the given particule type - 
                           Available types : I, I2, I-, I+, e
        """
        
        #if(type(charge)!=int):
        #    warnings.warn("Charge argument type *int* expected but got {}.".format(type(charge)))
        if(not part_type in self.types):
            warnings.warn("Type not recognized. You should choose amongst {}.".format(self.types))
        if(len(speed)==2):
            speed = MyVector(speed.x, speed.y, 0) # so it's compatible with 3D MyVector, even if we don't use it right now.
        
        self.charge = charge
        self.mass = mass
        self.pos = pos
        self.speed = speed
        self.part_type = part_type
        self.radius = radius
        self.status = status
        if(verbose):
            print("Created a " + self.to_string())
    
    # ------------------------ Particules functions -------------------------- #
    def lose_charge(self):
        # see if there is another possibility
        if(self.part_type=="I-"):
            self.part_type = "I"
            self.charge = 0.0
    # ------------------------ Utils functions ------------------------------ #

    def multiply_speed(self, scalar):
        """ Multiply the speed by a scalar. Useful for updating speed when colliding with a wall or another particule.

        Args:
            scalar (float): the other scalar by which we multiply the speed.
        """
        assert(type(scalar) == type(1) or type(scalar) == type(1.0))
        self.speed = self.speed.__mul__(scalar) 

    def rotate_speed_2D(self, theta):
        """ Rotate the supposed 2D speed by theta (around +z axis)
        Theta in degree.
        Args:
            theta (float): the given angle in degree.
        """
        assert(type(theta) == type(1) or type(theta) == type(1.0))
        new_2D_speed = MyVector(self.speed.x, self.speed.y).rotate(theta)
        self.speed = MyVector(new_2D_speed.x, new_2D_speed.y, self.speed.z)
        
    # ------------------------ Getter and Setter ----------------------------- #

        # mass
    def get_mass(self):
        return self.mass
    
    def set_mass(self, new_mass):
        self.mass = new_mass

        # charge
    def get_charge(self):
        return self.charge
    
    def set_charge(self, new_charge):
        #assert(type(new_charge)==int)
        self.charge = new_charge
        
    #def get_true_charge(self):
    #    return self.charge *self.e # another possibility is to save the charge as charge*e directly
                                  # and avoid recomputing each time this value (it can be costly) 

        # position
    def get_pos(self):
        return self.pos

    def set_pos(self, new_pos):
        assert type(new_pos) == MyVector
        self.pos = new_pos
    
        # position
    def get_2D_pos(self):
        return MyVector(self.pos.x, self.pos.y)

    def set_2D_pos(self, new_pos):
        assert type(new_pos) == MyVector
        self.pos = MyVector(new_pos.x, new_pos.y, self.pos.z)

        # speed
    def get_speed(self):
        return self.speed

    def set_speed(self, new_speed):
        assert type(new_speed) == MyVector
        self.speed = new_speed

        # 2D speed
    def get_2D_speed(self):
        return MyVector(self.speed.x, self.speed.y)

    def set_2D_speed(self, new_speed):
        assert type(new_speed) == MyVector
        self.speed = MyVector(new_speed.x, new_speed.y, self.speed.z)

        # particule type
    def get_part_type(self):
        return self.part_type

    def set_part_type(self, new_part_type):
        self.part_type = new_part_type

        # radius
    def get_radius(self):
        return self.radius

    def set_radius(self,new_radius):
        self.radius = new_radius

        # status
    def get_status(self):
        return self.status

    def set_status(self, new_status):
        self.status = new_status

    # ----------------- to String -------------------- #

    def to_string(self):
        return "{} particule of charge {} C and mass {} kg in position ({},{}) m with speed ({},{}) m/s" \
            .format(self.part_type, self.charge, self.mass, round(self.pos.x,3), round(self.pos.y,3), round(self.speed.x,3), round(self.speed.y,3))