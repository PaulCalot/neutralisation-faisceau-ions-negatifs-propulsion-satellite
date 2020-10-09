# local import

from  .vector import Vector

NUCLEON_MASS = 1.672e-27 # kg
ELECTRON_MASS = 9.11e-31

def get_mass_part(electrons_nb, protons_number, neutrons_number):
    return (neutrons_number+protons_number)*NUCLEON_MASS+electrons_nb*ELECTRON_MASS

IODINE_MASS = get_mass_part(53, 53, 88) # Iodine : 53 protons, 53 electrons, 131 nucleons => 88 neutrons

class Particule(object):
    
    e = 1.602e-19 # C

    def __init__(self, charge = 0, mass = IODINE_MASS, pos = Vector(0,0), speed = Vector(0,0), \
        part_type = "I", verbose = True):
        """ Create a particule, by default an neutral atom of mass ... (Iode).

        Args:
            charge (int): the charge of the particules will then be charge * e 
            mass (float): mass of the particule
            x (float): position along x axis
            y (float): position along x axis
            vx (float): speed along y axis
            vy (float): speed along y axis
            type (String): a name for the given particule type - 
                           Available types : I, I2, Im, I2p, e
        """
        
        # asserting we are indeed giving to charge an int
        assert(type(charge)==int)

        self.charge = charge
        self.mass = mass
        self.pos = pos
        self.speed = speed
        self.part_type = part_type

        if(verbose):
            print("Created a " + self.to_string())
    
    # ------------------------ Utils functions ------------------------------ #

    def multiply_speed(self, scalar):
        """ Multiply the speed by a scalar. Useful for updating speed when colliding with a wall or another particule.

        Args:
            scalar (float): the other scalar by which we multiply the speed.
        """
        assert(type(scalar) == type(1) or type(scalar) == type(1.0))
        self.speed = self.speed.__mul__(scalar) 

    def rotate_2D(self, theta):
        """ Rotate the supposed 2D speed by theta (around +z axis)

        Args:
            theta (float): the given angle
        """
        assert(type(theta) == type(1) or type(theta) == type(1.0))
        self.speed = self.speed.rotate(theta)
    # ------------------------ Getter and Setter ----------------------------- #

        # mass
    def get_mass(self):
        return self.mass
    def set_mass(self, new_mass):
        self.mass = new_mass

        # charge
    def get_charge(self):
        return self.charge*self.e # another possibility is to save the charge as charge*e directly
                                  # and avoid recomputing each time this value (it can be costly) 
    def set_charge(self, new_charge):
        assert(type(new_charge)==int)
        self.charge = new_charge

        # position
    def get_pos(self):
        return self.pos

    def set_pos(self, new_pos):
        assert type(new_pos) == Vector
        self.pos = new_pos
    
        # speed
    def get_speed(self):
        return self.speed

    def set_speed(self, new_speed):
        assert type(new_speed) == Vector
        self.speed = new_speed

        # particule type
    def get_part_type(self):
        return self.part_type

    def set_type(self, new_part_type):
        self.part_type = new_part_type


    # ----------------- to String -------------------- #

    def to_string(self):
        return "{} particule of charge {} C and mass {} kg in position ({},{}) m with speed ({},{}) m/s" \
            .format(self.part_type, self.charge, self.mass,self.pos.x, self.pos.y, self.speed.x, self.speed.y)