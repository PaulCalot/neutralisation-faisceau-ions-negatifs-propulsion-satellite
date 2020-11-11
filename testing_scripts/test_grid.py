from modules.particules import Particule
from modules.grid import Grid
from modules.vector import MyVector

my_grid = Grid(10,20,[10,20])

print("\nDefining several particules :")
part1 = Particule()
part2 = Particule(pos = MyVector(0.5,0.5))
part3 = Particule(pos = MyVector(1.5, 0))

print("\nTrying to get the charge of a particule ... {} [OK]".format(part1.get_charge()))

print("\nAdding those particules :")
my_grid.add(part1)
my_grid.add(part2)
my_grid.add(part3)

my_grid.to_string() 

print("\nChanging part 1 & 3 pos :")
old_pos1 = part1.get_pos()
part1.set_pos(MyVector(0.2,0))
old_pos3 = part3.get_pos()
part3.set_pos(MyVector(0.5,0.6))

my_grid.to_string()

print("\nUpdating part 1 & 3 in grid :")
my_grid.update(part1, old_pos1)
my_grid.update(part3, old_pos3)

my_grid.to_string()

print("\nTrying to have the closest particules ...")
my_parts = my_grid.get_closest_particules(part1)
for part in my_parts:
    print(part.to_string())

print("\nRemoving part 1 & 2 & 3 :")
my_grid.remove(part1)
my_grid.remove(part2)
my_grid.remove(part3)

my_grid.to_string()

print("\n[OK]")

print("\nTesting vectors ...")

v1 = MyVector(0.5,0.5)
v2 = MyVector(1.0,0.5)
alpha = 0.9

v3 = v1 + alpha*v2
print("V3 = ({},{})".format(v3.x, v3.y))