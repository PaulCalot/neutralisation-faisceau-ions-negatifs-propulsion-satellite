from modules.vector import MyVector

v1 = MyVector(-1,-1,0)
print(v1.normalize())

saving_period = 5
for k in range(1000):
    if(k%saving_period==0):
        print(k)