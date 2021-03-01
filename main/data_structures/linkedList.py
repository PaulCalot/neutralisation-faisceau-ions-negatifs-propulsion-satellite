    
# creating the linkedlist we'll use to store our particules inside the grid
    # https://www.codefellows.org/blog/implementing-a-singly-linked-list-in-python/
class Node(object):

    def __init__(self, data=None, next_node=None):
        self.data = data
        self.next_node = next_node

    def get_data(self):
        return self.data

    def get_next(self):
        return self.next_node

    def set_next(self, new_next):
        self.next_node = new_next


class LinkedList(object):
    def __init__(self, head=None):
        self.head = head
        self.size = 0
        if(head!=None):
            self.size = 1

    def get_head(self):
        return self.head

    def insert(self, data):
        self.size += 1
        new_node = Node(data)
        new_node.set_next(self.head)
        self.head = new_node
    
    # depreciated since we store the size now
    def size(self):
        current = self.head
        count = 0
        while current:
            count += 1
            current = current.get_next()
        return count

    def search(self, data):
        current = self.head
        found = False
        while current and found is False:
            if current.get_data() == data:
                found = True
            else:
                current = current.get_next()
        if current is None:
            raise ValueError("Data not in list")
        return current

    def delete(self, data):
        current = self.head
        previous = None
        found = False
        while current and found is False:
            if current.get_data() == data:
                found = True
                self.size -= 1
            else:
                previous = current
                current = current.get_next()
        if current is None:
            raise ValueError("Data not in list")
        if previous is None:
            self.head = current.get_next()
        else:
            previous.set_next(current.get_next())
    
    def get_pair(self, i, j):
        current = self.head
        assert(i<self.size and j<self.size)
        min_ij = min(i,j)
        max_ij = max(i,j)
        for k in range(min_ij):
            current = current.get_next()
        pair1 = current
        for k in range(max_ij-min_ij):
            current = current.get_next()
        pair2 = current
        return pair1.get_data(), pair2.get_data()
        
    # -------------- getter and setter ---------------- #
    def get_size(self):
        return self.size

    def to_list(self):
        current = self.head
        L = []
        while current :
            L.append(current.get_data())  
            current = current.get_next()
        return L