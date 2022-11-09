class POSet(object):
    def __init__(self, node=None):
        self.relations = frozenset()
        self.elements = set()

        if node is not None:
            self.relations = frozenset([(node, node)])
            self.elements.add(node)

    def series(self, poset):
        new_poset = POSet()
        # Union of the elements in the posets
        new_poset.elements = self.elements.union(poset.elements)
        # Cross product
        new_poset.relations = frozenset((a, b) for a in self.elements for b in poset.elements)
        # Union of relation sets
        new_poset.relations = new_poset.relations.union(self.relations.union(poset.relations))
        return new_poset

    def parallel(self, poset):
        new_poset = POSet()
        # Union of elements in posets
        new_poset.elements = self.elements.union(poset.elements)
        # Union of relation sets
        new_poset.relations.add(self.relations.union(poset.relations))
        return new_poset

    def __repr__(self):
        return "elements:\n{}\nrelation:\n{}".format(self.elements, self.relations)

class Node(object):
    """This is an noncomparable empty object representing the items in a poset.
    I've left it empty so that it can be modified for later implementations."""

    def __repr__(self):
        return "Node<{}>".format(id(self))
