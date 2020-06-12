# Takes a list of lists or tuples and turns it into a single list
def flatten(alist):
    new_list = []
    for item in alist:
        if isinstance(item, (list, tuple)):
            new_list.extend(flatten(item))
        else:
            new_list.append(item)
    return new_list

#Takes a list of tuples, and then turns each tuple into a frozenset. Returns a list containing one copy of each set represented by a tuple.
def unique_tuples(ls):
    s = set()
    for t in ls:
        s.add(frozenset(t))
    return [tuple(t) for t in s]
