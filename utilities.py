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

# Takes iterable ls1 and list ls2 such that ls1 is a sublist of ls2 and returns a sorted version of ls1 according to the ordering in ls2
def sort_subset_by_list(ls1, ls2):
    return sorted(ls1, key=lambda x: ls2.index(x))

def prettyPrint(formatObj, *args, depth=0, verbose):
    if not verbose:
        return
    if len(args) == 0:
        print(("\t" * depth) + str(formatObj))
    else:
        print(("\t" * depth) + formatObj.format(*args))

def causalEffectStr(x, y):
    return "P_{{{:}}}({:})".format(",".join(x), ",".join(y))

def conditionalProbStr(y, x):
    return "P({:}|{:})".format(",".join(y), ",".join(x))
