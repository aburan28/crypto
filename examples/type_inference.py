"""A tiny type-inference engine, explained simply.

THE WHOLE IDEA IN ONE BREATH
----------------------------
We have a function whose parameter/return types may contain "generics"
(placeholders like T1, T2). When the function is *called* with real
argument types, we want to figure out the real return type.

Example:
    function type:   [T1, T2, int, T1] -> [T1, T2]
    called with:     [int, char, int, int]
    answer:          [int, char]

We do it in two passes:
    1. UNIFY  : walk params next to the real args, and learn what each
                generic means.        ->  {T1: int, T2: char}
    2. SUBSTITUTE : take the return type and plug those meanings in.

That's it. Everything below is just bookkeeping for those two passes.


A TYPE IS ONE OF THREE THINGS
-----------------------------
    primitive : a concrete type, e.g. "int", "char", "float"
    generic   : a placeholder, e.g. "T1", "T2"  (anything starting with "T")
    tuple     : an ordered group of other types, e.g. [int, T1]
"""


class Node:
    """One node in a type tree (a primitive, a generic, or a tuple)."""

    def __init__(self, value=None, children=None):
        self.value = value          # the name, for primitives/generics
        self.children = children     # the parts, for tuples

    # --- "what kind of type am I?" -------------------------------------
    @property
    def is_generic(self):
        # A generic is a bare name that starts with "T", like T1.
        return self.value is not None and self.value.startswith("T")

    @property
    def is_primitive(self):
        # A primitive is any other bare name, like int / char / float.
        return self.value is not None and not self.value.startswith("T")

    @property
    def is_tuple(self):
        # A tuple has no name of its own, just a list of child types.
        return self.children is not None

    def __str__(self):
        if self.is_tuple:
            return "[ " + ", ".join(str(c) for c in self.children) + " ]"
        return self.value


# Short, named constructors so the rest of the code reads in English.
def generic(name):       return Node(value=name)        # T1, T2, ...
def prim(name):          return Node(value=name)        # int, char, ...
def tup(*children):      return Node(children=list(children))


class Func:
    """A function type: a list of parameter types and one return type."""

    def __init__(self, params, returns):
        self.params = params
        self.returns = returns

    def __str__(self):
        params = "; ".join(str(p) for p in self.params)
        return f"[{params}] -> {self.returns}"


# ---------------------------------------------------------------------------
# PASS 1: UNIFY  -- learn what each generic stands for.
# ---------------------------------------------------------------------------
def unify(param, arg, subs):
    """Match a parameter type against a real argument type.

    `subs` is the dictionary we are filling in, e.g. {"T1": <int node>}.
    Raises if the two types can't possibly be the same shape.
    """
    # A generic learns its meaning from the argument.
    if param.is_generic:
        if arg.is_generic:
            raise AssertionError("cannot match a generic with another generic")
        seen = subs.get(param.value)
        if seen is None:
            subs[param.value] = arg              # first sighting: record it
        elif str(seen) != str(arg):
            raise AssertionError(f"{param.value} is both {seen} and {arg}")
        return

    # A primitive must match the exact same primitive.
    if param.is_primitive:
        if not (arg.is_primitive and arg.value == param.value):
            raise AssertionError(f"expected {param.value}, got {arg}")
        return

    # A tuple must match a tuple of the same length, element by element.
    if param.is_tuple:
        if not arg.is_tuple or len(param.children) != len(arg.children):
            raise AssertionError(f"cannot match {param} with {arg}")
        for p, a in zip(param.children, arg.children):
            unify(p, a, subs)
        return

    raise AssertionError("invalid node")


# ---------------------------------------------------------------------------
# PASS 2: SUBSTITUTE -- rebuild a type with generics replaced by what we learned.
# ---------------------------------------------------------------------------
def substitute(node, subs):
    """Return a copy of `node` with every generic swapped for its real type."""
    if node.is_generic:
        return subs.get(node.value, generic(node.value))   # leave unknown ones as-is
    if node.is_primitive:
        return prim(node.value)
    if node.is_tuple:
        return tup(*(substitute(c, subs) for c in node.children))
    raise AssertionError("invalid node")


def infer_return(func, args):
    """Given a function type and the real argument types, return the real return type."""
    if len(func.params) != len(args):
        raise AssertionError("wrong number of arguments")
    subs = {}
    for param, arg in zip(func.params, args):   # pass 1: fill in T1, T2, ...
        unify(param, arg, subs)
    return substitute(func.returns, subs)        # pass 2: plug them into the return type


# ===========================================================================
# Tests
# ===========================================================================
if __name__ == "__main__":
    print("Node toStr")
    print(prim("int"))                                   # int
    print(tup(prim("int"), prim("char")))                # [ int, char ]
    print(tup(prim("int"), tup(prim("int"), prim("char"))))  # [ int, [ int, char ] ]

    print("\nFunc toStr")
    print(Func([prim("int"), prim("int")], prim("int")))           # [int; int] -> int
    print(Func([prim("int"), tup(prim("int"), prim("char")), prim("int")], prim("int")))

    example = Func(
        [prim("int"), tup(prim("int"), generic("T1")), generic("T2")],
        tup(prim("char"), generic("T2"), tup(prim("float"), prim("float"))),
    )
    print(example)
    assert str(example) == "[int; [ int, T1 ]; T2] -> [ char, T2, [ float, float ] ]"

    print("\nInfer return")
    # [T1, T2, int, T1] -> [T1, T2]   called with   [int, char, int, int]
    f = Func([generic("T1"), generic("T2"), prim("int"), generic("T1")],
             tup(generic("T1"), generic("T2")))
    result = infer_return(f, [prim("int"), prim("char"), prim("int"), prim("int")])
    print(result)
    assert str(result) == "[ int, char ]"

    print("\nAll tests passed.")
