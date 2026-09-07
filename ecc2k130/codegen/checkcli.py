"""Every Modal entry-point parameter has to survive the round trip to the CLI.

Modal builds the flag as "--" + name.replace("_", "-"), and click derives the
keyword it hands back by lowercasing that and turning "-" into "_".  So a
capital letter never comes back: --minBlocksList arrives as minblockslist and
the call dies with "unexpected keyword argument".  The invariant is therefore
simply that an entry-point parameter is lowercase, and it is worth checking
because the failure is invisible until someone runs that entry point -- bench,
autotune, search and fanout were all unrunnable this way while validate and
merge, which happen to have only lowercase parameters, worked fine.

The repository is otherwise camelCase; this is the one boundary where an
outside tool dictates the spelling.
"""

import ast
import sys


def entryPointParams(path):
    tree = ast.parse(open(path).read())
    for node in ast.walk(tree):
        if not isinstance(node, ast.FunctionDef):
            continue
        for dec in node.decorator_list:
            if isinstance(dec, ast.Call) and getattr(dec.func, 'attr', '') == 'local_entrypoint':
                # Every kind of named parameter, not just positional-or-keyword:
                # a keyword-only one is exactly where somebody would add the
                # next flag, and missing it would make this check useless in
                # the one case it exists for.
                a = node.args
                for arg in list(getattr(a, 'posonlyargs', [])) + list(a.args) + list(a.kwonlyargs):
                    yield node.name, arg.arg


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'modal_app.py'
    bad = [(fn, p) for fn, p in entryPointParams(path) if p != p.lower()]
    for fn, p in bad:
        flag = '--' + p.replace('_', '-')
        print('%s: parameter %r becomes %s, which click hands back as %r'
              % (fn, p, flag, flag.lstrip('-').replace('-', '_').lower()))
    if bad:
        print('%d entry-point parameter(s) cannot be passed on the command line' % len(bad))
        return 1
    print('modal entry points: %d parameters, all round-trip'
          % sum(1 for _ in entryPointParams(path)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
