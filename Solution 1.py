#!/usr/bin/env python3

## Default executable of a SAT solver (do not change this)
defSATsolver="z3"

## Change this to an executable SAT solver if z3 is not in your PATH or else
## Example (Linux): SATsolver="/home/user/z3-4.13/bin/z3"
## You can also include command-line options if necessary
SATsolver=defSATsolver

import sys
from subprocess import Popen
from subprocess import PIPE
import re
import random
import os
import shutil

gVarNumberToName = ["invalid"]
gVarNameToNumber = {}

def closed_range(start, stop, step=1):
    dir = 1 if (step > 0) else -1
    return range(start, stop + dir, step)

def varCount():
    global gVarNumberToName
    return len(gVarNumberToName) - 1

def allVarNumbers():
    return closed_range(1, varCount())

def varNumberToName(num):
    global gVarNumberToName
    return gVarNumberToName[num]

def varNameToNumber(name):
    global gVarNameToNumber
    return gVarNameToNumber[name]

def addVarName(name):
    global gVarNumberToName
    global gVarNameToNumber
    gVarNumberToName.append(name)
    gVarNameToNumber[name] = varCount()

# def printClause(clause):
#     print(map(lambda x: "%s%s" % (x < 0 and eval("'-'") or eval ("''"), varNumberToName(abs(x))) , clause))

def getVarNumber(**kwargs):
    return varNameToNumber(getVarName(**kwargs))

def getVarName(**kwargs):
    v = kwargs['v']
    c = kwargs['c']
    return f"x({v},{c})"

def genVarNames(**kwargs):
    # Extracting the number of vertices and colors
    n = kwargs['n']
    k = kwargs['k']

    # Generating the variables as tuples of (Vertex number, color number)
    for v in closed_range(1, n):
        for c in closed_range(1, k):
            name = getVarName(v=v, c=c)
            addVarName(name)

def genClauses(**kwargs):
    clauses = []

    n = kwargs['n']
    k = kwargs['k']
    edges = kwargs['edges']

    # Each vertex has at least one color
    for v in closed_range(1, n):
        clauses.append([getVarNumber(v=v, c=c) for c in closed_range(1, k)])

    # Each vertex has at most one color
    for v in closed_range(1, n):
        for c1 in closed_range(1, k):
            for c2 in closed_range(c1 + 1, k):
                clauses.append([-getVarNumber(v=v, c=c1), -getVarNumber(v=v, c=c2)])

    # No adjacent vertices share the same color
    for (u, v) in edges:
        for c in closed_range(1, k):
            if u != v:
                clauses.append([-getVarNumber(v=u, c=c), -getVarNumber(v=v, c=c)])

    return clauses

## A helper function to print the cnf header (do not modify)
def getDimacsHeader(clauses):
    cnt = varCount()
    n = len(clauses)
    str = ""
    for num in allVarNumbers():
        varName = varNumberToName(num)
        str += "c %d ~ %s\n" % (num, varName)
    for cl in clauses:
        print("c ", end='')
        for l in cl:
            print(("!" if (l < 0) else " ") + varNumberToName(abs(l)), "", end='')
        print("")
    print("")
    str += "p cnf %d %d" % (cnt, n)
    return str

## A helper function to print a set of clauses in CNF (do not modify)
def toDimacsCnf(clauses):
    return "\n".join(map(lambda x: "%s 0" % " ".join(map(str, x)), clauses))

## A helper function to print only the satisfied variables in human-readable format (do not modify)
def printResult(res):
    print(res)
    res = res.strip().split('\n')

    # If it was satisfiable, we want to have the assignment printed out
    if res[0] != "s SATISFIABLE":
        return

    # First get the assignment, which is on the second line of the file, and split it on spaces
    # Read the solution
    asgn = map(int, res[1].split()[1:])
    # Then get the variables that are positive, and get their names.
    # This way we know that everything not printed is false.
    # The last element in asgn is the trailing zero and we can ignore it

    # Convert the solution to our names
    facts = map(lambda x: varNumberToName(abs(x)), filter(lambda x: x > 0, asgn))

    # Print the solution
    print("c SOLUTION:")
    for f in facts:
        print("c", f)

def extractVertexColoring(res):
    # Parsing the result found by the solver into a list of lines
    res = res.strip().split('\n')

    # If no solution found, return
    if res[0] != "sat":
        print("c No solution found.")
        return

    # Keeping from the found solution only the variables that are valid (value = 1)
    asgn = map(int, res[1].split()[1:])
    true_vars = filter(lambda x: x > 0, asgn)

    # Using only the preserved variables numbers, go back to tuples (Vertex number, color number)
    coloring = {}
    for var_num in true_vars:
        # Converting variables number to name, in the form (Vertex number, color number)
        var_name = varNumberToName(abs(var_num))
        # Expecting format: x(v,c) -> See function getVarName
        match = re.match(r"x\((\d+),(\d+)\)", var_name)
        # If the regular expression correctly extracted the tuple data, add it to the dictionary of coloring
        if match:
            v = int(match.group(1))
            c = int(match.group(2))
            coloring[v] = c

    print("c FINAL COLORING:")
    for v in sorted(coloring.keys()):
        print(f"c Vertex {v} → Color {coloring[v]}")

## This function is invoked when the python script is run directly and not imported
if __name__ == '__main__':

    # Verifying that the solver path is correctly set
    path = shutil.which(SATsolver.split()[0])
    if path is None:
        if SATsolver == defSATsolver:
            print("Set the path to a SAT solver via SATsolver variable on line 9 of this file (%s)" % sys.argv[0])
        else:
            print("Path '%s' does not exist or is not executable." % SATsolver)
        sys.exit(1)

    # Initializing kwargs
    kwargs = {}

    # Verifying that user provided the graph's TXT
    if len(sys.argv) != 2:
        print("Usage: %s <input_file>" % sys.argv[0])
        sys.exit(1)

    # Parsing the graph's TXT
    input_file = sys.argv[1]
    with open(input_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
        k = int(lines[0])
        n, m = map(int, lines[1].split())
        edges = [tuple(map(int, line.split())) for line in lines[2:]]

    # Storing parsed data into kwargs
    kwargs['n'] = n
    kwargs['k'] = k
    kwargs['edges'] = edges

    # Generating the variables names and the clauses
    genVarNames(**kwargs)
    clauses = genClauses(**kwargs)

    # Helpers function provided by TAs (used to generate Z3 readable format)
    head = getDimacsHeader(clauses)
    cnf = toDimacsCnf(clauses)

    # Here we create a temporary cnf file for SATsolver
    fl = open("tmp_prob.cnf", "w")
    fl.write("\n".join([head, cnf]) + "\n")
    fl.close()

    # Run the SATsolver
    solverOutput = Popen([SATsolver + " tmp_prob.cnf"], stdout=PIPE, shell=True).communicate()[0]
    res = solverOutput.decode('utf-8')
    printResult(res)
    extractVertexColoring(res)