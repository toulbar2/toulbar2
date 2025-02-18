import pytoulbar2
import numpy as np
import sys

def read_dimacs_graph(filename,k):
    """Reads a graph file in the DIMACS format and returns a k-coloring CFN.

    Args:
        filename: The name of the DIMACS file to read.
        k: number of colors

    Returns:
        A CFN representing the k-coloring instance.
    """

    diff = np.eye(k).flatten()

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('c'):  # Comment line
                continue
            elif line.startswith('p'):  # Problem line
                _, _, num_vertices, num_edges = line.split()
                num_vertices = int(num_vertices)
                num_edges = int(num_edges)
                for i in range(num_vertices):
                    CFN.AddVariable(f'x{i+1}',list(range(k)))
            elif line.startswith('e'):  # Edge line
                _, vertex1, vertex2 = line.split()
                vertex1 = int(vertex1)
                vertex2 = int(vertex2)
                print(vertex1,vertex2)
                CFN.AddFunction([f'x{vertex1}',f'x{vertex2}'],diff)

    return CFN

CFN = pytoulbar2.CFN()
if len(sys.argv) != 3:
    throw("Needs a filename and number of colours as argument!")

k = int(sys.argv[2])
filename = sys.argv[1].strip()

CFN = read_dimacs_graph(filename,k)

CFN.SetUB(1)
print(CFN.Solve())
