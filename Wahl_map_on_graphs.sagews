︠00662b9f-1bd9-48e5-a702-be67461df9d1r︠
'''
Code by Angela Hanson. Last updated 9/12/22.
'''
from __future__ import absolute_import
from sympy import *
import sys
from sage.all import *
from sage.graphs import graph
from sage.graphs.generic_graph import *


###### This section of the code evaluates the Wahl map on a given graph. ######
def pair_up_cycles(cycleSet):
    '''This will produce pairs of cycles from the set of cycles in the order in which they
    appear in the set. The pairs are returned as a list of pairs where each pair is a list
    of two cycles. Cycles are represented as a list of integers which correspond to vertex
    labels.'''
    pairs=list()

    for i in range(len(cycleSet)):
        for j in range(i+1,(len(cycleSet))):
            pairs.append([cycleSet[i],cycleSet[j]])

    return(pairs)


def cycle_orientation(v, neighbor_pt, cycle):
    '''This will take a vertex v and one of its neighbors (neighbor_pt) and calculates
    whether the cycle runs into vertex v along that neighbor or out of the vertex along that
    neighbor. The output of this function is 1 (into) or -1 (out of).'''

    if cycle.index(v) == 0 and cycle.index(neighbor_pt) == len(cycle)-1:
        '''This accounts for if the neighbor is the last point in the cycle list and v is
        the first in the cycle list, so neighbor_pt comes before v.'''
        return(1)
    elif cycle.index(v) == len(cycle)-1 and cycle.index(neighbor_pt) == 0:
        '''This accounts for if the neighbor is the first point in the cycle list and v is
        the last in the cycle list, so neighbor_pt comes after v.'''
        return(-1)
    elif cycle.index(v) == cycle.index(neighbor_pt)-1:
        '''Here v and its neighbor are not at the end or beginning of the cycle list, and v
        comes before neighbor_pt.'''
        return(-1)
    elif cycle.index(v)-1 == cycle.index(neighbor_pt):
        '''Here v and its neighbor are not at the end or beginning of the cycle list, and v
        comes after neighbor_pt.'''
        return(1)
    else:
        '''This is a catch-all in case neighbor_pt and v are not actually both in this cycle
        consecutively.'''
        return(0)


def edge_orientation(v, neighbor_pt):
    '''This will take a vertex v and one of its neighbors (neighbor_pt) and calculates
    whether the edge orientation runs into vertex v along the neighbor or out of the vertex
    along the neighbor. The output of this function is 1 (into) or -1 (out of). The edge
    orientation is, by default, directed from lowest vertex label to highest. Since vertices
    are labelled by the natural numbers without repetition, this orientation is
    well-defined.'''

    if v < neighbor_pt:
        '''Since v is smaller than neighbor_pt, it comes before, so the edge between them is
        oriented out of v.'''
        return(-1)
    elif v > neighbor_pt:
        '''Since v is larger than neighbor_pt, it comes after, so the edge between them is
        oriented into v.'''
        return(1)
    else:
        print("Error with finding neighbor edge value. Check vertex labelling.")


def w_vertex(v, cycle1, cycle2, theGraph):
    '''This takes a vertex v and pair of cycles, cycle1 and cycle2, and calculates the Wahl
    map of the wedge of those two cycles on v. This involves calculating the determinant of
    a 2x2 matrix where the first row corresponds to cycle1 and the second to cycle2 and then
    the columns correspond to two of the neighboring vertices to v. Since these graphs are
    cubic, there are exactly three neighbors to each vertex. We will select two neighbors
    according to the permutation phi which is (123), where 1, 2, and 3 correspond to the
    first, second, and third neighbors of v in ascending order. This does not effect the
    rank of the Wahl map because rank is independent of vertex orienation. It returns the
    determinant value times the edge orientation values for each neighboring edge used in
    the determinant calculation. This product can equal -1, 0, or 1.'''
    neighborlist = theGraph.neighbors(v)
    n1 = min(neighborlist[0],neighborlist[1],neighborlist[2])

    '''Here we find the second neighbor based on the vertex orientation phi. Since phi has
    three elements for the three neighbors, we want the next neighbor to be next in the list
    modulo 3.'''
    index1 = neighborlist.index(n1)
    index2 = Mod(index1+1,3)
    n2 = neighborlist[index2]

    '''This starts the determinant calculation where aij is the ith row and jth column
    entry. The value comes from the cycle orientation from v to its neighbor. So a12 is the
    orientation value from v to n2 in cycle1. If the relavant neighbor is not in the
    relavant cycle, that matrix entry is 0.'''
    if (v in cycle1) and (v in cycle2):
        if (n1 in cycle1):
            a11 = cycle_orientation(v,n1,cycle1)
        else:
            a11 = 0
        if (n2 in cycle1):
            a12 = cycle_orientation(v,n2,cycle1)
        else:
            a12 = 0
        if (n1 in cycle2):
            a21 = cycle_orientation(v,n1,cycle2)
        else:
            a21 = 0
        if (n2 in cycle2):
            a22 = cycle_orientation(v,n2,cycle2)
        else:
            a22 = 0

        '''This is where we calculate the edge orientation value of 1 or -1 for towards or
        away from v, respectively.'''
        delta1 = edge_orientation(v, n1)
        delta2 = edge_orientation(v, n2)
        det = delta1*delta2*(a11*a22 - a12*a21)
    else:
        '''If v is in neither cycle1 nor cycle2, then the matrix determinant is
        trivially 0.'''
        det = 0

    return(det)



def w_edge(edge, cycle1, cycle2, theGraph):
    '''This takes an edge=(v1 and v2) and pair of cycles, cycle1 and cycle2, and calculates
    the Wahl map of the wedge of those two cycles on that edge. This involves calculating
    the determinant of a 2x2 matrix where the first row corresponds to cycle1 and the second
    to cycle2 and then the first column corresponds to a neighboring edge to v1 of edge and
    the second column to a neighboring edge of v2 of edge. Since these graphs are cubic,
    there are two neighbors to each vertex. We will select the neighbor with minimal vertex
    label for consistency. This does not effect the rank of the Wahl map since it is the
    same as the choice of vertex orientation. It returns the determinant value times the
    edge orientation values for each neighboring edge used in the determinant calculation.
    This product can equal -2, -1, 0, 1, or 2.'''
    endpoints = list(edge)
    v1 = endpoints[0]
    v2 = endpoints[1]

    '''Here we find the neighbors of each endpoint, v1 and v2, including the other endpoint
    of the given edge. Since these graphs are trivalent, there are three vertex neighbors.'''
    neighborlist1 = theGraph.neighbors(v1)
    neighborlist2 = theGraph.neighbors(v2)

    '''This chooses the neighbors n1 for v1 and n2 for v2 according to the permutation phi.
    To make sure the permutation does not return the other endpoint of the edge we are
    referencing, we make the other endpoint the starting vertex in the permutation and
    choose the neighbor that is next in the orientation.'''
    index1 = neighborlist1.index(v2)
    index2 = neighborlist2.index(v1)
    '''Here we find the next neighbor of v1, after v2, based on the vertex orientation phi.
    Since phi has three elements for the three neighbors, we want the next neighbor to be
    next in the list modulo 3.'''
    index3 = Mod(index1+1,3)
    n1 = neighborlist1[index3]
    '''Here we find the next neighbor of v2, after v1, based on the vertex orientation phi.
    Since phi has three elements for the three neighbors, we want the next neighbor to be
    next in the list modulo 3.'''
    index4 = Mod(index2+1,3)
    n2 = neighborlist2[index4]

    '''This starts the determinant calculation where aij is the ith row and jth column
    entry. The value comes from the cycle orientation from v1 to n1 and from v2 to n2. So
    a12 is the orientation value from v2 to n2 in cycle1. If the relavant neighbor is not
    in the relavant cycle, that matrix entry is 0.'''
    if (v1 in cycle1) and (n1 in cycle1):
        a11 = cycle_orientation(v1,n1,cycle1)
    else:
        a11 = 0
    if (v2 in cycle1) and (n2 in cycle1):
        a12 = cycle_orientation(v2,n2,cycle1)
    else:
        a12 = 0
    if (v1 in cycle2) and (n1 in cycle2):
        a21 = cycle_orientation(v1,n1,cycle2)
    else:
        a21 = 0
    if (v2 in cycle2) and (n2 in cycle2):
        a22 = cycle_orientation(v2,n2,cycle2)
    else:
        a22 = 0

    '''This is where we calculate the edge orientation value of 1 or -1 for towards or away
    from v1 and v2, respectively.'''
    delta1 = edge_orientation(v1, n1)
    delta2 = edge_orientation(v2, n2)

    return(delta1*delta2*(a11*a22 - a12*a21))


def generate_row(cycle1, cycle2, theGraph):
    '''This function takes in two cycles, cycle1 and cycle2, and returns a vector which
    contains the Wahl map values on each vertex in ascending order and then on each edge in
    lexicographic order.'''
    vect=list()

    for v in theGraph.vertices():
        entry = w_vertex(v, cycle1,cycle2, theGraph)
        vect.append(entry)
    for e in theGraph.edges():
        entry = w_edge(e, cycle1, cycle2, theGraph)
        vect.append(entry)

    return(vector(ZZ,vect))


def build_matrix(G):
    '''Here the the graph G is used to calculate the matrix of Wahl map value vectors for
    all cycle wedge pairs.'''

    '''The list of cycles we use is the cycle basis determined by SageMath.'''
    pairs=pair_up_cycles(G.cycle_basis())

    '''We need an empty matrix of the proper dimension to fill with Wahl map values. This
    matrix has as many rows as cycle pairs and as many columns as vertices and edges.'''
    edgeMatrix=matrix(ZZ,len(pairs),len(G.vertices())+len(G.edges()))

    '''Here we find the vector of Wahl map values for each cycle wedge pair and insert them
    in the matrix as rows.'''
    for i in range(len(pairs)):
        edgeMatrix.set_row(i,generate_row(pairs[i][0],pairs[i][1],G))

    return(edgeMatrix)


def is_surj(G):
    '''This takes in a graph and prints whether or not the Wahl map is surjective on that
    graph.'''
    M=build_matrix(G)

    if M.rank() == M.ncols():
        print("The Wahl map is surjective on this graph.")
    else:
        print("The Wahl map is NOT surjective on this graph.")
        print(M.rank(), M.ncols())


def is_V_surj(G):
    '''This takes in a graph and prints whether or not the vertex map V is surjective on
    that graph. Note that the Wahl map on a graph equals V+E where V is the map restricted
    to acting on the vertices, and E is the map restricted to acting on the edges.'''
    M=build_matrix(G)

    '''Here we build a list of indices for the rows and columns of the image of V.'''
    r = list()
    for i in range(M.nrows()):
        r.append(i)
    c = list()
    for j in range(len(G.vertices())):
        c.append(j)

    '''This builds a submatrix of the Wahl matrix of values to only keep the columns
    corresponding to vertex values.'''
    vM = M.matrix_from_rows_and_columns(r, c)

    if vM.rank() == vM.ncols():
        print("The vertex map is surjective on this graph.")
    else:
        print("The vertex map is NOT surjective on this graph.")
        print(vM.rank(), vM.ncols())


###### This section of the code will generate random cubic graphs of minimum girth 5. ######
def generate(numVertices,girth=5):
    '''This function intakes the number of vertices desired and returns a cubic graph of
    with that many vertice and girth at least 5.'''
    n = numVertices
    '''We need an empty graph to add edges and vertices to.'''
    G = Graph()

    '''This double checks that there are at least the girth number of vertices.'''
    if n < girth:
        print("You have not given enough vertices to build a graph of this girth.")
    else:
        '''Here we build the vertex set labeled 0 to n-1 and form the primary cycle 0 up to
        n-1 and back to 0. This makes all vertices valence 2 so we add exactly more edge per
        vertex in the next step to be cubic.'''
        for v in range(n-1):
            G.add_edge(v,v+1)
        G.add_edge(n-1,0)

        remaining_vertices = [i for i in range(0,n)]
        for v1 in remaining_vertices:
            '''This stops any infinite loops looking for girth big enough when it is not
            possible. We do this by working from a list specific to v1.'''
            possible_neighbors = remaining_vertices.copy()

            '''We remove vertices adjacent to v1 in the original cycle so that we remain
            simple. We also remove vertices within girth-1 of v1 to avoid obvious girth
            violations.'''
            for j in range(girth):
                if Mod(v1+j,n) in possible_neighbors:
                    possible_neighbors.remove(Mod(v1+j,n))
                if Mod(v1-j,n) in possible_neighbors:
                    possible_neighbors.remove(Mod(v1-j,n))

            '''We want to keep looking for a neighbor for v1 until it is degree 3 or we run
            out of vertices to pair with it.'''
            while G.degree(v1) < 3 and len(possible_neighbors) > 0:

                '''Now we pick the second vertex to form an edge with.'''
                v2_index = sage.misc.prandom.randint(0,len(possible_neighbors)-1)
                v2 = possible_neighbors[v2_index]

                G.add_edge(v1,v2)
                if G.girth() < girth or G.degree(v2) > 3:
                    '''If girth is too small, we remove the edge and drop v2 from the
                    potential neighbors of v1.'''
                    G.delete_edge(v1,v2)
                    possible_neighbors.remove(v2)
                else:
                    '''If girth is high enough, we remove v2 from the vertex list so that
                    valence of v2 stays at 3.'''
                    remaining_vertices.remove(v2)

        if G.is_regular() and G.girth()==5:
            return(G)
        else:
            return(Graph())



def sample_random_girth5(n,vcount,noniso=True,check_planarity=False,check_triconnected=False):
    '''Here we choose a number n of random cubic graphs we want to test, and it returns a
    list of n cubic graphs. If noniso=True, then these graphs will be non-isomorphic. Otherwise,
    the graphs will be truly random, with possible repetition. We can specify the number of
    vertices or have a number selected at random. We also have the option to check for non-planar,
    3-vertex-connected graphs or not. By default this function does not make those checks.'''
    sample = list()
    graphset = list()
    total_number = 0

    while total_number < n:
        G = generate(vcount)
        new = True
        if G != Graph():

            '''This skips graphs that are the same up to isomorphism as ones we seen before
            in a single test run'''
            if noniso:
                for H in sample:
                    if H.is_isomorphic(G):
                        new = False

            if new == True:
                sample.append(G)

                if check_planarity and check_triconnected:
                    if G.is_planar() == False and G.is_triconnected():
                         '''These are 3-vertex-connected and non-planar graphs which we
                         generally want. They are also of girth at least 5 which seems to be
                         necessary for surjectivity.'''
                         graphset.append(G)
                         total_number += 1
                elif check_planarity:
                    if G.is_planar() == False:
                        '''These are girth at least 5, non-planar graphs, but they may not
                        be 3-connected.'''
                        graphset.append(G)
                        total_number += 1
                elif check_triconnected:
                    if G.is_triconnected():
                        '''These are girth at least 5, 3-connected graphs, but they may not
                        be non-planar.'''
                        graphset.append(G)
                        total_number += 1
                else:
                    '''These are graphs of girth at least 5.'''
                    graphset.append(G)
                    total_number += 1

    return(graphset)

def add_to_file(examples):
    '''This function creates a file and stores the graphs in the list of examples. It also
    records which of them the Wahl map is surjective on.'''
    lines = list()

    '''Here we build a list of which of the graphs the Wahl map is surjective on.'''
    isSet = list()
    for g in examples:
        M=build_matrix(g)
        if M.rank() == M.ncols():
            isSet.append(examples.index(g))

    print(str(len(isSet))," out of ",str(len(examples)))
    '''This puts how many graphs the Wahl map is surjective on, how many are in the example
    list, and the indexes of the graphs the Wahl map is surjective on.'''
    lines.append(str(len(isSet)))
    lines.append(str(len(examples)))
    lines.append(str(isSet))

    for ex in examples:
        '''This is where the examples are recorded as an index, its girth, and its list of
        edges.'''
        lines.append(str(examples.index(ex)))
        lines.append(str(ex.girth()))
        lines.append(str(ex.edges()))

    with open('RandomExamples.sagews', 'w') as f:
        '''Here is where the file is created and written in with the information above.'''
        f.write('\n'.join(lines))


add_to_file(sample_random_girth5(50,275))
print("DONE")
︡c11cf849-aede-417d-9822-6ec8e2abdac2︡{"stdout":"'\\nCode by Angela Hanson. Last updated 9/12/22.\\n'"}︡{"stdout":"\n"}











