# Wahl_map
This repository contains a SageMath program associated with my doctoral thesis, "Surjectivity of the Wahl map," at the University of Kentucky from May 2023. This program evaluates the Wahl map on a given graph and generates random cubic graphs of girth at least 5. 

All graphs should be entered as Graph({vertex1:[neighbors1], vertex2:[neighbors2], vertex3:[neighbors3]}) or some other acceptable graph object according to SageMath.

To evaluate the Wahl map on a graph G, use build_matrix(G). To then check if the output is surjective, use is_surj(G). The output to is_surj() is a printed statement. Similarly, you can check if the Wahl map restricted to the vertices is surjective with is_V_surj(G). It will also print a statement.

In order to generate n non-isomorphic, random cubic graphs with m vertices of girth at least 5, run the command sample_random_girth5(n,m). If you want to make these graphs truely random for probabalistic purposed, add the entry noniso=False at the end of the input list. The output will be a list of graphs.

To test and store examples, run add_to_files(examples) where examples is a list of graphs.
