.. _tuto_bmp:

======================
Block modeling problem
======================

.. include:: menu_backto.rst

Brief description
=================

This is a clustering problem, occuring in social network analysis. The problem is to divide a given graph G into k clusters such that the interactions between clusters can be summarized by a k*k 0/1 matrix M: if M[i,j]=1 then all the nodes in cluster i should be connected to all the nodes in cluster j in G, else if M[i,j]=0 then there should be no edge between the nodes in G. The goal is to find a k-clustering and the associated matrix M minimizing the number of erroneous edges.
`A Mattenet, I Davidson, S Nijssen, P Schaus. Generic Constraint-Based Block Modeling Using Constraint Programming. CP 2019, pp656-673, Stamford, CT, USA <https://www.jair.org/index.php/jair/article/download/12280/26656>`_.

CFN model
=========

We create N variables for every node of the graph with domain size k. We add k*k Boolean variables for representing M. For all triplets of two nodes u, v, and one matrix cell M[i,j], we have a ternary cost function which returns a cost of 1 if node u is assigned to cluster i, v to j, and M[i,j]=1 but (u,v) is not in G, or M[i,j]=0 and (u,v) in G. In order to break symmetries, we constrain the first k-1 node variables to be assigned to cluster index less than or equal to their index

Data
====

You can try a small example :download:`simple.mat <../../../web/TUTORIALS/simple.mat>` with optimum value equal to 0 for 3 clusters.

Perfect solution for the small example with k=3 (Mattenet et al, CP 2019)

.. image:: ../../../web/IMAGES/simple.png
   :height: 250px

More examples with 3 clusters (Stochastic Block Models `[Funke and Becker, Plos One 2019] <https://doi.org/10.1371/journal.pone.0215296>`_)

.. image:: ../../../web/IMAGES/SBM.png
   :height: 250px

See other examples, such as `PoliticalActor <https://www.ifip.com/Partitioning_Political_Actor.html>`_ and more, here :
:download:`100.mat <../../../web/TUTORIALS/100.mat>` |
:download:`150.mat <../../../web/TUTORIALS/150.mat>` |
:download:`200.mat <../../../web/TUTORIALS/200.mat>` |
:download:`30.mat <../../../web/TUTORIALS/30.mat>` |
:download:`50.mat <../../../web/TUTORIALS/50.mat>` |
:download:`hartford_drug.mat <../../../web/TUTORIALS/hartford_drug.mat>` |
:download:`kansas.mat <../../../web/TUTORIALS/kansas.mat>` |
:download:`politicalactor.mat <../../../web/TUTORIALS/politicalactor.mat>` |
:download:`sharpstone.mat <../../../web/TUTORIALS/sharpstone.mat>` |
:download:`transatlantic.mat <../../../web/TUTORIALS/transatlantic.mat>`.

Python model generator
======================

The following code using python3 interpreter will generate the corresponding cost function network (e.g. "python3 blockmodel.py simple.mat 3"). Download the AWK script :download:`sol2block.awk <../../../web/TUTORIALS/sol2block.awk>` to pretty print the results (e.g., "python3 blockmodel.py simple.mat 3 | toulbar2 --stdin=cfn -s=3 | awk -f ./sol2block.awk").

:download:`blockmodel.py<../../../web/TUTORIALS/blockmodel.py>`

.. literalinclude:: ../../../web/TUTORIALS/blockmodel.py

We improve the previous model by sorting node variables by decreasing out degree and removing the lower triangular matrix of M if the input graph is undirected (symmetric adjacency matrix).

:download:`blockmodel2.py<../../../web/TUTORIALS/blockmodel2.py>`

.. literalinclude:: ../../../web/TUTORIALS/blockmodel2.py
