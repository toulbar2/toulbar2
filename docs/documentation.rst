.. _documentation:

=============
Documentation
=============

User and reference manuals
==========================

  - :download:`User manual<../doc/userdoc.pdf>`
  - :download:`Reference manual<../doc/refman.pdf>`

Input formats
=============

The available file formats (possibly compressed by gzip or xz, e.g., .cfn.gz, .wcsp.xz) are :

  - Cost Function Network format (:download:`.cfn<../doc/CFNformat.pdf>` file extensions)
  - Weighted Constraint Satisfaction Problem (:download:`.wcsp<../doc/wcspformat.pdf>` file extension)
  - Probabilistic Graphical Model (`.uai <http://www.cs.huji.ac.il/project/PASCAL/fileFormat.php>`_ / .LG file extension ; the file format .LG is identical to .UAI except that we expect log-potentials)
  - Weigthed Partial Max-SAT (`.cnf/.wcnf <http://www.maxsat.udl.cat/08/index.php?disp=requirements>`_ file extension)
  - Quadratic Unconstrained Pseudo-Boolean Optimization (:download:`.qpbo<../doc/QpboFormat.txt>` file extension)
  - Pseudo-Boolean Optimization (`.opb <http://www.cril.univ-artois.fr/PB16/format.pdf>`_ file extension)

Some examples
=============

- A simple 2 variables maximization problem `maximization.cfn <https://github.com/toulbar2/toulbar2/raw/master/validation/default/maximization.cfn>`_ in JSON-compatible CFN format, with decimal positive and negative costs.                 
  
- Random binary cost function network :download:`example.wcsp<../web/EXAMPLES/example.wcsp.xz>`, with a specific variable ordering :download:`example.order<../web/EXAMPLES/example.order>`, a tree decomposition :download:`example.cov<../web/EXAMPLES/example.cov>`, and a cluster decomposition :download:`example.dec<../web/EXAMPLES/example.dec>`
  
- Latin square 4x4 with random costs on each variable :download:`latin4.wcsp<../web/EXAMPLES/latin4.wcsp.xz>`
  
- `Radio link frequency assignment CELAR <http://miat.inrae.fr/schiex/Doc/Export/CELAR.ps.gz>`_ instances :download:`scen06.wcsp<../web/EXAMPLES/scen06.wcsp.xz>`, :download:`scen06.cov<../web/EXAMPLES/scen06.cov>`, :download:`scen06.dec<../web/EXAMPLES/scen06.dec>`, :download:`scen07.wcsp<../web/EXAMPLES/scen07.wcsp.xz>`
  
- `Earth observation satellite management SPOT5 <https://link.springer.com/content/pdf/10.1023/A:1026488509554.pdf>`_ instances :download:`404.wcsp<../web/EXAMPLES/404.wcsp.xz>` and :download:`505.wcsp<../web/EXAMPLES/505.wcsp.xz>` with associated tree/cluster decompositions :download:`404.cov<../web/EXAMPLES/404.cov>`, :download:`505.cov<../web/EXAMPLES/505.cov>`, :download:`404.dec<../web/EXAMPLES/404.dec>`, :download:`505.dec<../web/EXAMPLES/505.dec>`
  
- Linkage analysis instance :download:`pedigree9.uai<../web/EXAMPLES/pedigree9.uai.xz>`
  
- Computer vision superpixel-based image segmentation instance :download:`GeomSurf-7-gm256.uai<../web/EXAMPLES/GeomSurf-7-gm256.uai.xz>`
  
- `Protein folding <http://miat.inrae.fr/degivry/Schiex14a.pdf>`_ instance :download:`1CM1.uai<../web/EXAMPLES/1CM1.uai.xz>`
  
- Max-clique DIMACS instance :download:`brock200_4.clq.wcnf<../web/EXAMPLES/brock200_4.clq.wcnf.xz>`
  
- Graph 6-coloring instance :download:`GEOM40_6.wcsp<../web/EXAMPLES/GEOM40_6.wcsp.xz>`

Many more instances available `here <http://genoweb.toulouse.inra.fr/~degivry/evalgm>`_ and  `there <https://forgemia.inra.fr/thomas.schiex/cost-function-library>`_.

Command line arguments
======================

See *'Available options'* below :

.. literalinclude:: ../HELP

