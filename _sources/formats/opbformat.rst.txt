.. _opb_format:
  
=================
OPB format (.opb)
=================

The OPB file format is used to express pseudo-Boolean satisfaction and optimization models. 
These models may only contain :math:`0/1` Boolean variables. The format is defined by an optional objective function followed by a set of linear constraints.
Variables may be multiplied together in the objective function, but currently not in the constraints due to some restriction in the reader.
The objective function must start with the **min:** or **max:** keyword followed by **coef_1 varname_1_1 varname_1_2 ... coef2 varname_2_1 ...** and end with a **;**.
Linear constraints are composed in the same way, ended by a comparison operator (**<=**, **>=**, or **!=**) followed by the right-hand side coefficient and **;**.
Each coefficient must be an integer beginning with its sign (**+** or **-** with no extra space).
Comment lines start with a \*.

An example with a quadratic objective and 7 linear constraints is: ::

  max: +1 x1 x2 +2 x3 x4;
  +1 x2 +1 x1 >= 1;
  +1 x3 +1 x1 >= 1;
  +1 x4 +1 x1 >= 1;
  +1 x3 +1 x2 >= 1;
  +1 x4 +1 x2 >= 1;
  +1 x4 +1 x3 >= 1;
  +2 x1 +2 x2 +2 x3 +2 x4 <= 7;

Internally, all integer costs are multiplied by a power of ten depending on the -precision option. 
For problems with big integers, try to reduce the precision (*e.g.*, use option -precision 0).

