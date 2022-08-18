.. _tuto_bcp:

========================
Board coloration problem
========================

.. include:: menu_backto.rst

Brief description
=================

Given a rectangular board with dimension n*m, the goal is to color the cells such that any inner rectangle included inside the board doesn't have all its corners colored with the same color.
The goal is to minimize the number of colors used.

For examples, this is not a valid solution of the 3*4 problem, because the red and blue rectangles have both their 4 corners having the same color:

.. image:: ../../../web/IMAGES/tchess.png

On the contrary the following coloration is a valid solution of the 3*4 problem because every inner rectangle inside the board does not have a unique color for its corners:

.. image:: ../../../web/IMAGES/fchess.png

CFN basic model
===============

We create n*m variables, one for each square of the board, with domain size equal to n*m representing all the possible colors. We also create one variable for the number of colors.

We create hard quaterny constraints for every rectangle inside the board with cost equal to 0 if the 4 variables have different values and a forbidden cost if not.

We then create hard binary constraints between the variable of the number of colors for each cell to fix the variable for the number of colors as an upper bound.

Then we create a soft constraint on the number of colors to minimize it.


Python model
============

The following code using pytoulbar2 library solves the board coloration problem with the first two arguments being the dimension n and m of the board (e.g. "python3 boardcoloration.py 3 4").

:download:`boardcoloration.py<../../../web/TUTORIALS/boardcoloration.py>`

.. literalinclude:: ../../../web/TUTORIALS/boardcoloration.py

