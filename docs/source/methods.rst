.. methods::

Methods
=======

Enrichment of histone marks
---------------------------

The number of MNase reads per nucleosome is a confounding variable for the number of histone reads on the same nucleosome. Thus, to minimize spurious enrichment signals, for each nucleosome we normalized the number of histone reads, :math:`x_j` by the expected number of histone reads due to the number of MNase reads alone:

.. math::

   r_j = \frac{ x_j } { E(X|n_j)}

Here, :math:`E(X|n_j)` correspond to expected number of histone reads, :math:`X`, given a number :math:`n_j`, of MNase reads. 

.. math::

   E(X|n_j) = \frac{1}{||J(n_j)||} \sum_{j \in J(n_j)} x_j
