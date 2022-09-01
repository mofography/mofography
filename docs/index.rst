.. mofography documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to mofography's documentation!
=========================================================




.. _what:

What is mofography?
=========================
add text here



.. _how:

How does it work?
=========================
We present a novel approach using advanced image processing techniques to
represent the cavities of a material as a *pore graph,* where the *nodes*
represent regions of arbitrary shape segmented around the local distance maxima,
and *edges* represent the connectivity of these regions. Regions connect through
one or more pore windows. For example, the pore graphs of ZIF-8 (cage-like) and
NU-1000 (micro and meso channels) is shown in Figure 1.



These nodes and edges can then be annotated with properties relevant to
adsorption or pore characterization, which are then aggregated to create
MOF-level features. Thus, the pore graph allows us to retain relevant
information, without having to work with large 3D images of MOF cavities. Using
pore graphs, we can understand and represent MOF pores comprehensively, compare
the cavities of two MOFs towards unsupervised ML (E.g., clustering of MOFs),
regress pore-influenced properties such as adsorption uptake (supervised ML) and
identify key pore characteristics necessary for achieving specific separations.


.. _why:

Why do we need it ?
==============================
Most of the state-of-the-art approaches to computational materials discovery of
MOFs via machine learning (ML), can be summarized in terms of three key
components. A way to sample the phase space of MOFs (databases, material
generators etc.), a way to make predictions (a regression model, a data lookup
table etc.), and a way to select/optimize materials to maximize the desired
property(s) (gradient descent, genetic algorithms, Bayes optimization etc.). How
we represent materials plays a key-role in all these steps. In adsorption, the
shape, size, and connectivity of the cavities often play a bigger role than the
chemistry of the individual structural building blocks (SBUs) that make up the
MOF. Arriving at numerical features that capture the 3D confining environment
of a nanoporous materials remains a challenging task till date.



.. _use:

why should you try it ?
=========================

-  One stop shop for analysis, features, computing and visualization for you materials
-  Pure python, so it plays well with your other favorite libraries like ASE, Pymatgen, RDKit etc.
-  Very easy to extend. If something you need is not there, you can add you soubroutines to it, without worrying about nested modules, inheritance or naything like that. In fact I would encourage you to do exactly that and fine-tune it to your needs. We hope that the code organization is sensible.




.. _installation:

Installation
=========================
add install insruct here


.. _citation:

Citation
=====================
add citation here


.. _core:

Core concepts
===================
These are short, hands-on tutorials explaining the core concepts in mofography.
In each example, there is a short introduction, followed by example code for calculation and visualization using mofography.
These concepts are organized in increasing order of complexity, and often rely on the earlier ones.

.. toctree::
   :caption: Core concepts
   :maxdepth: 1

   Tutorial 0: Structure Graph <irmof_example_23aug2020.ipynb>
   Tutorial 1: Distance Grid  <irmof_example_23aug2020.ipynb>
   Tutorial 2: Region segmentation  <irmof_example_23aug2020.ipynb>
   Tutorial 3: Applying Periodic Boundary Conditions <irmof_example_23aug2020.ipynb>
   Tutorial 4: Region Adjacency Graph (RAG) or Pore Graph <irmof_example_23aug2020.ipynb>
   Tutorial 5: Connections and Windows <irmof_example_23aug2020.ipynb>
   Tutorial 6: Paths, dimensionality and isomorphism  <irmof_example_23aug2020.ipynb>

.. _pop:

Analysis of popular nanoporous materials
===========================================
These are more complex, more realistic tutorials which leverage the basic tutorials to create nice workflows for 
accomplishing specific goals. Note that in most cases, the framework structure can be swapped out for yours. If it doesn't work as 
intended, perhaps there is scope for improvement on our side. Please raise a detailed issue, and someone will respond.

.. toctree::
   :caption: Complex examples
   :maxdepth: 1
   
   Analysing the cavities of IRMOF-1 <framework_visualization_example.ipynb>  
   Detecting micro, meso and intersections of NU-1000 <irmof_example_23aug2020.ipynb>
   Detecting the cation sites in LTA <irmof_example_23aug2020.ipynb>
   Detecting rings in LTA using structure graph <irmof_example_23aug2020.ipynb>
   Visualize framework in fun, custom ways <framework_visualization_example.ipynb>
   Decomposing an adsorption isotherm based on regions <irmof_example_23aug2020.ipynb>
   Estimating blocked pockets from not-perfect isotherm <irmof_example_23aug2020.ipynb>
   Grabbing data in bulk from mofdb <irmof_example_23aug2020.ipynb>
   Create and apply a mofography workflow using prefect and dask <irmof_example_23aug2020.ipynb>



Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
