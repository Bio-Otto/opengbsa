Generalized Born Models
=======================

The choice of Generalized Born (GB) model determines how the effective Born radii are calculated and how the polar solvation energy is estimated. OpenGBSA provides access to several standard GB models implemented in OpenMM.

Model Descriptions
------------------

OBC1 (Onufriev-Bashford-Case I)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A fundamental improvement over the original GB model, the OBC models rescale the effective Born radii to better match Poisson-Boltzmann results.

*   **Key Feature**: Rescales Born radii to account for interstitial dielectrics.
*   **Citation**: Onufriev, A., Bashford, D., & Case, D. A. (2000). Modification of the generalized Born model suitable for macromolecules. *The Journal of Physical Chemistry B*, 104(6), 1541-1548.
*   **OpenMM Tag**: ``OBC1``

OBC2 (Onufriev-Bashford-Case II)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Recommended Default**. An optimization of OBC1 with different scaling parameters (:math:`\alpha, \beta, \gamma`). It is widely regarded as one of the most accurate analytical GB models for protein-ligand binding affinities.

*   **Key Feature**: Optimized parameters (:math:`\alpha=1.0, \beta=0.8, \gamma=4.85`) that reduce errors for larger proteins compared to OBC1.
*   **Citation**: Onufriev, A., Bashford, D., & Case, D. A. (2004). Exploring protein native states and large-scale conformational changes with a modified generalized born model. *Proteins: Structure, Function, and Bioinformatics*, 55(2), 383-394.
*   **OpenMM Tag**: ``OBC2``

GBn (GB-Neck)
~~~~~~~~~~~~~
Addressed a specific deficiency in earlier generic GB models regarding the "neck" region between two heteroatoms.

*   **Key Feature**: Improved hydration free energies for molecular shapes with crevices or "necks".
*   **Citation**: Mongan, J., Simmerling, C., & McCammon, J. A. (2007). Generalized Born model with a simple, robust molecular volume correction. *Journal of Chemical Theory and Computation*, 3(1), 156-169.
*   **OpenMM Tag**: ``GBn``

GBn2
~~~~
A further refinement of the GBn model, parameterized with a larger dataset.

*   **Citation**: Nguyen, H., Roe, D. R., & Simmerling, C. (2013). Improved generalized Born solvent model parameters for protein simulations. *Journal of Chemical Theory and Computation*, 9(4), 2020-2034.
*   **OpenMM Tag**: ``GBn2``

HCT (Hawkins-Cramer-Truhlar)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
One of the earliest pairwise descreening approximations. It is computationally efficient but often considered less accurate than OBC models for complex biomolecular systems.

*   **Citation**: Hawkins, G. D., Cramer, C. J., & Truhlar, D. G. (1995). Pairwise solute descreening of solute charges from a dielectric medium. *Chemical Physics Letters*, 246(1-2), 122-129.
*   **OpenMM Tag**: ``HCT``

Salt Screening
--------------

The presence of mobile ions in the solvent is modeled using Debye-Hückel screening. The screening effect reduces the range of electrostatic interactions.

The screening parameter :math:`\kappa` (inverse Debye length) is calculated from the salt concentration :math:`C` (in Molar):

.. math::

   \kappa \approx 0.304 \sqrt{C} \quad \text{(at 298 K for monovalent salt)}

In OpenGBSA, this is controlled by the ``salt_concentration`` parameter in the configuration file.
