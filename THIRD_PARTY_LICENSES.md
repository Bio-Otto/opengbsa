# Third-Party Licenses

OpenGBSA itself is MIT licensed (see [`LICENSE`](LICENSE)). This file
documents third-party components vendored or bundled with the conda
package that carry a different license.

## TprParser (GPLv3)

Native `.tpr` (GROMACS binary topology) parsing is provided by an
optional, vendored [`TprParser`](https://pypi.org/project/TprParser/)
wheel (`conda-recipe/wheels/`), installed by the conda package's
post-link script only when a matching wheel exists for the target
platform and Python version.

`TprParser` is licensed under the
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html)
by its author, Yujie Liu. It is distributed as a separate, optional
component -- not statically or dynamically linked into `opengbsa`'s
own code -- and its absence does not disable any other feature.
Without it, native `.tpr` input is unavailable, but all other input
formats (Amber `.prmtop`, GROMACS `.top`, PDB) and every analysis
feature work unaffected.

If you redistribute a build of OpenGBSA that bundles this wheel, you
are also redistributing GPLv3-licensed code and should comply with its
terms for that component.
