# 1.1.9
* Correct gene annotation used in regionplot when using build 37
* Fix of partly incorrect coloring of association points near the end or beginning of a chromosome.

# 1.1.8
* Missing shade added for chr22 on Manhattan plot when chrX is not included in the dataset.

# 1.1.7
* Placeholder added when there are no genes in the region on display
* Fix for the manhattan function to handle one chromosome datasets
* xaxis_label added as an input argument to the manhattan plot function

# 1.1.6
* Updated documentation

# 1.1.5
* Fix for chrX when displaying genes in regionplot

# 1.1.4
* Default significance threshold changed from 5e-9 to 5e-8
* Updated README

# 1.1.3
* Fix to make topr compatible with the next version of dplyr

# 1.1.2

* the scale parameter for annotate_with_vline was fixed
* the Manhattan plot was changed to the more conventional two-color display
* the argument theme_grey can be set to TRUE to make the Manhattan plot look like it does in older versions of topr

# 1.1.1

* Minor bug fix for plot colors and legends when plotting multiple datasets
* region argument added to the manhattan function

# 1.1.0

* Update to increase control over plot aesthetics and variant labeling
* minor bug fixes
* Includes the option to choose between genome builds 38 (GRCh38.p13) and 37 (GRCh37.p13).

# 1.0.0

* First official release of topr

# 0.0.1

* First version including locuszoom plot

# 0.0.0

* First version submitted to github
