
#simplex
<img src="https://upload.wikimedia.org/wikipedia/commons/e/e7/Tetrahedron-4-3D-balls.png" align="right" width="150" height="150" />

Because studying ecological complexity doesn't have to be complicated. 

###Purpose
**simplex** performs three tasks:

1. Assembles individual-based models from combinations of parameters and processes to simulate stochastic eco-evolutionary dynamics across many ecological and environmental conditions.

2. Quantifies, tracks, and records detailed information from genomes and physiological states of individuals to the aggregate properties of the entire ecosystem.

3. Provides quantitative tools to perform statistical analyses.

##Suggested software
**simplex** was developed on the free Enthought Canopy Python distribution (version 1.5.5) available here: https://store.enthought.com/

**simplex** implements unit testing using pytest version 2.8.1; see: http://pytest.org/latest/getting-started.html#getstarted

##Unit tests
Source code for unit testing, available in the **tools/unit_tests** directory.
The following unit testing is currently implemented:

* Several tests on 15 biodiversity metrics

To come:

* Tests on remaining biodiversity metrics
* Tests on functions that simulate ecological processes

##The ODD protocol
The ODD protocol is an accepted standard for describing individual-based models.
We descibe **simplex** as close as reasonable according to an ODD protocol.

Grimm, V. *et al*. (2006) A standard protocol for describing individual-based and agent-based models. Ecological Modeling. **198,** 115-126. *not OA*

## Using simulated data in R
Though coded in Python (and sometimes in Cython for speed), the output of **simplex** can be imported into Python and R environments as dataframes. R Markdown scripts are provided in the analyses folder.