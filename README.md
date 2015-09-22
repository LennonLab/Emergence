#hydrobide

Individual based modeling of ecological and evolutionary dynamics in complex fluid and non-fluid systems. 

### Author: Kenneth J. Locey

#Purpose
The code in this repository was developed to simulate individual-based, probabilistic, eco-evolutionary, information intensive models. Let's break that down...

###Nomenclature
* **Individual-based:** Every particle is simulated. Particles can be organisms, resources, inert tracers, or dispersal barriers.

* **Probabilistic:** Changes to particles, populations, and communities are brought about through random sampling. This reflects the stochastic nature of individual responses to environmental conditions, competition, predator-prey interactions, dispersal, etc. Making all processes probabilistic also always greater freedom in how patterns, dynamics, and combinations of traits arise.

* **Eco-evolutionary:** Incorporates or will soon incorporate principles and processes from community ecology, biogeography, ecological niche and resource limitation theory, trade-offs from life history theory, population genetics, ecosystem science and nutrient stoichiometry, mass-balance, ecological and evolutionary neutral and nearly-neutral theory, and of course, the theory of evolution by natural selection.

* **Self-assembling**
The model is really a platform for modeling, because the user only sets the ranges of values for many parameters and processes the model can potentially include. But, when run, the model initiates with random combinations of parameters and processes, including:
	* Whether mutation and/or speciation occur
	* Whether reproduction is via binary fission or is sexual
	* If sexual, then whether reproduction is haploid or diploid
	* Whether immigration occurs
	* Whether flow is continuous or interrupted
	* Whether disturbance occurs
	* Whether the environment is fluid or static
	* Whether there are dispersal barriers and how large they are
	* Rate of flow 
	* Number and types of inflowing resources
	* Regional (outside) community structure
	* Aspects of fluctuating flow
	* Ecological disturbance

* **Data intensive** 
Information on all particles, groups of particles, and the entire systems is tracked, recorded, or quantified through time.

## Data that hydrobide tracks
**The following is tracked on every organism:**

* time in the system
* species ID and individual ID
* spatial location in three dimensions
* closest available resources and competitors
* levels of multiple endogenous resources (resource specific cell quotas).
* species and individual-level metabolic maintenance cost
* species and individual-level active dispersal ability
* maximum theoretical growth rate
* pedigree, i.e., direct lineage 
* resource-specific efficiency values
* what the individual is a predator of
* what the individual is prey to
* byproducts of metabolism

**The following is tracked on every resource particle:**

* time in the system
* particle ID
* location in three dimensions
* whether the particle is Nitrogen, Carbon, or Phosphorus
	* 	The modeling is increasingly incorporating stoichiometry
* what type of N, C, or P the resource is
* size of the particle

**The following is tracked on inert tracers:**

* time in the system
* location in three dimension
* particle ID

**Model variables:**

1. length
2. width
3. height
4. flow rate
5. immigration rate
6. number of resources
7. inflowing resource concentration
8. motion
9. metacommunity diversity
10. community max. growth rate
11. specific max. growth rate
12. community max. maintenance
13. specific max. maintenance
14. community active dispersal rate
15. specific active dispersal rate
16. community specific resource use efficiency
17. Subpopulation variation in:
	* Growth, maintenance, dispersal, and resource use
18. n sources of Carbon
19. m sources of Nitrogen
20. p sources of Phosphorus
21. Amplitude and frequency of fluctuating flow
22. Pulsing flow
23. Stochastic disturbance via decimation (randomly removing 1/10th of the community at random intervals)

21. *speciation rate
22. *interguild predation
23. *intraguild predation
24. *mutualism
25. *secondary resource production

***Still building**

## Using simulated data in R
Though coded in Python (and sometimes in Cython for speed), the output of hydrobide can be imported into Python and R environments as dataframes. R Markdown scripts are provided in the **analyses** folder.


