#Properties of simplex
* **Individual-based:** Every particle is simulated. Particles can be organisms with unique physiology, resource particles with unique stoichiometry, inert tracer particles to track fluid dynamics, or dispersal barriers that create spatial aggregation.

* **Probabilistic:** Changes to particles, populations, and communities occur through random sampling. This sampling can be biased or unbiased (like a biased or fair coin) and reflects the often unpredictable nature of ecological systems. Making all processes probabilistic also allows patterns, dynamics, and combinations of traits to arise more freely and for a larger number of possibilities to be explored.

* **Eco-evolutionary:** Incorporates principles and processes from community ecology, biogeography, ecological niche and resource limitation theory, trade-offs from life history theory, population genetics, ecosystem science and nutrient stoichiometry, mass-balance, ecological and evolutionary neutral and nearly-neutral theory, and of course, the theory of evolution by natural selection.

* **Self-assembling**
**simplex** is a platform for free-modeling where the user only sets the ranges of values for parameters and processes that a model will potentially include. When run, **simplex** initiates a model initiates with random combinations of parameters and processes, including:
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

## Data that simplex tracks
**The following is tracked on every organism:**

* time each particle spends in the system (aka residence time)
* species ID
* individual ID
* 2D spatial location
* closest available resources and competitors
* levels of 3 endogenous resources (resource specific cell quotas).
* species specific metabolic maintenance cost
* individual-level metabolic maintenance cost
* species-level maximum dispersal rate
* individual-level maximum dispersal rate
* species-specific maximum theoretical growth rate
* pedigree, i.e., direct lineage 
* species-specific resource use efficiency, for each resource
* products of individual-level metabolism

**The following is tracked on every resource particle:**

* time in the system
* particle ID
* 2D spatial location
* whether the particle is Nitrogen, Carbon, or Phosphorus
* what type of Nitrogen, Carbon, or Phosphorus the resource is
* size of the particle

**The following is tracked on inert tracers:**

* time in the system
* location in three dimension
* particle ID

**Model variables:**

* width
* height
* flow rate
* immigration rate
* number of resources
* inflowing resource concentration
* motion
* metacommunity diversity
* community max. growth rate
* specific max. growth rate
* community max. maintenance
* specific max. maintenance
* community active dispersal rate
* specific active dispersal rate
* community specific resource use efficiency
* Subpopulation variation in:
	* Growth, maintenance, dispersal, and resource use
* c sources of Carbon
* n sources of Nitrogen
* p sources of Phosphorus
* Amplitude and frequency of fluctuating flow
* Pulse in flow
* Stochastic disturbance via decimation (randomly removing 1/10th of the community at random intervals)
* speciation rate
* *interguild predation
* *intraguild predation
* *mutualism
* *secondary resource production

***Unbuilt or rebuilding**