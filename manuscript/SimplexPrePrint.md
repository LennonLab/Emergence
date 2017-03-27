---
title: Simplex: A modeling platform for the simultaneous emergence of ecological relationships

author(s):
    - Kenneth J. Locey, Jay. T. Lennon
date: \today{}
geometry: margin=2cm
output: pdf_document
header-includes:
    - \usepackage{setspace}
    - \doublespacing
    - \usepackage{lineno}
    - \linenumbers
---

# Introduction

**P1**
Ecology is the study of interactions and relationships between organisms and their environment.
This broad definition encompasses a multitude of subfields, each with its own principles, processes, and theories.
From the microscopic scales of molecular ecology and microbial ecology to the continental and global scales of macroecology, ecology spans all scales of space and abundance on Earth.
Likewise, from the geological scales of paleoecology to the contemporary scales of urban ecology, ecology spans all scales of time within which organisms live and evolve.
While ecologists have sometimes sought to unify the observations of their subfields under common principles, ecologists have yet to ask whether the combination of their paradigms allows otherwise disparate ecological patterns to emerge simultaneously.

**P2**

Ecologists often argue the importance and perspectives of their own subfields and methodologies, e.g., controlled small-scale experiments vs. powerful global-scale statistical patterns.


**P3**



**P4**



**P5**



Modeling is an elementary approach to understanding ecological systems, the influence of ecological processes, and the predictability of ecological patterns and dynamics.
Before modern computing, ecological models were almost exclusively equation-based representations of highly simplified systems (Black and McKane 2012). Since the advent of personal computing, ecological models have been increasingly constructed to handle greater complexity and to explicitly simulate ecological processes.

In ecology, simulation-based models are often used to examine analytically challenging or intractable scenarios.
Examples are markov models that simulate stochastic demographic changes and early ecological null models, both of which operate on matrices (Gotelli and Entsminger 2001, Hubbell 2001).
Ecologists have also simulated growth and interactions among individuals for over two decades with individual-based models (IBMs) (DeAngelis and Gross 1992, Rosindell et al. 2015).
IBMs explicitly encode rules of how individuals change and interact. 
Once the rules are encoded, population to ecosystem-level dynamics can emerge as an IBM simulates over time and spatially explicit environments.

A body of literature including original research, comprehensive reviews, and textbooks reveal the use, advantages, and challenges of ecological IBMs (Grimm and Railsback 2005).
IBMs can provide degrees of ecological realism, individual variability, and spatial heterogenetiy that are unattainable with other models. IBMs also offer the potential for realistic and unanticipated ecological dynamics and patterns to emerge from individual-level interactions.
However, IBMs offer challenges that include greater computational complexity, the difficulty of explicitly encoding ecological theory, and the use of ecological problem solving (Matthews et al. 2007, Grimm and Railsback 2005).

Despite their frequent use, the power of ecological IBMs has yet been leveraged to the greatest advantage.
Few, if any, include fluid dynamics or combine fluid dynamics with growth and active dispersal.
While growth, sensing, and decision making are often modeled, few IBMs integrate physiology, evolution, community ecology, biogeography, and sampling theory. 
Ecological IBMs are more often constructed to model a specific system than to synthesize general theories and mechanisms of ecology and evolution (but see Rosindell et al. 2015). 
Yet, IBMs can allow researchers to simulate and track information from the level of genomes and internal physiology to the stoichiometry of resource particles and distributions of abundance among species and across space.

Here, we present an IBM platform for exploring the simultaneous emergence of relationships and patterns under combined ecological paradigms. 
We refer to the platform as Simplex, referring to the emergence of complex systems from simple conditions.
A 'simplex' is also a geometric concept and a generalized notion of triangle abstracted to n-dimensions.
This geometric definition is also fitting because the Simplex platform is based on the general importance of three phenomena. These are the multiplicative interaction of stochastic variables (i.e., lognormal dynamics), the ubiquitous influence of energetic constraints, and the importance of ecological selection.
Simplex also accomplishes three primary tasks. 
First, Simplex simulates a broad range of ecological conditions.
Second, Simplex records information from individuals to the ecosystem level.
Third, Simplex provides computing code for analyzing simulation output.
Below, we provide detailed explanation of how Simplex works, the data it quantifies and tracks, the theories and principles Simplex integrates, and the analyses that can be conducted using the code we provide.


# Methods

## Platform description
Simplex source code is written in Python, a versatile high-level programming language that has many scientific and plotting libraries.
Here, we describe the concepts and capabilities of Simplex largely according to the ODD protocol (Overview, Design concepts, Details), which is standard for describing IBMs (Grimm et al. 2006). 
Because Simplex is a platform for studying the simultaneous emergence of ecological patterns and is provided with modifiable source code and analysis files, Simplex is not simply a single IBM for use with a limited set of hypotheses.
Detailed descriptions of Simplex source files, functions, and analysis code can be found on the public GitHub source code repository (https://github.com/LennonLab/simplex).

### Purpose
Simplex accomplishes three primary objectives. 
First, Simplex combines aspects of ecological paradigms of metabolic scaling, macroecology, neutral theory, resource limitation theory, ecological energetics, and chemostat theory. 
Second, Simplex reveals how realistic trait syndromes and patterns of metabolic scaling, biodiversity, and community ecology can simultaneously emerge from ecological selection on initially random conditions.
Third, Simplex reveals the general importance of multiplicative interactions of random variables (i.e., lognormal dynamics), energetic constraints, and ecological selection.
Simplex also accomplishes three proximate objectives. 
First, Simplex assembles and runs IBMs from random combinations of system variables and species traits.
Second, Simplex stores the outputs of many iterations including animations of models. 
Third, Simplex provides R and python code for analyzing simulation data.

### Entities & their state variables
**Individual organisms**
-- Simplex simulates life history processes of growth, dispersal, reproduction, and basal respiration at the individual level. 
Individuals are distinguished by collections of elements within dictionaries, i.e., data objects that hold key-value pairs. 
For example, the dictionary holding information on individual organisms is structured as follows:

	iDict = {'ind1' : 'sp' = 0, 
					  'x'  = 1.2, 
					  'y'  = 3.5, 
					  'sz' = 813.2, 
					  'q'  = 1, 
					  'st' = 'a'; 
			 'ind2' : ...}

where 'ind1' is the key and the variables following the colon are the values for the species ID, x-coordinate, y-coordinate, body size, amount of endogenous resources, and metabolic state (active or dormant).
Individuals undergo changes when randomly sampled from the dictionary.

**Species**
-- Each species is characterized by the individuals that share a common set of traits, such as maximum growth rate, metabolic maintenance cost. 
Species information is stored in dictionaries, again as key-value pairs. 

	spDict = {'1' : 'gr' = 0.8, 
					'di'  = 0.5, 
					'rp'  = 0.3, 
					'mt' = 0.2, 
					'mf'  = 0.1, 
					'ef' = [0.1, 0.2, 0.3]; 
			  '2' : ...}

where the species with ID of '1' has an intrinsic growth rate of 0.8, an active dispersal rate of 0.5, a 0.3 probability of randomly resuscitating from a metabolically dormant state, a basal mass specific metabolic rate of 0.2, and resource specific use efficiencies of 0.1, 0.2, and 0.3.

When sampling individuals, the information about their species is gained by accessing the species dictionary.

**Resource particles**
-- Simplex simulates the movement and consumption of individual resource particles. These particles can vary over several orders of magnitude in size and belong to three resource types. As with individual organisms and species, information about individual resource particles is stored in dictionaries.  

	rDict = {'1' : t = 1, 
				   'v'  = 0.5, 
				   'x'  = 0.3, 
				   'y' = 0.2, 
			 '2' : ...}
	

where 't' is the resource type, 'v' is the size of the particle, and 'x' and 'y' are the x and y coordinates.

### System level state variables
Each Simplex model begins with random choices for the values of:

* width in arbitrary units
* height in arbitrary units
* flow through rate in units of distance moved per time step by the environmental matrix; a minimum of 0.

### Spatial and temporal scales

**Spatial extent**
-- The environment of Simplex is square and two dimensional, and can vary along each axis from 1 to any number of arbitrary units. 
All particles move in decimal units the limit of which is determined by Python's decimal precision. 
This means that individual particles can occupy practically infinite locations.

**Temporal extent**
-- Simplex models can run for any number of time steps and record data at any number of time steps.

### Process overview and scheduling
**Model assembly**
-- The Simplex user runs the main python program (i.e., main.py). 
The main program chooses random values for system-level state variables including whether disturbance, immigration, speciation, fluid dynamics, etc. will occur and at what rates.
The main program also imports modules (i.e., groups of functions) for diversity metrics, spatial analysis, the initiation of output files, and for simulating ecological life history processes (immigration, maintenance, death, growth, consumption, disturbance, passive dispersal, active dispersal, resource flow, resource inflow, and metabolic state transitions).

**Simulating life history**
-- Simplex models begin simulation immediately after assembly. 
The order of life history processes is randomized at each time step.
At each time step, each individual is given the chance to consume, grow, reproduce, starve, die, and to disperse.
Simplex models simulate growth, reproduction, and death via weighted random sampling. 
This approach is intended to simulate the partly probabilistic and partly deterministic nature of ecological processes.

*Immigration:* 
Individuals can be allowed to enter at any point in the environment, which can be adjusted in the 'immigration.py' file.
Species identities of inflowing propagules are chosen at random from a uniform distribution rather than an ecologically realistic source pool (e.g., log-series or lognormal distribution).
The reason for this is two-fold.
First, pulling immigrants from a uniform distribution maximizes the starting diversity of each model.
Second, beginning with a uniform distribution allows more realistic distributions of species abundances to emerge.

*Active dispersal:* 
Individuals can actively move against the force of flow, at random, or in specified directions.
Preferences for particular modes of movement can be specifed in the 'active_dispersal.py' file.

*Passive dispersal:*
Individuals can optionally be moved passively (e.g., as planktonic organisms) through the system at rates determined solely by the rate of flow through or by the combination of flow rate and active dispersal.

*Consumption:* 
Sampled individuals increase their levels of endogenous resources by feeding on resource particles.
These endogenous resources (a.k.a. cell quotas) can be used to add structural biomass and to pay the energetic costs of life history processes. 
Individual consume resources according to their species specific consumption rates for three resource types.
The number of simulated resource types can be changes in the source code files (consume.py, immigration.py, and resource_inflow.py). 
The simulation of consumption can be modified in the 'consume.py' file.

*Growth:* 
Sampled individuals grow in size by integrating endogenous resources as structural biomass.
Individuals grow according to species specific rates of growth ranging between 0.1% and 100% increase in size per time step.
Individuals' endogenous resources are decreased in direct proportion.

*Reproduction:* 
Individuals reproduce with a probability determined by the combination of their endogenous resources and growth rate. 
Sampled individuals reproduce via clonal reproduction with the possibility of 
mutation.
The endogenous resources of the mother individual is evenly divided between two 
daughter individuals. 
Unless in the case of speciation, the daughters are given a unique individual ID and the species ID of the mother.

*Speciation:* 
Speciation is simulated within Simplex occurs as a discrete event, i.e., where clonal reproduction produces a new species.
Speciation is accompanied by mutations in the values of one or more species traits. 
This approach allows for diversity to arise within the system, which the environment can then select on.

*Death:* Individuals sampled at random will die if their levels of endogenous resources or ability to draw resources from structural biomass falls below the minimum metabolic requirements.
Dead individuals are removed from the system, i.e., scavenging and recycling do not currently occur in Simplex.

*Emigration:* Individuals are considered to have left the environment when they pass beyond edges of the environment.

**Simulating resource processes**

*Supply/inflow:*
Resource particles can be allowed to enter at any point in the environment, which can be adjusted in the 'resource_inflow.py' file.
The size and type of each inflowing resource particle is chosen at random from a uniform distribution, and can be modified in the same 'resource_inflow.py' file.

*Resource dispersal:*
Resource particles can be moved passively through the system at rates determined by the rate of flow through.
Resource dispersal can be turned off in the 'main.py' file.

*Resource depletion:*
Resource particles are depleted through consumption, which can be partial or complete.

### Design concepts

**Basic principles** 
*Ecological selection on random variation:* 
Simplex operates according to a basic principle of evolution, i.e., natural selection on random variation.
More specifically, Simplex simulates ecological selection or environmental selection, which refers to natural selection without respect to sexual selection.
Simplex assembles models from random combinations of system-level, species-level, and individual-level variables.
Simplex then allows the environmental characteristics (e.g., flow rate, resource supply, spatial extent) to select on these random trait combinations.

*Energy-limited life history:*
All Simplex models impose energetic costs to growth and activity.
These energetic costs are directly proportional to life history parameters.
For example, the energetic cost of dispersal is the product of dispersal rate and dispersal distance.
This intuitively means that the energetic cost of dispersal is multiplied (or compounded) across distance such that moving a distance of x requires half the energy as moving a distance of 2x.
In the same way, growing at a rate of x is half as costly as growing at a rate of 2x.
Simplex does not explicitly force allometric relationships because one goal of Simplex is to allow allometric relationships to emerge from variation in measureable variables (e.g., growth rate, dispersal rate).

*Lognormal dynamics:*
Multiplicative interactions of random variables underpin one of the most successful models of complex systems, i.e., the lognormal (Crow et al. 1988).
The lognormal was introduced to ecologists by Preston (1962) as a statistical description of how abundance varies among species.
The lognormal has been one of the two most successful models of the species abundance distribution for macroscopic plants and animals, and was recently used to form a macroecological theory of microbial diversity.
By "multiplicative interactions" one simply means that two or more variables or processes have multiplicative (i.e., synergistic) interactions.
Such interactions are common in ecology (Putnam 1993) and complex systems and result in non-linear changes (e.g., multiplicative population growth, energetic costs multiplied across distance).
By "random variables" one simply means two independent processes or constraints with degrees of stochastic change.
Simplex explicitly simulates lognormal dynamics. For example, energetic costs are multiplied across dispersal distance and magnitudes of biological change (e.g., growth, consumption, transition between metabolic states). 
These energetic costs are determined by the values of randomly chosen species traits (i.e., random variables).
Finally, the occurrence of each life history process is simulated in random order to avoid systemic programmatic biases and order-dependence among processes.

*Simultaneous emergence:*
One of the most popular advantages of using IBMs is the potential to study the emergence of complex patterns and relationships from otherwise simple individual-level processes and interactions.
Simplex is aimed at allowing emergence in two ways that IBMs are rarely employed.

First, Simplex allows realism to emerge from several orders of magnitude in random variation of starting conditions. 
Simplex is intended to run thousands of models, where each is initiated with random conditions and combinations of species traits.
In this way, Simplex models initiate with highly unrealistic ecological communities with highly unrealistic combinations of species traits, and then allows for realism to develop over time steps as a result of appropriately modeled processes, energetic costs, and lognormal dynamics.
This approach allows the user of Simplex to avoid one of the greatest challenges to ecological modeling, i.e., the circularity of documenting patterns, relationships, and trade-offs that are practically forced to occur.
The code made available for the analysis of Simplex data is intended to examine the variation and central limiting behavior of thousands of models, some of which can defy ecological intuition.

Second, the Simplex user can examine how multiple ecological relationships and patterns emerge together. 
A hallmark of powerful ecological theory is the ability to predict multiple relationships and patterns, yet few ecological theories and models predict more than one or two patterns.
In contrast, Simplex includes python code to analyze the simultaneous emergence of patterns of metabolic scaling, species abundance distributions, species-area relationships, diversity-abundance scaling relationships, Taylor's scaling law, the influence of ecosystem residence time on diversity and abundance, diversity-productivity relationships, among others.

**Other design concepts**
*Hypotheses:* 
These are entirely up to the user to formulate and test according to the capabilities and analytical tools of Simplex source code.

*Learning:* 
Currently, there is no aspect of individual-based learning in Simplex.

*Prediction:*
Individuals in simplex do not have the ability to anticipate conditions.

*Sensing:*
Individuals can sense and move towards resource particles.

*Interaction:* 
Individuals interact through excluding each other from resources, i.e., 
preemption. 
There is no explicit communication.

*Observation:*
An unlimited number of Simplex models can be run to examine trends and variation in ecological patterns.
The following is recorded for each Simplex model:

* Values of randomly chosen input variables
	* length
	* width
	* flow rate
	* maximum resource particle size
	* maximum number of inflowing resource particles
	* number of inflowing resource types
	* initial community size 
* Total individual abundance ($N$) of the dormant and active fractions of ecological communities
* Species richness ($S$) of the dormant and active fractions of ecological communities
* Mean abundance-weighted specific growth rate for the dormant and active fractions of ecological communities
* Mean abundance-weighted species-specific dispersal rate for the dormant and active fractions of ecological communities
* Mean abundance-weighted species-specific resource consumption rates for the dormant and active fractions of ecological communities
* Mean abundance-weighted species-specific rate of random metabolic state transitions (dormant to active) for the dormant and active fractions of ecological communities
* Mean abundance-weighted species-specific rate of decrease in basal metabolic rate accomplished by transitioning to a dormant state, for the dormant and active fractions of ecological communities
* Species turnover
	* Whittaker's $\beta$ 
* Species evenness
	* Smith and Wilson's evenness, $E_{var}$
	* Simpson's evenness, $E_{1/D}$ 
* Species diversity
	* Shannon's diversity, $H'$ 
	* Simpson's diversity, $D_{1/D}$
* Dominance
	* Absolute, $N_{max}$
	* Relative, $N_{max}/N$
* Productivity of
	* Individuals
	* Biomass

These data are stored as comma separated value (.csv) files, the data from which can be directly imported into an R or Python environment.

## Output data
Simplex generates three files of output data. 
Each Simplex model quantifies and writes output data for every $n^{th}$ time step, where $n$ can be designated by the user. 
The three .csv files are: 

*SimData.csv*
-- Each column of this file corresponds to a piece of data about the system that was modeled (e.g., flow rate, total abundance, species richness, species turnover, resource supply and diversity, rate of disturbance, etc.).
Most analyses in Simplex source code are conducted on the data in this file.

*SAR.csv*
-- A file holding species-area relationships (SARs) from Simplex models.
SARs quantify the rate at which species are encountered with increasing area of a sample, study, landscape, etc.
The SAR is among the most intensively and long-studied patterns in ecology and is one of two patterns commonly predicted by biodiversity theories (Lomolino 2000, Hubbell 2001, Harte 2011).

*RADs.csv*
-- A file holding rank-abundance distributions (RADs) from Simplex models.
Also referred to as species-abundance distributions (SADs), rank-abundance curves (RACs), and Whittaker plots, RADs are vectors of the abundances of species in a community.
Along with the SAR, RADs are one of the most intensively studied and commonly predicted ecological patterns (Hubbell 2001, McGill et al. 2007, Harte 2011).

# Discussion

## UNDER CONSTRUCTION

The Simplex platform combines ecological paradigms and general mechanisms to allow many ecological patterns and relationships to emerge simultaneously. 
In this way, the relatively simple source code of Simplex allows patterns of community ecology, macroecology, and biodiversity science to emerge from individual-level changes and relatively few explicit constraints. 
Simplex assembles individual-based models from random combinations of system-level variables, species traits, and process-related rates. 
This large degree of initially random conditions allows users to explore a large degree of variation in ecological conditions and to avoid the circularity of building overly-informed models. 
In this way, Simplex allows the user to examine the simultaneous emergence of ecological relationships under wide-ranging conditions and combinations of ecological paradigms.

Simplex is designed to explore and test the simultaneous emergence of ecological relationships and other patterns. 
The usefulness of this is manifold. 
First, ecologists often have little knowledge of how ecological patterns are related or either interdependent or even perhaps mutually exclusive. 
Simplex is perhaps ecology's first tool for targeting such questions. The ability to examine simultaneous emergence can also be used to great effect when testing or developing ecological theory. While many  ecological theories predict only one or a few patterns, strong theories should unify and predict multiple patterns (McGill et al. 2007). 
Likewise, strong tests should evaluate the ability of ecological theories to explain and predict multiple patterns (e.g., Xiao et al. 2015).

Simplex has also been designed to combine concepts from several ecological paradigms into a single platform that offers source code that is easy to read, use, and modify. 
By default, Simplex includes the stochastic dynamics inherent to some ecological theories (e.g., neutral theory, stochastic geometry, stochastic resource limitation theory). 
Likewise, Simplex includes the resource-limited growth dynamics of resource limitation theory and chemostat theory, and the inherent species sorting and environmental filtering (i.e., ecological selection) of community ecology. 
Finally, Simplex includes the lognormal dynamics of a recent macroecological theory, i.e., multiplicative interactions of stochastic variables, in addition to energetic constraints that underpin trade-offs of life history theory (e.g., r vs. K selection).

Simplex offers the first individual-based modeling platform for metabolic scaling theory. 
Metabolic scaling predicts that the magnitude of basal metabolic rate ($B$) increases to a fractional power of body mass ($M$), i.e., $B$ = $M^{z}$. 
This relationship is most often observed to follow a 3/4 power scaling, i.e., $z$ = 3/4. 
The 3/4 scaling law has become one of the most universal and highly supported statistical relationships in the biological sciences. 
MTE and many related studies have used this scaling law to predict aspects of metabolic power, population dynamics, community ecology, ecosystem function, and trophic interactions (*refs*). 
The reasoning behind the 3/4 power scaling law rests on fractal branching networks responsible for the delivery of resources to cells and tissues. 
However, scientists still argue whether a 2/3 scaling based on surface area to volume ratios is more accurate and appropriate (*refs*) and even whether a 3/4 and 2/3 scaling can arise from the same mechanism (*refs*).
Finally, some have shown that the scaling of metabolic rate to body size is nearly isometric (0.9 < $z$ < 1) for some taxa (*ref*).




5) LOOKING TO THE FUTURE -- Future developments will include additional ecological dynamcis such as predator-prey, mutualism, and parasitism. 
Improvements to Simplex will also include an accounting of nutrient stoichiometry, and Phosphorus, as well as degrees of biocomplexity. For example, given the chemical structure of phosphate and the size of a resource particle, the particle could be assigned $x$ units of phosphorus as well as the biocomplexity value of phosphate (estimated as a form of Shannon's entropy). 
Likewise, improvements and future versions of Simplex will provide increasing numbers of files for statistical analysis of simulated outputs.

