#HYDRO-BIDE

Python code authored by Ken Locey for a project using individual based modeling and an array of spatial dynamics to examine the influence of residence time on biodiversity. This code accompanies a manuscript and generates the figures to be published. Yay, for reproducible science!

##Description
**Hydrobide** a probabilistic individual-based model for exploring how aspects of abundance and diversity respond to physical constraints such as residence time, biological constraints like maximum growth rate, and ecological complexity arising from fluid and non-fluid dynamics, resource diversity, open/closed system dynamics, and predation. Let's parse some of that out:

* **Probabilistic:** changes to populations and communities across time are brought about through random sampling, either unweighted (i.e. neutral) or weighted to reflect the probabilistic nature of individual responses to the environmental filter, competition dynamics, predator-prey dynamics, dispersal, etc.

* **Individual-based:** Every individual organism, resource particle, predator, and inert tracer is tracked, to include the amount of time each has spent in the system. For example, each individual organism is simulated and has an individual ID, a species ID, an evolutionary heritiage, an amount of internal resources (cell quota), species-specific dispersal capacities and resource efficiency values for all resources in the system, and x, y, z spatial coordinates.

* **Aspects of abundance and diversity:** 
	* Total abundance: the number of individuals in a population, community, or ecosystem)
	* Productivity: number of individuals and/or biomass produced. Biomass is calculated from the model by summing individual cell quotas; don't go crazy, it's just a proxy because it's just a model and the organisms in it don't weigh anything ;)
	* Taxonomic richness: the number of taxa in a community or ecosystem
	* Turnover: Also referred to as beta-diversity; a measure of heterogeneity in the taxa found among samples. Due to the nature of the modeling, individual-level turnover is also possible.
	* Evenness: Similarity in abundance among taxa.

## Using simulated data in R
The output of hydrobide can be imported into Python and R environments as dataframes. R Markdown scripts are provided.
