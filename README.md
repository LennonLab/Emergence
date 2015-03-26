#HYDRO-BIDE

Python code authored by Ken Locey for a project using individual based models to explore the influence of residence time on aspects of biodiversity. This code accompanies a manuscript and generates the figures to be published. Yay, for reproducible science!

##Description
**Hydrobide** is a stochastic individual-based model for exploring how aspects of abundance and diversity respond to physical constraints such as residence time, biological constraints like maximum growth rate, and ecological complexity. Let's parse some of that out:

* **Stochastic:** changes to populations and communities across time are brought about through random sampling, either unweighted (i.e. neutral) or weighted to reflect the probabilistic nature of individual responses to the environmental filter, competition, and non-competitive dynamics.

* **Individual-based:** each individual is simulated and has, at least, an individual ID, a species ID, an amount of internal resources (cell quota), and x and y spatial coordinates.

* **Aspects of abundance and diversity:** 
	* Total abundance: the number of individuals in a population, community, or ecosystem), 
	* Taxonomic richness: the number of taxa in a community or ecosystem
	* Turnover: Also referred to as beta-diversity; a measure of heterogeneity in the taxa found among samples
	* Evenness: Similarity in abundance among taxa


