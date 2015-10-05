#Results

## Unit tests

simplex passed all units tests for 15 diversity indices, ensuring that each index returns either the correct calculated value, or 'NaN' if given any values that cannot be used (e.g., negative numbers, string characters, empty lists).

## Speed & Memory

**Results 100 randomly assembled models.** On average, simplex models required 96.60 +- 32.70 seconds to run, had an average total abundance (*N*) of 31367.0 +- 878.92, and required 160.5 +- 4.28 megabytes of memory.
The longest any model took to complete was 48.95 minutes. This model required 261MB of memory and also had generated the greatest total abundance (*N* = 68,808). The shortest any model took to complete was 5.8 seconds with a total abundance of 0 individuals. The least amount of memory any model required was 90MB. These and other analyses from the 100 randomly assembled models can be run using the SpeedMemory.Rmd file located in the results/TestResults directory.

## Output data
simplex generates six files as its output data. They are: 

*SimData.csv*--Formatted as an R data frame, where each row is a run of a randomly assembled model, and each column holds a piece of data about the system that was modeled (e.g., flow rate, total abundance, species richness, species turnover, rate of disturbance, etc.).  

*IndRTD.csv*, *ResRTD.csv*, *TracerRTD.csv*--Each line of these three files contains the amount of time that each individual, resource particle, or inert tracer (from beginning to end) spent in the system. Each line is a run of a randomly assembled model.

*RADs.csv*--Each line is a run of a randomly assembled model and contains the rank-abundance vector at the point when the model was stopped.

*Species.csv*--Each line is a run of a randomly assembled model and contains a vector of species labels corresponding to RADs.csv.

All output data were correctly formatted and placed within results/simulated_data/examples directory.
Each of the six output data files was able to be imported into the RStudio environment using the Exploratory.Rmd Rmarkdown file provided in the "GitHub/simplex/results/analyses/Rmd/" path. The user can use the Exploratory.Rmd file to craft a .Rmd file for their own specific analyses.

## Analysis of saved output

Running the Exploratory.Rmd file conducts the following analyses:






## Animation





