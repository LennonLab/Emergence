#Results

## Unit tests

simplex passed all units tests for 15 diversity indices, ensuring that each index returns either the correct calculated value, or 'NaN' if given any values that cannot be used (e.g., negative numbers, string characters, empty lists).

## Speed & Memory

Simplex models do not complete until a specific number of tracer particles have left the system. Because simplex models can range from flowing very quickly to hardly flowing at all (e.g., moving between 10% and 0.0001% of the landscape each time step) simulations can potentially take several minutes or more to complete. 

Running simplex on Mid 2010 MacBook Pro (OS X 10.9.5) with a 2.4 GHz Intel Core 2 Duo processor and 4GB of Memory, we found that completion of 100 models run at the slowest rate of environmental flow (i.eg. 0.0001% of the landscape at each time step) took, on average ...minutes. and carried in memory 100MB of information.

## Output data

**SimData.csv**--

**IndRTD.csv**--

**RADs.csv**--

**ResRTD.csv**--

**Species.csv**--

**TracerRTD.csv**--






## Animation





## Analysis of saved output