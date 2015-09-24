# Methods
**Note:** Some numbers in parentheses correspond to references; a reference section is being built.

## Data
###P1
**S1:** We used 15,535 sites of communities of bacteria, archaea, and microscopic fungi.  
**S2:** 14,962 of these sites were from the Earth Microbiome Project (EMP) (14) obtained on 22 August, 2014.  
**S3:** Sample processing, sequencing and amplicon data are standardized and performed by the EMP and all are publicly available at www.microbio.me/emp.  
**S4:** The EMP data consist of open and closed reference datasets, which are defined in the QIIME tutorial (http://qiime.org/tutorials/otu_picking.html).  
**S5:** QIIME defines closed-reference as a classification scheme where any reads that do not hit a sequence in a reference collection are excluded from analysis.  
**S6:** In contrast, open-reference refers to a scheme where reads that do not hit a reference collection are subsequently clustered de novo and represent unique but unclassified taxonomic units.  
**S7:** Our main results are based on closed-reference data, due to the greater accuracy of the approach and because unclassified sequences were excluded from other microbial datasets (below).	
###P2
**S1:** We also used 4,303 sites from the Data Analysis and Coordination Center (DACC) for the National Institutes of Health (NIH) Common Fund supported Human Microbiome Project (HMP).  
**S2:** These data consisted of samples taken from 15 or 18 locations (including the skin, gut, vagina, and oral cavity) on each of 300 healthy individuals.  
**S3:** The v3-v5 region of the 16S rRNA gene was sequenced for each sample.  
**S4:** We excluded sites from pilot phases of the HMP as well as time-series data; see http://hmpdacc.org/micro_analysis/microbiome_analyses.php. for details on HMP sequencing and sampling protocols.

###P3
**S1:** We also included 1,319 non-experimental PCR-targeted rRNA amplicon sequencing projects from the Argonne National Laboratory metagenomics server MG-RAST (16).  
**S2:** Represented in this compilation were samples from arctic aquatic systems (130 sites; MG-RAST id: mgp138), hydrothermal vents (123 sites; MG-RAST id: mgp327) (37), freshwater lakes in China (187 sites; MG-RAST id: mgp2758) (38), arctic soils (44 sites; MG-RAST id: mgp69) (39), temperate soils (84 sites; MG-RAST id: mgp68) (40), bovine fecal samples (16 sites; MG-RAST id: mgp14132), human gut microbiome samples not part of the HMP project (529 sites; MG-RAST id: mgp401) (41), a global-scale dataset of indoor fungal systems (128 sites) (42), and freshwater, marine, and intertidal river sediments (34 sites; MG-RAST id: mgp1829). 

###P4**S1:** The use of MG-RAST allowed us to choose common parameter values for percent sequence similarity (i.e. 97% for species-level) and taxa assignment including a maximum e-value (probability of observing an equal or better match in a database of a given size) of 10-5, a minimum alignment length of 50 base pairs, and minimum percent sequence similarities of 95, 97, and 99% to the closest reference sequence in MG-RAST’s M5 rRNA database (37-42). **S2:** Quantifying dominance, evenness, rarity, and richness. 
**S3:** We calculated or estimated aspects of diversity (dominance, evenness, rarity, richness) for each site in our data compilation. 
**S4:** All analyses can be reproduced or modified for further exploration by using code, data, and following directions provided here: https://github.com/LennonLab/MicroMETE.  

## MaxEnt predictions of the SAD
### METE
**P1**  
**S1:** The maximum entropy theory of ecology (METE) (Harte et al. 2008, 2009, Harte 2011) is based on two empirical inputs: species richness (*S*) and total abundance (*N*).   
**S2:** These, along with an inferred rate of community-level metabolism (*E*), form the state variables of METE.  
**S3:** Four constraints are produced from these state variables.  
**S4:** These are the average number of individuals per species (*N*/*S*), the average per species metabolic flux (*E*/*S*), and the constraints that no species has more than *N* individuals or a greater total metabolic rate than *E*.  
**S5:** *E* is later integrated out of the SAD prediction.  

**P2**  
**S1:** The prediction of METE is based on a joint conditional probability distribution that describes the distribution of individuals (*n*) over species and of metabolism (*ε*) over individuals within a species (Harte et al. 2008, Harte 2011). 
**S2:** Entropy of the distribution is then maximized according to the method of Lagrangian multipliers (Jaynes 2003, Harte 2011).  
**S3:** The SAD is then derived by integrating out energy and dropping terms that are vanishingly small. This process then yields the log-series SAD (Fisher et al. 1943).  
**S4:** The log-series distribution is among the oldest and most successful SAD models but has generally lacked a convincing first-principle explanation from either an ecological or statistical perspective.  
**S4:** In this case, METE predicts the shape of which is dependent only on the values of *S* and *N*:

###insert equation 1

where β is defined by the equation 

###insert equation 2


### Broken-stick 
**S1:** While some other MaxEnt models produce similar, if not, identical (Pueyo et al. 2007, Dewar and Porté 2008, Frank 2011) predictions for the SAD, MaxEnt models based on different assumptions can yield very different predictions (Haegeman and Etienne 2010).   
**S2:** One example is the simultaneous discrete Broken-stick model of MacArthur (1960), which as pointed out by Haegeman and Etienne (2010) is simply the geometric distribution with mean *N*/*S*.   
**S3:** Unlike the log-series, the broken-stick model predicts a relatively even distribution which is often a poor fit to empirical SADs (Hubbell 2001).  
**S4:** The broken-stick gives equal weight to all ordered configurations of *S* species whose abundances sum to *N*, the equation for which is:

###insert equation 3

## Testing MaxEnt predictions
**S1:** Both METE (which predicts a log-series distribution) and the Broken-stick (i.e., the geometric distribution) produce predictions for the rank-abundance form of the SAD.   
**S2:** This form of the SAD is simply a vector of species abundances ranked from greatest to least.  
**S3:** Both predictions yield the same value of *S* that is given as the empirical input.  
**S4:** This means that the observed and predicted SADs can be directly compared using regression analyses to reveals the percent variation explained by each model (METE, Broken-stick).  
**S5:** We generated the predicted forms of the SAD using the source code of White et al. (2012) (https://github.com/weecology/white-etal-2012-ecology) and the macroecotools repository (https://github.com/weecology/macroecotools), which contains functions for fitting maximum-likelihood forms of species abundance models in addition to functions for other macroecological analyses.  
**S6:** Using that source code, we calculated the modified coefficient of determination (r-square) around the 1-to-1 line (as per White et al. 2012, Locey and White 2013, Xiao et al. 2014).