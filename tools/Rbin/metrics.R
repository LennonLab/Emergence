################################################################################
#                                                                              #
# Functions for calculating metrics of diversity, evenness, rarity, etc.       #
# These are not included in other diversity packages, e.g., Vegan              #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Ken Locey                                                        #
#                                                                              #
################################################################################
#                                                                              #
# Recent Changes:                                                              #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1.                                                                   #
#         2. Add warnings                                                      #
#                                                                              #
################################################################################

#### A function to generate observed richness
S.obs <- function(x = ""){ rowSums(x > 0) * 1}

sp.turnover <- function(site1, site2){
  
  site1 <- site1[site1 != '']
  site1 <- site1[!is.na(site1)]
  site2 <- site2[site2 != '']
  site2 <- site2[!is.na(site2)]
  #if(length(site1) | length(site2) == 0){
  #  return -1
  #  }
  gamma = union(site1, site2)         # Gamma species pool
  s     = length(gamma)                                   # Gamma richness
  a.bar = mean(c(length(site1), length(site2)))   # Mean sample richness
  b.w   = round(s/a.bar - 1, 4)
  return(b.w)
  }

