# Introduction

##P1: Ecological modeling
**S1:** Modeling has been a primary approach to understanding ecological systems for nearly a century.  
**S2:** Before modern computing, ecological models were almost exclusively equation-based representations of highly simplified systems.  
**S3:** Examples, are basic Lotka-volterra predator-prey models, simple models of community structure, and mark-recapture models.  
**S4:** However, since the advent of personal computing, ecological models have been increasingly constructed to handle greater analytical complexity and to explicitly simulate ecological processes.

##P2: Individual-based modeling
**S1:** In ecology, simulation-based models are often used to examine the outcomes of analytically intractable scenarios.  
**S2:** Examples are markov models that simulate demographic changes and early ecological null models, both of which primarily operate on matrices (refs).  
**S3:** At a more basic level, ecologists have been simulating growth and interactions among individuals for over two decades with individual-based models (DeAngelis and Gross 1992).  
**S4:** In short, individual-based models (IBMs) require explicitly encoded rules of how individuals change and interact.  
**S5:** Once the rules are encoded, population to ecosystem-level dynamics can emerge as an IBM simulates over time and spatially explicit landscapes.

##P3: Advantages and challenges of IBMs  
**S1:** A body of literature including research, reviews, and textbooks reveal the use, advantages, and challenges of IBMs in ecology (Grimm 1999, Grimm and Railsback 2005).  
**S2:** IBMs can provide degrees of ecological realism, individual-variability, and spatial heterogenetiy that are unattainable with analytical frameworks and other simulation models.  
**S3:** Another primary advantage is the potential for realistic and unanticipated ecological dynamics and patterns to emerge from individual-level interactions.  
**S4:** However, IBMs also offer challenges including greater computational complexity and the difficulty of encoding ecological theories that are not explicitly individual-based.  

##P4: Advancing ecological IBMs
**S1:** Despite their increased use, the power of IBMs has yet been leveraged to the greatest advantage.  
**S2:** Few, if any, IBMs include fluid dynamics or combine flowing environments with active dispersal.  
**S3:** Likewise, while growth and response to stimuli are often modeled in detail, few, if any, IBMs integrate physiology, evolution, community ecology, biogeography, and sampling theory.  
**S4:** More often than not, ecological IBMs are constructed to model a highly specified system as the result of deterministic rules.  
**S5:** IBMs are less used as in silico environments for testing and synthesizing ecological and evolutionary theory.  
**S6:** Yet, IBMs can allow researchers to simulate and track information from the level of genomes and internal physiology to the stoichiometry of resource particles and the aggregate states of an ecosystem.

##P5: Hydrobide
**S1:** We constructed source code for a general IBM platform that simulates a broad range of ecological conditions, records information from individuals to the ecosystem, and provides computing code for analyzing output.  
**S2:** We refer to this platform as hydrobide to denote two of its primary features: the inclusion of fluid dynamics and the central role of birth (b), immigration (i), death (d), and emigration (e) processes.  
**S3:** We constructed hydrobide with the aim of simulating ecologically complex scenarios, from which the results could be used to test (*i*) predictions of ecological theory and (*ii*) the general influence of constraints and processes.  
**S4:** Here, we provide detailed explanation of how hydrobide works, the data it quantifies and tracks, the theories and principles hydrobide integrates, and the analyses that can be conducted using the code we provide.