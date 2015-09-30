# Introduction

Modeling is a primary approach to understanding ecological systems. 
Before modern computing, ecological models were almost exclusively equation-based representations of highly simplified systems. 
Examples, are basic Lotka-volterra predator-prey models, simple models of community structure, and mark-recapture models. 
Since the advent of personal computing, ecological models have been increasingly constructed to handle greater complexity and to explicitly simulate ecological processes.

In ecology, simulation-based models are often used to examine analytically challenging or intractable scenarios. 
Examples are markov models that simulate stochastic demographic changes and early ecological null models, both of which primarily operate on matrices (Gotelli and Entsminger 2001, Hubbell 2001). 
At a more basic level, ecologists have been simulating growth and interactions among individuals for over two decades with individual-based models (DeAngelis and Gross 1992). 
In short, individual-based models (IBMs) explicitly encode rules of how individuals change and interact. 
Once the rules are encoded, population to ecosystem-level dynamics can emerge as an IBM simulates over time and spatially explicit environments.

A body of literature including original research, comprehensive reviews, and textbooks reveal the use, advantages, and challenges of IBMs in ecology (e.g., Grimm 1999, Grimm and Railsback 2005). 
IBMs can provide degrees of ecological realism, individual-variability, and spatial heterogenetiy that are unattainable with analytical frameworks and other simulation models. 
IBMs also offer the potential for realistic and unanticipated ecological dynamics and patterns to emerge from individual-level interactions. 
However, IBMs offer challenges that include greater computational complexity and the difficulty of explicitly encoding ecological theory.  

Despite their frequent use, the power of IBMs has yet been leveraged to the greatest advantage. 
Few, if any, ecological IBMs include fluid dynamics or combine fluid dynamics with growth and active dispersal.  
Likewise, while growth and response to stimuli are often modeled in detail, few IBMs integrate physiology, evolution, community ecology, biogeography, and sampling theory. 
More often than not, ecological IBMs are constructed to model a highly specified system as the result of deterministic rules.  
IBMs are less used as in silico environments for synthesizing ecological and evolutionary theory. 
Yet, IBMs can allow researchers to simulate and track information from the level of genomes and internal physiology to the stoichiometry of resource particles and the aggregate states of an ecosystem.

I constructed source code for a general IBM platform that simulates a broad range of ecological conditions, records information from individuals to the ecosystem, and provides computing code for analyzing output.  
**S2:** We refer to this platform as hydrobide to denote two of its primary features: the inclusion of fluid dynamics and the central role of birth (b), immigration (i), death (d), and emigration (e) processes.  
**S3:** We constructed hydrobide with the aim of simulating ecologically complex scenarios, from which the results could be used to test (*i*) predictions of ecological theory and (*ii*) the general influence of constraints and processes.  
**S4:** Here, we provide detailed explanation of how hydrobide works, the data it quantifies and tracks, the theories and principles hydrobide integrates, and the analyses that can be conducted using the code we provide.