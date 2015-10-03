---
output: pdf_document
---
# Introduction

Modeling is a elementary approach to understanding ecological systems, the influence of ecological processes, and the predictability of ecological patterns and dynamics. 
Before modern computing, ecological models were almost exclusively equation-based representations of highly simplified systems (Black and McKane 2012, Otto and Troy 2007, Grant and Swannack 2008). 
Since the advent of personal computing, ecological models have been increasingly constructed to handle greater complexity and to explicitly simulate ecological processes.

In ecology, simulation-based models are often used to examine analytically challenging or intractable scenarios. 
Examples are markov models that simulate stochastic demographic changes and early ecological null models, both of which operate on matrices (Gotelli and Entsminger 2001, Hubbell 2001). 
Ecologists have also simulated growth and interactions among individuals for over two decades with individual-based models (IBMs) (DeAngelis and Gross 1992, Rosindell et al. 2015). 
IBMs explicitly encode rules of how individuals change and interact. 
Once the rules are encoded, population to ecosystem-level dynamics can emerge as an IBM simulates over time and spatially explicit environments.

A body of literature including original research, comprehensive reviews, and textbooks reveal the use, advantages, and challenges of ecological IBMs (e.g., Grimm 1999, Grimm and Railsback 2005). 
IBMs can provide degrees of ecological realism, individual variability, and spatial heterogenetiy that are unattainable with other models. 
IBMs also offer the potential for realistic and unanticipated ecological dynamics and patterns to emerge from individual-level interactions. 
However, IBMs offer challenges that include greater computational complexity, the difficulty of explicitly encoding ecological theory, and the use of ecological problem solving (Matthews et al. 2007, Grimm and Railsback 2005).

Despite their frequent use, the power of ecological IBMs has yet been leveraged to the greatest advantage. 
Few, if any, include fluid dynamics or combine fluid dynamics with growth and active dispersal. 
While growth, sensing, and decision making are often modeled, few IBMs integrate physiology, evolution, community ecology, biogeography, and sampling theory. 
Ecological IBMs are more often constructed to model a specific system than to synthesize general theories of ecology and evolution (but see Rosindell et al. 2015). 
Yet, IBMs can allow researchers to simulate and track information from the level of genomes and internal physiology to the stoichiometry of resource particles and distributions of abundance among species and across space.

Here, I present source code and a detailed description for an ecological IBM platform (called simplex) that simulates a broad range of ecological conditions, records information from individuals to the ecosystem, and provides computing code for analyzing output. 
This platform incorporates computational fluid dynamics (Succi 2001), life history processes and concepts and patterns from biogeography (Hubbell 2001, Bell 2001), evolution and aspects of selection (Hartl and Clark 1997), and integrates nutrient-limited growth and physiology (Pirt 1965, Droop 1983) among dimension of ecology and evolutionary science. 
I constructed simplex with the aim of simulating ecologically complex scenarios, from which the results could be used to test (*i*) predictions of ecological theory and (*ii*) the general influence of constraints and processes. 
Below, I provide detailed explanation of how hydrobide works, the data it quantifies and tracks, the theories and principles hydrobide integrates, and the analyses that can be conducted using the code we provide.