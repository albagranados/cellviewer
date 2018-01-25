# Spatial data analysis of single molecule localization microscopy

# R code

Point pattern data describes data in which random events are observed over some domain, with the number and locations of these events (points) being random. 

Statistical spatial point pattern analysis is concerned with describing short-range correlations and spacings (or shortest distances) among points for stochastic dependency assesment.

In a marked point processes each point is assigned a mark. This is the case in dual color dSTORM, where each molecular event is assigned  a qualitative mark (or label) corresponding to the molecule type.

Methods:

- Summary characteristics: assessing dependency by correlations: Multivariate or inter-type K-function; Cross (or partial) pair correlation function

Remember: Summary functions do not completely characterise the point process. Correlation is not causation! If CSR -> no correlation; <- is not true!
