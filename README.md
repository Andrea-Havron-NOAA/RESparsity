# RESparsity

Stock assessment models have historically used a non-centered parameterzation approach when estimating variability around recruitment. In the Bayesian context, the non-centered parameterization is often preferred due to better mixing and convergence in MCMC sampling. Little research has looked into this issue within the MLE context. With respect to time series analysis, it is hypothesized that the centered approach results in sparse covariance matrices resulting in computational efficiency gains. Additionally, the Moving Average (MA) formulation has been little explored in stock assessments and a comparison with AR1 will be informative for improving random effect formulations in the future. Part of the reason is it's ineffiecncy. Investigating how to improve the sparisty of MA1 could lead to increased adoption of this approach. NOAA’s next generation stock assessment modeling platform, the Fisheries Integrated Modeling System (FIMS) aims to incorporate best available science while meeting the historical needs of fishery assessments, therefore changes to formations are best justified through scientific analysis and documentation. This manuscript aims to answer the following questions:

1. Is there any difference between centered / non-centered parameterizations in maximum likelihood (TMB) world?
   * Answer: no, even with short time series, and high observation sd : process sd, all formulations give the same answer
2. If there is a difference, can we provide recommendations around TMB/FIMS using simulated data and / or case studies? 
3. Can we formulate the MA1 so that it is efficient with the laplace approximation?

Links: (co-authors only)
[Project Proposal](https://docs.google.com/document/d/1yXtSZJWwB_VZeUHBTA4y-8gTtps1mipRjXdAi7fhZ8U/edit?tab=t.0#heading=h.m5crj94thbbs)
[Meeting Notes](https://docs.google.com/document/d/1pZx1sor-ab2wq5UTF5hCFT0Phi7ktnagjXMP50nIErM/edit?tab=t.0)
