# Seabird_trends

This repo includes everything necessary to re-create the seabird population trend analyses that were conducted for the 2019 State of Canada's Birds. 
It's a mix of hierarchical GAM models to track populations among colonies, and/or among long-term monitoring plots nested within colonies. 
The GAMs were parameterised using code from Crainiceanu et al. 2005, or in the case of the WEst coast plots, using base models built on the defaul basis functions from the jagam() function in the R-package mgcv.

All of these models could be easily transformed to run in brms. There's nothing particularly complicated about the models.

If one of the required datasets is missing, please let me know ASAP. I have everything stored locally, but past-me sucked at file-management even more that present-me does, so I have far more files saved than I could includ here.



