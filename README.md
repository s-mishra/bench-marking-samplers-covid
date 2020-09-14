# bench-marking-samplers-covid
Code from our covid19 model to be used for benchmarking various samplers.

To run the code from terminal do following
```
cd gitrepo
# run europe script
Rscript src/base-europe.R -F
# run usa script
Rscript src/base-usa.R -F
```
To run the code from rstudio
```
# source europe script, make sure you make full true in options
source base-europe.R
# source usa script, make sure you make full true in options
source base-usa.R
```
## Notice
 :warning: This code shouldn't be used for inferring any estimates. It is mainly for benchmarking purposes. For actual code look [here](https://github.com/ImperialCollegeLondon/covid19mode)

:warning: This code is released with no support. We try our best to look at issues and pull request but can't help people with setup most of the time. We have docker images and conda environment file to make it easy for you to get started with the setup, any other approach assumes user can handle their computing environments appropriately.

:warning: This model is in active development and so parameter name and behaviours, and output file formats will change without notice.

:warning: As with any mathematical model, it is easy to misconfigure inputs and therefore get meaningless outputs. The development team only endorses outputs it has itself generated.
