# raptor-galaxy

This readme explains how to add a tool to bioconda and galaxy (by using planemo).

## Bioconda
Bioconda has so called recipes for every available tool which are saved in their git repo https://github.com/bioconda/bioconda-recipes .
To get your tool into bioconda you have to create a PR with a recipe for your tool.
You can use the recipe of raptor as a template: https://github.com/bioconda/bioconda-recipes/tree/master/recipes/raptor

## Install and setup of planemo
The next 3 steps will install planemo into the directory `.venv-planemo`.
```
$ virtualenv .venv-planemo; source .venv-planemo/bin/activate
$ pip install "pip>=7"
$ pip install planemo
```
Each time you open a new console and want to use planemo, it is needed to activate the planemo environment:
```
$ source .venv-planemo/bin/activate
```
For more information on how to install planemo checkout the project https://github.com/galaxyproject/planemo .


## Account on the toolshed
There exists the normal toolshed at https://toolshed.g2.bx.psu.edu/ and the test toolshed at https://testtoolshed.g2.bx.psu.edu/.
The steps for both are the same. We will show here how to use the testtoolshed.
- create an account on https://testtoolshed.g2.bx.psu.edu
- run `$ planemo config_init` to create a planemo config
- insert shed_username and api key into ~/.planemo.yml

## Creating a new description for a tool
Lets assume you want to add a new subcommand to raptor called `newsubcommand`.
- copy raptor-build.xml to raptor-newsubcommand.xml
- adjust the xml file to own needs
- run `$ planemo test raptor-newsubcommand.xml` to check tests
- run `$ planemo lint raptor-newsubcommand.xml` to lint

## Publishing your tool
Make sure to bump the version number of the tool for every new release (see xml tag  <tool version="..."> ).
- run `$ planemo shed_update --shed_target testtoolshed path/to/this/repo` to publish the tools of this repository
