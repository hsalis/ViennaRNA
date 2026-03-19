## Here, we define several environment variables required by our unit tests
##

# include current directory into our Python path variable
PYTHONPATH='/Users/howardsalis/ViennaRNA/tests:/Users/howardsalis/ViennaRNA/tests:${PYTHONPATH}'

export PYTHONPATH

# include path to the built executables to check their functionality later on
PATH=/Users/howardsalis/ViennaRNA/src/bin:${PATH}

export PATH

# the path to various data files required to run our tests
export DATADIR=/Users/howardsalis/ViennaRNA/tests/data

# current ViennaRNA Package version
export CURRENT_VERSION=2.7.2

# set DIFF variable
export DIFF=/usr/bin/diff

# set results directories
export RNAFOLD_RESULTSDIR=/Users/howardsalis/ViennaRNA/tests/RNAfold/results
export RNAALIFOLD_RESULTSDIR=/Users/howardsalis/ViennaRNA/tests/RNAalifold/results
export RNACOFOLD_RESULTSDIR=/Users/howardsalis/ViennaRNA/tests/RNAcofold/results

# misc/ directory
export MISC_DIR=/Users/howardsalis/ViennaRNA/misc
