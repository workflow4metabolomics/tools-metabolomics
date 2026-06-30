#!/bin/bash

# set path
scriptdir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# run test
Rscript $scriptdir/filter_wrap.R dataMatrix_in "$scriptdir/test-data/in_DM1.tabular" sampleMetadata_in "$scriptdir/test-data/in_SM1.tabular" variableMetadata_in "$scriptdir/test-data/in_VM1.tabular" Numeric "TRUE" num_file "variable" parm_col "rt" Interval "lower" low_value "1.2" Factors "TRUE" qual_file "sample" factor_col "Time" factors_value "3" dataMatrix_out "$scriptdir/out_DM1.tabular" sampleMetadata_out "$scriptdir/out_SM1.tabular" variableMetadata_out "$scriptdir/out_VM1.tabular"

# test diff
diff $scriptdir/out_DM1.tabular $scriptdir/test-data/out_DM1.tabular || exit 2
diff $scriptdir/out_SM1.tabular $scriptdir/test-data/out_SM1.tabular || exit 2
diff $scriptdir/out_VM1.tabular $scriptdir/test-data/out_VM1.tabular || exit 2