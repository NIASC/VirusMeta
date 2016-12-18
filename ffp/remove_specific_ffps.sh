#!/bin/sh

#removes specified .ffp files from merged_ffp directrory generated by step2. Run before step 3

export merged_ffp=$1
export remove_ffp_list_file=$2

#And remove all block files
awk -v working_dir=$merged_ffp '{print "rm  "working_dir"/"$1".ffp"}' $remove_ffp_list_file > $merged_ffp/remove_specific_ffp.sh
chmod +x $merged_ffp/remove_specific_ffp.sh
$merged_ffp/remove_specific_ffp.sh

