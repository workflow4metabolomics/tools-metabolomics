#!/bin/bash
# vi: fdm=marker

# Constants {{{1
################################################################

SCRIPT_PATH=$(dirname $BASH_SOURCE)
SCRIPT_NAME=$(basename $BASH_SOURCE)
MTBLSDWNLD=$SCRIPT_PATH/../mtbls-dwnld
ASPERA_PUBLIC_TOKEN=Xz68YfDe
[[ -z $(which ascp) ]] || ASCP=ascp
WGET=wget
ASPERA_N_TRIES=10

# Test wget whole study {{{1
################################################################

test_wget_whole_study() {

local study=MTBLS164

# Remove previously downloaded directory
rm -rf $study

# Download
expect_success $MTBLSDWNLD -g $study || return 1

# Test
expect_folder "$study" || return 1
expect_file "$study/i_Investigation.txt" || return 1
expect_files_in_folder "$study" '^._.*\.t.*$' || return 1
expect_other_files_in_folder "$study" '^._.*\.t.*$' || return 1
}

# Test wget zipped study {{{1
################################################################

test_wget_zipped_study() {

	local study=MTBLS164

	# Remove previously downloaded zip
	rm -f $study.zip

	# Download
	expect_success $MTBLSDWNLD -gMc $study || return 1

	# Test
	expect_non_empty_file "$study.zip" || return 1
}

# Test wget metadata only {{{1
################################################################

test_wget_metadata_only() {

	local study=MTBLS164

	# Remove previously downloaded directory
	rm -rf $study

	# Download
	expect_success $MTBLSDWNLD -gM $study || return 1

	# Test
	expect_folder "$study" || return 1
	expect_file "$study/i_Investigation.txt" || return 1
	expect_no_other_files_in_folder "$study" '^._.*\.t.*$' || return 1
}

# Get private study path {{{1
################################################################

get_private_study_path() {

	local method=$1
	local study=$2

	# Get path
	path=TEST_MTBLS_DWNLD_${method}_${study}_PRIVATE_PATH
	path=$(echo $path | tr '[:lower:]' '[:upper:]')
	expect_str_not_null "${!path}" "No path to test private download of $study with method $method. Please set environment variable $path." >/dev/null || return 1

	echo "${!path}"
}

# Get private study token {{{1
################################################################

get_private_study_token() {

	local method=$1
	local study=$2

	# Get token
	token=TEST_MTBLS_DWNLD_${method}_${study}_PRIVATE_TOKEN
	token=$(echo $token | tr '[:lower:]' '[:upper:]')
	expect_str_not_null "${!token}" "No token to test private download of $study with method $method. Please set environment variable $token." >/dev/null || return 1

	echo "${!token}"
}

# Test wget private study {{{1
################################################################

test_wget_private_study() {

	local study=MTBLS1

	# Remove previously downloaded directory
	rm -rf $study

	# Get study info
	token=$(get_private_study_token wget $study)
	path=$(get_private_study_path wget $study)

	if [[ -n $token && -n $path ]] ; then
		# Download
		expect_success $MTBLSDWNLD -gpM -t $token $path || return 1

		# Test
		expect_folder "$study" || return 1
		expect_file "$study/i_Investigation.txt" || return 1
	fi
}

# Test ascp whole study {{{1
################################################################

test_ascp_whole_study() {

	local study=MTBLS164

	# Remove previously downloaded directory
	rm -rf $study

	# Download
	expect_success_after_n_tries $ASPERA_N_TRIES $MTBLSDWNLD -ag -t "$ASPERA_PUBLIC_TOKEN" $study || return 1

	# Test
	expect_folder "$study" || return 1
	expect_file "$study/i_Investigation.txt" || return 1
	expect_files_in_folder "$study" '^._.*\.t.*$' || return 1
	expect_other_files_in_folder "$study" '^._.*\.t.*$' || return 1
}

# Test ascp default key {{{1
################################################################

test_ascp_default_key() {

	local study=MTBLS164

	# Remove previously downloaded directory
	rm -rf $study

	# Download
	expect_success_after_n_tries $ASPERA_N_TRIES $MTBLSDWNLD -aMg $study || return 1

	# Test
	expect_folder "$study" || return 1
	expect_file "$study/i_Investigation.txt" || return 1
	expect_files_in_folder "$study" '^._.*\.t.*$' || return 1
	expect_other_files_in_folder "$study" '^._.*\.t.*$' || return 1
}

# Test ascp metadata only {{{1
################################################################

test_ascp_metadata_only() {

	local study=MTBLS1

	# Remove previously downloaded directory
	rm -rf $study

	# Download
	expect_success_after_n_tries $ASPERA_N_TRIES $MTBLSDWNLD -agM $study || return 1

	# Test
	expect_folder "$study" || return 1
	expect_file "$study/i_Investigation.txt" || return 1
	expect_no_other_files_in_tree "$study" '^._.*\.t.*$' || return 1
}

# Test ascp private study {{{1
################################################################

test_ascp_private_study() {

	local study=MTBLS1

	# Remove previously downloaded directory
	rm -rf $study

	# Get study info
	token=$(get_private_study_token ascp $study)
	path=$(get_private_study_path ascp $study)

	if [[ -n $token && -n $path ]] ; then

		# Download
		expect_success_after_n_tries $ASPERA_N_TRIES $MTBLSDWNLD -agpM -t $token -o $study $path || return 1

		# Test
		expect_folder "$study" || return 1
		expect_file "$study/i_Investigation.txt" || return 1
	fi
}

# Test wget factor slicing {{{1
################################################################

test_wget_factor_slicing() {

	local study=MTBLS174

	# Remove previously downloaded directory
	rm -rf ${study}

	# Check sliced study
	expect_success $MTBLSDWNLD -g -f Age=33 $study || return 1
	expect_folder "$study" || return 1
	expect_file "$study/i_Investigation.txt" || return 1
	expect_files_in_folder "$study" '^._.*\.t.*$' || return 1
	expect_other_files_in_tree "$study" '^._.*\.t.*$' || return 1
	nb_files=$(find $study -type f | wc -l)
	nb_files_expected=5
	expect_num_eq $nb_files $nb_files_expected "The sliced study should contain $nb_files_expected files. Found $nb_files_sliced_study." || return 1
}

# Test wget temp in relative output {{{1
################################################################

test_wget_temp_in_rel_output() {

	local study=MTBLS164
	local output_dir=MTBLS164_output

	# Remove previous folders
	rm -rf "$study" "$output_dir"

	expect_success $MTBLSDWNLD -g -T -o "$output_dir" "$study" || return 1
	expect_failure test -e "$study" || return 1
	expect_folder "$output_dir" || return 1
	expect_file "$output_dir/i_Investigation.txt" || return 1
	expect_failure test -e "$output_dir/$study" || return 1
}

# Test wget temp in absolute output {{{1
################################################################

test_wget_temp_in_abs_output() {

	local study=MTBLS164
	local output_dir=$(mktemp -t ${SCRIPT_NAME}_MTBLS164_output.XXXXXX)

	# Remove previous folders
	rm -rf "$study" "$output_dir"

	expect_success $MTBLSDWNLD -g -T -o "$output_dir" "$study" || return 1
	expect_failure test -e "$study" || return 1
	expect_folder "$output_dir" || return 1
	expect_file "$output_dir/i_Investigation.txt" || return 1
	expect_failure test -e "$output_dir/$study" || return 1
}

# Test ascp temp in relative output {{{1
################################################################

test_ascp_temp_in_rel_output() {

	local study=MTBLS1
	local output_dir=MTBLS1_output

	# Remove previous folders
	rm -rf "$study" "$output_dir"

	expect_success_after_n_tries $ASPERA_N_TRIES $MTBLSDWNLD -a -M -g -T -o "$output_dir" "$study" || return 1
	expect_failure test -e "$study" || return 1
	expect_folder "$output_dir" || return 1
	expect_file "$output_dir/i_Investigation.txt" || return 1
	expect_failure test -e "$output_dir/$study" || return 1
}

# Test ascp temp in absolute output {{{1
################################################################

test_ascp_temp_in_abs_output() {

	local study=MTBLS164
	local output_dir=$(mktemp -t ${SCRIPT_NAME}_MTBLS164_output.XXXXXX)

	# Remove previous folders
	rm -rf "$study" "$output_dir"

	expect_success_after_n_tries $ASPERA_N_TRIES $MTBLSDWNLD -a -M -g -T -o "$output_dir" "$study" || return 1
	expect_failure test -e "$study" || return 1
	expect_folder "$output_dir" || return 1
	expect_file "$output_dir/i_Investigation.txt" || return 1
	expect_failure test -e "$output_dir/$study" || return 1
}

# Main {{{1
################################################################

# Small studies to test:        Data files
# MTBLS164    2.30MB            one .xlsx and one .json
# MTBLS174    3.51MB            .zip files and one .json
# MTBLS116    5.72MB            .zip files and one .json
# MTBLS77     6.70MB            .zip files, one .docx and one .xlsx
# MTBLS342    8.09MB            .fid files, one .json
# MTBLS361    23.32MB           .zip files, one .xlsx and one .json
# MTBLS566    34.10MB           2 .mzML files.
# MTBLS1      148.82MB          .zip files and .nmrML files
# MTBLS321    492.24MB          mzXML
# MTBLS213    384.47MB          mzXML
# MTBLS433    204.41MB          mzML
# MTBLS381    776.72MB          mzData
# MTBLS10     180.49MB          CDF

# wget tests
test_context "Testing mtbls-dwnld with wget"
test_that "Download of whole study with wget works correctly." test_wget_whole_study
test_that "Download of study as a zip works correctly." test_wget_zipped_study
test_that "Download of metadata only with wget works correctly." test_wget_metadata_only
test_that "Download of private study with wget works correctly." test_wget_private_study
test_that "Factor slicing works." test_wget_factor_slicing
test_that "Download with all temporary files written into relative output directory works correctly." test_wget_temp_in_rel_output
test_that "Download with all temporary files written into absolute output directory works correctly." test_wget_temp_in_abs_output

# aspera test
if [[ -n $ASCP ]] ; then
	test_context "Testing mtbls-dwnld with Aspera ascp"
	test_that "Download of whole study with ascp works correctly." test_ascp_whole_study
	test_that "Download of study using default key with ascp works correctly." test_ascp_default_key
	test_that "Download of metadata only with ascp works correctly." test_ascp_metadata_only
	test_that "Download with all temporary files written into relative output directory works correctly." test_ascp_temp_in_rel_output
	test_that "Download with all temporary files written into absolute output directory works correctly." test_ascp_temp_in_abs_output
	test_that "Download of private study with ascp works correctly." test_ascp_private_study
fi
