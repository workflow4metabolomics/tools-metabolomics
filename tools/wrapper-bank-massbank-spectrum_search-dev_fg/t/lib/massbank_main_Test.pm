package lib::massbank_main_Test ;

use diagnostics; # this gives you more debugging information
use warnings;    # this warns you of bad practices
use strict;      # this prevents silly errors
use Exporter ;
use Carp ;

use Data::Dumper ;

our $VERSION = "1.0";
our @ISA = qw(Exporter);
our @EXPORT = qw( run_main_massbank_pl);
our %EXPORT_TAGS = ( ALL => [qw(run_main_massbank_pl)] );

#use lib '/Users/fgiacomoni/Inra/labs/perl/galaxy_tools/massbank_ws_searchspectrum' ;
#use lib::mapper qw( :ALL ) ;

sub run_main_massbank_pl {

	my ($input_file, $col_mz, $col_pcgroup, $col_int, $lineheader, $score_threshold, $mode, $instruments, $max, $unit, $tol, $cutoff, $server, $json, $csv, $xls) = @_ ;
	
	my $msg = `perl /Users/fgiacomoni/Inra/labs/perl/galaxy_tools/massbank_ws_searchspectrum/massbank_ws_searchspectrum.pl -masses $input_file -col_mz $col_mz -col_pcgroup $col_pcgroup -col_int $col_int -lineheader $lineheader -mode $mode -score_threshold $score_threshold -instruments $instruments -max $max -unit $unit -tolerance $tol -cutoff $cutoff -server $server -json $json -csv $csv -xls $xls`;
	print $msg ;
#-masses /Users/fgiacomoni/Inra/labs/tests/massbank_V02/pcgrp_annot.tsv
#-col_mz 1
#-col_pcgroup 14
#-col_int 7
#-lineheader 1
#-score_threshold 0.5
#-mode Positive
#-instruments all
#-max 2
#-unit unit
#-tolerance 0.3
#-cutoff 5
#-server JP
#-json /Users/fgiacomoni/Inra/labs/tests/massbank_V02/out.json
#-csv /Users/fgiacomoni/Inra/labs/tests/massbank_V02/out.csv
#-xls /Users/fgiacomoni/Inra/labs/tests/massbank_V02/out.xls

}


1 ;