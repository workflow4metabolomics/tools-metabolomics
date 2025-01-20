#!perl

## script  : XXX.pl

## Notes :
#	-> manage score sorting : Cleaned_pcGroups done but not in outputs !

#=============================================================================
#                              Included modules and versions
#=============================================================================
## Perl modules
use strict ;
use warnings ;
use Carp qw (cluck croak carp) ;

use threads;
use threads::shared;
use Thread::Queue;

use Data::Dumper ;
use Getopt::Long ;
use POSIX ;
use FindBin ; ## Allows you to locate the directory of original perl script

## Specific Perl Modules (PFEM)
use lib $FindBin::Bin ;
my $binPath = $FindBin::Bin ;
use lib::csv  qw( :ALL ) ;
use lib::conf  qw( :ALL ) ;

## Dedicate Perl Modules (Home made...)
use lib::massbank_api qw( :ALL ) ;
use lib::threader qw(:ALL) ;
use lib::mapper qw(:ALL) ;
use lib::writter qw(:ALL) ;
use lib::massbank_parser qw(:ALL) ;



## Initialized values
my ($help, $mzs_file, $col_mz, $col_int, $col_pcgroup, $line_header ) = ( undef, undef, undef, undef, undef,undef, undef ) ;
my ($server, $ion_mode, $score_threshold, $instruments, $max, $unit, $tol, $cutoff) = ( undef, undef, undef, undef, undef, undef, undef ) ;
my ($output_json, $output_tabular, $output_xlsx, $output_html ) = ( undef, undef, undef, undef ) ;

## Local values ONLY FOR TEST :
#my $server = 'JP' ;
#my $threading_threshold = 6 ;

#=============================================================================
#                                Manage EXCEPTIONS
#=============================================================================
&GetOptions ( 	"help|h"     		=> \$help,       # HELP
				"masses:s"			=> \$mzs_file,
				"col_mz:i"			=> \$col_mz,
				"col_int:i"			=> \$col_int, ## optionnal
				"col_pcgroup:i"		=> \$col_pcgroup,
				"lineheader:i"		=> \$line_header,
				"mode:s"			=> \$ion_mode, 
				"score_threshold:f"	=> \$score_threshold, 
				"instruments:s"		=> \$instruments, # advanced -> to transform into string with comma => done !
				"max:i"				=> \$max, # advanced
				"unit:s"			=> \$unit, # advanced
				"tolerance:f"		=> \$tol, 
				"cutoff:i"			=> \$cutoff, # advanced : intensity cutoff
				"server:s"			=> \$server, ## by default JP and # advanced
				"output_json:s"		=> \$output_json,
				"output_xlsx:s"		=> \$output_xlsx,
				"output_tabular:s"	=> \$output_tabular,
				"output_html:s"		=> \$output_html,
            ) ;
         
## if you put the option -help or -h function help is started
if ( defined($help) ){ &help ; }

#=============================================================================
#                                MAIN SCRIPT
#=============================================================================

## -------------- Conf file ------------------------ :
my ( $CONF ) = ( undef ) ;
foreach my $conf ( <$binPath/*.cfg> ) {
	my $oConf = lib::conf::new() ;
	$CONF = $oConf->as_conf($conf) ;
}

## -------------- HTML template file ------------------------ :
foreach my $html_template ( <$binPath/*.tmpl> ) { $CONF->{'HTML_TEMPLATE'} = $html_template ; }

## Main variables :
my ($pcs, $mzs, $into, $complete_rows, $pcgroups) = (undef, undef, undef, undef, undef) ;

## manage csv file containing list of masses (every thing is manage in jar)
if ( ( defined $mzs_file ) and ( $mzs_file ne "" ) and ( -e $mzs_file ) ) {
	
	## parse csv ids and masses
	my $is_header = undef ;
	my $ocsv = lib::csv->new() ;
	my $csv = $ocsv->get_csv_object( "\t" ) ;
	if ( ( defined $line_header ) and ( $line_header > 0 ) ) { $is_header = 'yes' ;    }
	$pcs = $ocsv->get_value_from_csv_multi_header( $csv, $mzs_file, $col_pcgroup, $is_header, $line_header ) ; ## retrieve pc values on csv
	$mzs = $ocsv->get_value_from_csv_multi_header( $csv, $mzs_file, $col_mz, $is_header, $line_header ) ; ## retrieve mz values on csv
	$into = $ocsv->get_value_from_csv_multi_header( $csv, $mzs_file, $col_int, $is_header, $line_header ) if ( defined $col_int ); ## retrieve into values on csv // optionnal in input files
	$complete_rows = $ocsv->parse_csv_object($csv, \$mzs_file) ; ## parse all csv for output csv build

	## manage input file with no into colunm / init into with a default value of 10
	if ( !defined $col_int ) {
		my $nb_mzs = scalar(@{$mzs}) ;
		my @intos = map {10} (0..$nb_mzs-1) ;
		my $nb_intos = scalar(@intos) ;
		if ($nb_intos == $nb_mzs) { $into = \@intos ;	}
		else { carp "A difference exists between intensity and mz values\n" }
	}
	
	## manage instruments string to array_ref
	if (defined $instruments ) {
		if ($instruments eq '') { ## in xml : can select nothing...
			$instruments = ['all'] ;
		}
		else {
			my @instruments = split(/,/, $instruments) ;
			$instruments = \@instruments ;
		}
	}
	
	
	## Build pcgroups with their features :
	my $omap = lib::mapper->new() ;
	$pcgroups = $omap->get_pcgroups($pcs, $mzs, $into ) ;
	my $pcgroup_list = $omap->get_pcgroup_list($pcs ) ;
	
#	print Dumper $pcgroups ;
	
	my $pc_num = 0 ;
	$pc_num = scalar(@{$pcgroup_list}) ;
	
	## manage a list of query pc_group dependant:
	if ($pcgroups) {
		## - - - - - - -  - - - - -  - - - -  - - - - - Multithreadind mode if pcgroups > 6 - - - - - - - - - - - - - - - - 
		if ($pc_num > $CONF->{'THREADING_THRESHOLD'}) {
			print $server."\n" ;
			print "\n------  ** ** ** Using multithreading mode ** ** ** --------\n\n" ;
			my $time_start = time ;
			
			our $NBTHREADS = $CONF->{'THREADING_THRESHOLD'} ;

#			use constant THREADS => 6 ;
			my $Qworks = Thread::Queue->new();
			my @threads = () ;
			my @queries = () ;
			my @Qresults = () ;
			
			foreach my $pc_group_id (keys %{$pcgroups}) {
				push (@queries, $pcgroups->{$pc_group_id}) if $pcgroups->{$pc_group_id} ;
			}
			
			for (1..$NBTHREADS) {
				my $oworker = lib::threader->new ;
			    push @threads, threads->create(sub { $oworker->searchSpectrumWorker($Qworks, $server, $ion_mode, $instruments, $max, $unit, $tol, $cutoff) ; } ) ;
			}
			
			$Qworks->enqueue(@queries);
			$Qworks->enqueue(undef) for 1..$NBTHREADS;
			push @Qresults, $_->join foreach @threads;

			
			my $time_end = time ;
			my $seconds = $time_end-$time_start ;
			print "\n------  Time used in multithreading mode : $seconds seconds --------\n\n" ;
			
#			print Dumper @Qresults ;
			
			## controle number of returned queries :
			my $massbank_results_num = 0 ;
			$massbank_results_num = scalar @Qresults ;
			
			if ( $massbank_results_num == $pc_num ) {
				## Map @Qresults with annotation hash : pcgroup_id in @Qresults (pcgroup2) // id in $pcgroups (pcgroup2)
				foreach my $result (@Qresults) {
					## manage annotation part
					if ($result->{'pcgroup_id'}) {
						if ($pcgroups->{$result->{'pcgroup_id'}}) {
							$pcgroups->{$result->{'pcgroup_id'}}{'annotation'} = $result ;
						}
						else { carp "Carefull : no mapping possible between massbank results and initial pcgroups data\n";}
					}
					else { carp "Carefull : no pcgroup id defined in massbank results\n"; }
					
					## manage massbank_ids part
					if ($result->{'res'}) {
						my @tmp_res = map {$_->{'id'}} @{$result->{'res'}} ;
						$pcgroups->{$result->{'pcgroup_id'}}{'massbank_ids'} = \@tmp_res ;					
					}
				}
			}
			else {
				croak "[ERROR] : problem between massbank results number and pcgroups number\n";
			}
		}
		## - - - - - - -  - - - - -  - - - -  - - - - - mono thread mode if pcgroups <= 6 - - - - - - - - - - - - - - - - 
		else {
			## connexion
#			print $server."\n" ;
			my $omassbank = lib::massbank_api->new() ;
			my $soap = $omassbank->selectMassBank($server) ;
			print "\n------  ** ** ** Using batch mode ** ** ** --------\n\n" ;
			my $time_start = time ;
			foreach my $pcgroup (keys %{$pcgroups}) {
				## searchSpectrum via SOAP
				print "Annot pcgroup n-$pcgroup\n" ;
				my $oquery = lib::massbank_api->new() ;
				my ($results, $num) = $oquery->searchSpectrum($soap, $pcgroups->{$pcgroup}{'id'}, $pcgroups->{$pcgroup}{'mzmed'}, $pcgroups->{$pcgroup}{'into'}, $ion_mode, $instruments, $max, $unit, $tol, $cutoff) ;
				$pcgroups->{$pcgroup}{'annotation'} = $results ;
	#			print Dumper $results ;
			}
			my $time_end = time ;
			my $seconds = $time_end-$time_start ;
			print "\n------  Time used in foreach mode: $seconds seconds --------\n\n" ;
		}
	}
	else {
		croak "The pcgroup object is not defined\n" ;
	}
#	print "Init pcGroups results are\n"  ;
#	print Dumper $pcgroups ;
	
} ## End of elsif "defined $mzs_file"
else {
	warn "[WARN] Can't use Massbank WS service without an existing input tabular file\n" ;
	&help ;
}

## Clean zone - use threshold on massbank entry returned score
my $omap = lib::mapper->new() ;
my $cleaned_pcgroups = $omap->filter_pcgroup_res($pcgroups, $score_threshold) ;

#print "Cleaned_pcGroups are\n" ;
#print Dumper $cleaned_pcgroups ;

## add min/max value of each mzmed in the pc_group
my $pcgroups_with_intervales = $omap->add_min_max_for_pcgroup_res($cleaned_pcgroups, $tol ) ;

#print "pcGroups_with_intervales are\n" ;
#print Dumper $pcgroups_with_intervales ;


## search in the local indexed db - - - TODO - - - 

## OR search new ones 

## get all unique Massbank Ids found
my $oids = lib::mapper->new() ;
my $all_massbank_ids = $omap->compute_ids_from_pcgroups_res($cleaned_pcgroups) ;

## get entries on the MassBank server by ID by pieces of 10
my $omapper = lib::mapper->new() ;
my $recordList = $omapper->get_massbank_records_by_chunk ($server, $all_massbank_ids, 10) ;
#print "\n\nRecords are\n" ;
#print Dumper $recordList ;
#print Dumper $all_massbank_ids ;

## foreach record - get id and peaks - create a object
my %records = ();
foreach (@$recordList) {
	## parse record handles
	my $parser = lib::massbank_parser->new() ;
	my $id = $parser->getIdFromString($_) ;
	$records{$id}{'peaks'} = $parser->getPeaksFromString($_) ;
	$records{$id}{'names'} = $parser->getChemNamesFromString($_) ;
	$records{$id}{'instrument_type'} = $parser->getInstrumentTypeFromString($_) ;
	$records{$id}{'precursor_type'} = $parser->getPrecursorTypeFromString($_) ;
	$records{$id}{'ms_type'} = $parser->getMsTypeFromString($_) ;
	$records{$id}{'formula'} = $parser->getFormulaFromString($_) ;
	$records{$id}{'exact_mz'} = $parser->getExactMzFromString($_) ;
	$records{$id}{'inchi'} = $parser->getInchiFromString($_) ;
}
#print Dumper %records ;

## Map pc_groups and records
my $well_annoted_pcGroups = $omapper->mapGroupsWithRecords($pcgroups_with_intervales, \%records) ;

#print Dumper $well_annoted_pcGroups ;

## Output writting :
my ( $massbank_matrix ) = ( undef ) ;

## XLS OUTPUT -- new format
if (  (defined $output_xlsx) and  (defined $well_annoted_pcGroups) and  (defined $mzs) and  (defined $pcs)  ) {
	my $owritter = lib::writter->new() ;
	$owritter->write_xls_skel(\$output_xlsx, $mzs, $pcs, $well_annoted_pcGroups, \%records) ;
}

## CSV OUTPUT
if (  (defined $output_tabular) and  (defined $well_annoted_pcGroups) and  (defined $pcs) and  (defined $mzs) ) {
	my $omapper = lib::mapper->new() ;
	if ( ( defined $line_header ) and ( $line_header == 1 ) ) { $massbank_matrix = $omapper->set_massbank_matrix_object('massbank', $pcs, $mzs, $well_annoted_pcGroups, \%records ) ; }
	elsif ( ( defined $line_header ) and ( $line_header == 0 ) ) { $massbank_matrix = $omapper->set_massbank_matrix_object(undef, $pcs, $mzs, $well_annoted_pcGroups, \%records ) ; }

	$massbank_matrix = $omapper->add_massbank_matrix_to_input_matrix($complete_rows, $massbank_matrix) ;
	my $owritter = lib::writter->new() ;
	$owritter->write_csv_skel(\$output_tabular, $massbank_matrix) ;
}

my $json_scalar = undef ;
## JSON OUTPUT
if (  (defined $output_json) and  (defined $well_annoted_pcGroups) and  (defined $mzs) and  (defined $pcs)  ) {
	my $omapper = lib::mapper->new() ;
	$json_scalar = $omapper->map_pc_to_generic_json($pcs, $well_annoted_pcGroups, \%records) ;
	my $owritter = lib::writter->new() ;
	$owritter->write_json_skel(\$output_json, $json_scalar) ;
}

## HTML OUTPUT -- TODO
if (  (defined $output_html) and  (defined $json_scalar)  ) {
	
#	print Dumper $json_scalar ;
	
	## Uses N mz and theirs entries per page (see config file).
	# how many pages you need with your input mz list?
	my $nb_pages_for_html_out = ceil( scalar(@{$mzs} ) / $CONF->{HTML_ENTRIES_PER_PAGE} )  ;
	
	## Search condition:
	my $search_condition = "Search params : Molecular specie = $ion_mode / delta ($unit) = $tol / Score threshold = $score_threshold and max hit = $max per pcgroup" ;
	
	my $oHtml = lib::mapper->new() ;
	my ($tbody_object) = $oHtml->set_html_tbody_object( $nb_pages_for_html_out, $CONF->{HTML_ENTRIES_PER_PAGE} ) ;
	($tbody_object) = $oHtml->add_mz_to_tbody_object($tbody_object, $CONF->{HTML_ENTRIES_PER_PAGE}, $mzs, $json_scalar) ;
	($tbody_object) = $oHtml->add_entries_to_tbody_object($tbody_object, $CONF->{HTML_ENTRIES_PER_PAGE}, $mzs, $json_scalar) ;
	
	my $oWritter = lib::writter->new() ;
	$oWritter->write_html_skel(\$output_html, $tbody_object, $nb_pages_for_html_out, $search_condition, $CONF->{'HTML_TEMPLATE'}, $CONF->{'JS_GALAXY_PATH'}, $CONF->{'CSS_GALAXY_PATH'}) ;
}
else {
	warn "[WARN] The html output file or the json iss not defined\n" ;
}







#====================================================================================
# Help subroutine called with -h option
# number of arguments : 0
# Argument(s)        :
# Return           : 1
#====================================================================================
sub help {
	print STDERR "
massbank_ws_searchspectrum.pl

# massbank_ws_searchspectrum.pl is a script to use SOAP massbank webservice and send specific queries about spectra searches. 
# Input : a list of mzs, intensities, pcgroup.
# Author : Franck Giacomoni
# Email : franck.giacomoni\@clermont.inra.fr
# Version : 1.0
# Created : 20/01/2017
USAGE :		 
		massbank_ws_searchspectrum.pl -help OR
		
		massbank_ws_searchspectrum.pl 
			-masses [name of input file] -col_id -col_mz -col_int -col_pcgroup -lineheader
			-mode [ion mode : Positive, Negative or Both ]
			-score_threshold [Ignore Massbank results with a score lower than the defined threshold]
			-instruments [array of string: all or values obtained by getInstrumentTypes method]
			-max [0 is all results or int]
			-unit [unit or ppm]
			-tolerance [Tolerance of values of m/z of peaks: 0.3 unit or 50 ppm]
			-cutoff [Ignore peaks whose intensity is not larger than the value of cutoff. Default: 50)]
			-server [name of the massbank server : EU or JP only]
			-output_json [ouput file for JSON]
			-output_xls [ouput file for XLS]
			-output_tabular [ouput file for TABULAR]
		
		";
	exit(1);
}

## END of script - F Giacomoni 

__END__

=head1 NAME

 XXX.pl -- script for

=head1 USAGE

 XXX.pl -precursors -arg1 [-arg2] 
 or XXX.pl -help

=head1 SYNOPSIS

This script manage ... 

=head1 DESCRIPTION

This main program is a ...

=over 4

=item B<function01>

=item B<function02>

=back

=head1 AUTHOR

Franck Giacomoni E<lt>franck.giacomoni@clermont.inra.frE<gt>
Yann Guitton 

=head1 LICENSE

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=head1 VERSION

version 1 : 05 / 01 / 2016

version 2 : ??

=cut