package lib::massbank_mapper_Test ;

use diagnostics; # this gives you more debugging information
use warnings;    # this warns you of bad practices
use strict;      # this prevents silly errors
use Exporter ;
use Carp ;

use Data::Dumper ;

our $VERSION = "1.0";
our @ISA = qw(Exporter);
our @EXPORT = qw( max_from_double_with_deltaTest min_from_double_with_deltaTest get_massbank_records_by_chunkTest compute_ids_from_pcgroups_resTest filter_pcgroup_resTest get_pcgroup_listTest get_pcgroupsTest set_massbank_matrix_objectTest add_massbank_matrix_to_input_matrixTest);
our %EXPORT_TAGS = ( ALL => [qw(max_from_double_with_deltaTest min_from_double_with_deltaTest get_massbank_records_by_chunkTest compute_ids_from_pcgroups_resTest filter_pcgroup_resTest get_pcgroup_listTest get_pcgroupsTest set_massbank_matrix_objectTest add_massbank_matrix_to_input_matrixTest)] );

use lib '/Users/fgiacomoni/Inra/labs/perl/galaxy_tools/massbank_ws_searchspectrum' ;
use lib::mapper qw( :ALL ) ;

sub get_pcgroupsTest {
	my ( $pcs, $mzs, $ints ) = @_;
	my $omap = lib::mapper->new() ;
	my $pcgroups = $omap->get_pcgroups($pcs, $mzs, $ints) ;
#	print Dumper $pcgroups ;
	return ($pcgroups) ;
}


sub get_pcgroup_listTest {
	my ($pcs) = @_ ; 
	my $omap = lib::mapper->new() ;
	my $pcgroup_list = $omap->get_pcgroup_list($pcs) ;
#	print Dumper $pcgroup_list ;
	return ($pcgroup_list) ;
}

sub set_massbank_matrix_objectTest {
	my ($header, $init_pcs, $pcgroups) = @_ ; 
	my $omap = lib::mapper->new() ;
	my $matrix = $omap->set_massbank_matrix_object($header, $init_pcs, $pcgroups) ;
#	print Dumper $matrix ;
	return ($matrix) ;
}

sub add_massbank_matrix_to_input_matrixTest {
	my ($input_matrix_object, $massbank_matrix_object ) = @_ ; 
	my $omap = lib::mapper->new() ;
	my $matrix = $omap->add_massbank_matrix_to_input_matrix($input_matrix_object, $massbank_matrix_object ) ;
#	print Dumper $matrix ;
	return ($matrix) ;
}

## SUB TEST for filter_pcgroup_res
sub filter_pcgroup_resTest {
    # get values
    my ( $pcgroups, $threshold ) = @_;
    my $cleaned_pcgroups = () ;
    
    
    my $omap = lib::mapper->new() ;
    $cleaned_pcgroups = $omap->filter_pcgroup_res($pcgroups, $threshold) ;
    
    return($cleaned_pcgroups) ;
}
## End SUB

## SUB TEST for 
sub compute_ids_from_pcgroups_resTest {
    # get values
    my ( $pcgroups ) = @_;
    my $ids = () ;
    
    my $var2 = lib::mapper->new() ;
    $ids = $var2->compute_ids_from_pcgroups_res($pcgroups) ;
    
    return($ids) ;
}
## End SUB

## ## SUB TEST for 
sub min_from_double_with_deltaTest {
    # get values
    my ( $double, $type, $delta ) = @_;
    
    my $var2 = lib::mapper->new() ;
    my ($min, $max) = $var2->min_and_max_from_double_with_delta($double, $type, $delta) ;
#    print $min ;
    return($min) ;
}
## End SUB

## ## SUB TEST for 
sub max_from_double_with_deltaTest {
    # get values
    my ( $double, $type, $delta ) = @_;
    
    my $var2 = lib::mapper->new() ;
    my ($min, $max) = $var2->min_and_max_from_double_with_delta($double, $type, $delta) ;
#    print $max ;
    return($max) ;
}
## End SUB

## SUB TEST for 
sub get_massbank_records_by_chunkTest {
    # get values
    my ( $server, $ids, $chunk_size ) = @_;
    
    my $var2 = lib::mapper->new() ;
    my $records = $var2->get_massbank_records_by_chunk ($server, $ids, $chunk_size) ;
    
#    print Dumper $records ;
    
    return($records) ;
}
## End SUB


1 ;