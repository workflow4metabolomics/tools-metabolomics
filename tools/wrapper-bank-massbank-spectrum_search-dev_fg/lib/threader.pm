package lib::threader ;

use strict;
use warnings ;
use Exporter ;
use Carp ;
use threads;
use threads::shared;
use Thread::Queue;
use diagnostics ;
use Data::Dumper ;
use Carp qw (cluck croak carp) ;
use LWP::UserAgent;
use LWP::Simple ; ## Lib de protocole HTTP de download
use SOAP::Lite + trace => qw(fault); ## SOAP for web service version 0.67
import SOAP::Data qw(name);

use Data::Dumper ;

use vars qw($VERSION @ISA @EXPORT %EXPORT_TAGS);

our $VERSION = "1.0";
our @ISA = qw(Exporter);
our @EXPORT = qw(threading_getRecordInfo);
our %EXPORT_TAGS = ( ALL => [qw(threading_getRecordInfo )] );

=head1 NAME

My::Module - An example module

=head1 SYNOPSIS

    use My::Module;
    my $object = My::Module->new();
    print $object->as_string;

=head1 DESCRIPTION

This module does not really exist, it
was made for the sole purpose of
demonstrating how POD works.

=head1 METHODS

Methods are :

=head2 METHOD new

	## Description : new
	## Input : $self
	## Ouput : bless $self ;
	## Usage : new() ;

=cut

sub new {
    ## Variables
    my $self={};
    bless($self) ;
    return $self ;
}
### END of SUB

=head2 METHOD searchSpectrumWorker

	## Description : work with searchSpectrum method in threading mode
	## Input : $Qworks
	## Output : $results
	## Usage : my ( $results ) = searchSpectrumWorker( $Qworks ) ;
	
=cut
## START of SUB
sub searchSpectrumWorker {
	my $self = shift;
	my ($Qworks, $server, $ion_mode, $instruments, $max, $unit, $tol, $cutoff) = @_ ;
	my @results = () ;
	my @fake = () ;

	my $omassbank = lib::massbank_api->new() ;
	my $soap = $omassbank->selectMassBank($server) ;
	
	while(my $pcgroup = $Qworks->dequeue) {
#		print Dumper $pcgroup ;
	    my $oquery= lib::massbank_api->new() ;
		my ($result, $num) = eval{$oquery->searchSpectrum($soap, $pcgroup->{'id'}, $pcgroup->{'mzmed'}, $pcgroup->{'into'}, $ion_mode, $instruments, $max, $unit, $tol, $cutoff) ; } or die;
#		print "The query send to massbank return $num entries...\n" ;
#		print Dumper $result ;

		if (!defined $num) { $num = 0  ; }
		if ($num >= 0 ) {
			push @results, $result ;
		}
		else {
			push @fake, $pcgroup ;
		}
    }
	return (@results) ;
}
## END of SUB

   
=head2 METHOD threading_getRecordInfo

	## Description : prepare parallel threads - DEPRECATED
	## Input : $soap, $list
	## Output : $results
	## Usage : my ( $results ) = threading_getRecordInfo( $soap, $list ) ;
	
=cut
## START of SUB
sub threading_getRecordInfo {
	## Retrieve Values
    my $self = shift ;
    my ( $osoap, $list ) = @_ ;
    
    my @results = () ;
    my $i = 0 ; # position in the ids list
    
    if ( ( defined $list ) ) {
    	
    	my $oquery = lib::massbank_api->new() ;
    	
    	for (my $i = 0; $i < (scalar @{$list}); $i++ ) {
    		my $thr = threads->create( sub { $oquery->getRecordInfo($osoap, $list->[$i]) } ) ;
    		push ( @results, $thr->join )  if $list->[$i] ;
    	}
    }
    else {
    	warn "Your input list of ids is undefined\n" ;
    }
    return(\@results) ;
}
## END of SUB


=head2 METHOD threading_searchSpectrum

	## Description : prepare parallel threads - DEPRECATED
	## Input : $soap, 
	## Output : $results
	## Usage : my ( $results ) = threading_searchSpectrum( $soap,  ) ;
	
=cut
## START of SUB
sub threading_searchSpectrum {
	
	## http://www.perlmonks.org/?node_id=735923
	## http://www.nntp.perl.org/group/perl.ithreads/2003/05/msg696.html
	## http://stackoverflow.com/questions/15222480/web-service-using-perl-wsdl-and-multi-threading-does-not-working
	## Retrieve Values
    my $self = shift ;
    my ( $osoap, $pcgroup_list, $pcgroups, $ion_mode, $instruments, $max, $unit, $tol, $cutoff ) = @_ ;
    
    my @results = () ;
    my $n = 6 ; # position in the ids list
    
    if ( ( defined $pcgroups ) ) {
    	
    	print Dumper $pcgroups ;
    	
    	my $oquery = lib::massbank_api->new() ;
    	
    	foreach my $pc (@{$pcgroup_list}) {
    		
    		print "\t---> Create a thread for pcgroup n-$pc\n" ;
    		
    		my $thr = threads->create(
    			sub { 
    					$oquery->searchSpectrum($osoap, $pcgroups->{$pc}{'mzmed'}, $pcgroups->{$pc}{'into'}, $ion_mode, $instruments, $max, $unit, $tol, $cutoff) ;
    			} 
    		) ;
    		push ( @results, $thr->join ) ;
    	} ## end foreach
    }
    else {
    	warn "Your input list of ids is undefined\n" ;
    }
    return(\@results) ;
}
## END of SUB


=head2 METHOD thread_and_queue_searchSpectrum

	## Description : prepare parallel and queuing threads - DEPRECATED
	## Input : $soap, 
	## Output : $results
	## Usage : my ( $results ) = thread_and_queue_searchSpectrum( $soap,  ) ;
	
=cut
## START of SUB
sub thread_and_queue_searchSpectrum {
#    my $self = shift;
#    my ( ) = @_;
#    
#    our $THREADS = 10;
#	my $Qwork = new Thread::Queue;
#	my $Qresults = new Thread::Queue;
#	
#	## Create the pool of workers
#	my @pool = map{
#	    threads->create( \&worker, $Qwork, $Qresults )
#	} 1 .. $THREADS;
#	
#	## Get the work items (from somewhere)
#	## and queue them up for the workers
#	while( my $workItem = getWorkItems() ) {
#	    $Qwork->enqueue( $workItem );
#	}
#	
#	## Tell the workers there are no more work items
#	$Qwork->enqueue( (undef) x $THREADS );
#	
#	## Process the results as they become available
#	## until all the workers say they are finished.
#	for ( 1 .. $THREADS ) {
#	    while( my $result = $Qresults->dequeue ) {
#	
#	        ## Do something with the result ##
#	        print $result;
#	    }
#	}
#	
#	## Clean up the threads
#	$_->join for @pool;
}
## END of SUB



1 ;


__END__

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

 perldoc XXX.pm

=head1 Exports

=over 4

=item :ALL is ...

=back

=head1 AUTHOR

Franck Giacomoni E<lt>franck.giacomoni@clermont.inra.frE<gt>

=head1 LICENSE

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=head1 VERSION

version 1 : xx / xx / 201x

version 2 : ??

=cut