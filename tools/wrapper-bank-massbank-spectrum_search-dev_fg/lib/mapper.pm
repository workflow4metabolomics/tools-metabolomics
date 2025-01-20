package lib::mapper ;

use strict;
use warnings ;
use Exporter ;
use Carp ;

use Data::Dumper ;

use vars qw($VERSION @ISA @EXPORT %EXPORT_TAGS);

our $VERSION = "1.0";
our @ISA = qw(Exporter);
our @EXPORT = qw( add_min_max_for_pcgroup_res get_massbank_records_by_chunk compute_ids_from_pcgroups_res filter_pcgroup_res get_pcgroup_list get_pcgroups set_massbank_matrix_object add_massbank_matrix_to_input_matrix map_pc_to_generic_json set_html_tbody_object add_mz_to_tbody_object add_entries_to_tbody_object);
our %EXPORT_TAGS = ( ALL => [qw( add_min_max_for_pcgroup_res get_massbank_records_by_chunk compute_ids_from_pcgroups_res filter_pcgroup_res get_pcgroup_list get_pcgroups set_massbank_matrix_object add_massbank_matrix_to_input_matrix map_pc_to_generic_json set_html_tbody_object add_mz_to_tbody_object add_entries_to_tbody_object)] );

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

=head2 METHOD get_pcgroups

	## Description : get and prepare pcgroup features (mzs, into, names) from input cvs parser
	## Input : $pcs, $mzs, $ints, $names
	## Output : $pcgroups
	## Usage : my ( $pcgroups ) = get_pcgroups( $pcs, $mzs, $ints, $names ) ;
	
=cut
## START of SUB
sub get_pcgroups {
    my $self = shift;
    my ( $pcs, $mzs, $ints ) = @_;
    
    my %pcgroups = () ;
    my $i = 0 ;
    
    ## Warn diff matrix dimension :
    my $num_pcs = scalar(@{$pcs}) ;
    my $num_mzs = scalar(@{$mzs}) ;
    my $num_ints = scalar(@{$ints}) ;
    
    if ( ($num_pcs == $num_mzs ) and ( $num_mzs == $num_ints ) ) {
		my @pcs = @{$pcs} ;
    	
    	foreach my $pc (@{$pcs}) {
	    	
	    	if ( ! $pcgroups{$pc} ) { $pcgroups{$pc}->{'id'} = $pc ;	$pcgroups{$pc}->{'annotation'} = {} ; $pcgroups{$pc}->{'massbank_ids'} = [] ; }
	    	
	    	push (@{$pcgroups{$pc}->{'mzmed'}}, $mzs->[$i]) if ($mzs->[$i]) ; ## map mzs by pcgroup

	    	if ($ints->[$i] > 0 ) { 	push (@{$pcgroups{$pc}->{'into'}}, $ints->[$i])  ; ## map into by pcgroup
	    	}
	    	elsif ($ints->[$i] == 0) {
	    		push (@{$pcgroups{$pc}->{'into'}}, $ints->[$i])  ; ## map into by pcgroup even value is 0
	    	}
	    	else {
	    		warn "Undefined value found in pcgroups array\n" ;
	    	}
	    	$i++ ;
	    }
    }
    else {
    	warn "The different ARRAYS (pcs, mzs, ints) doesn't have the same size : mapping is not possible \n!!"
    }
    return (\%pcgroups) ;
}
### END of SUB

=head2 METHOD get_pcgroup_list

	## Description : get and prepare unik pcgroup list from input cvs parsed list
	## Input : $pcs
	## Output : $list
	## Usage : my ( $list ) = get_pcgroup_list( $pcs ) ;
	
=cut
## START of SUB
sub get_pcgroup_list {
	my $self = shift;
    my ( $pcs ) = @_;
    
    my @pcgroup_list = () ;
    my $i = 0 ;
    
    my %hash = map { $_, 1 } @{$pcs} ;
 	@pcgroup_list = keys %hash;
 	@pcgroup_list = sort { $a <=> $b } @pcgroup_list ;
	
	return (\@pcgroup_list) ;
}

### END of SUB


=head2 METHOD filter_pcgroup_res

	## Description : This method filter the results returned by massbank with a user defined score threshold
	## Input : $pcgroups, $threshold
	## Output : $pcgroups
	## Usage : my ( $pcgroups ) = filter_pcgroup_res ( $pcgroups, $threshold ) ;
	
=cut
## START of SUB
sub filter_pcgroup_res {
    ## Retrieve Values
    my $self = shift ;
    my ( $pcgroups, $threshold ) = @_ ;

    my %temp = () ;
    
    if (!defined $threshold) {
    	$threshold = 0.5 ; ## default value
    }

	if ( (defined $pcgroups) and (defined $threshold) ) {
		%temp = %{$pcgroups} ;
		
		foreach my $pc (keys %temp) {
			
			if ( $temp{$pc}{'annotation'}{'num_res'} > 0 ) {
					my @filtered_annot = reverse(grep { $_->{'score'} >= $threshold if ($_->{'score'})  } @{$temp{$pc}{'annotation'}{'res'}}) ;
					my $new_num_res = scalar (@filtered_annot) ;
					my @ids = () ;
					foreach (@filtered_annot) { push (@ids, $_->{'id'} ) }
					$temp{$pc}{'annotation'}{'res'} =\@filtered_annot ;
					$temp{$pc}{'annotation'}{'num_res'}  = $new_num_res ;
					$temp{$pc}{'massbank_ids'} = \@ids ;
			}
			else {
				warn "No result found for this pcgroup $pc\n" ;
			}
		}
	} ## End IF
	else {
		warn "No pcgroup and threshold defined\n" ;
	}
    return (\%temp) ;
}
### END of SUB

=head2 METHOD add_min_max_for_pcgroup_res

	## Description : This method add min / max value for each mzmed contained in pcgroup
	## Input : $pcgroups
	## Output : $pcgroups
	## Usage : my ( $pcgroups ) = add_min_max_for_pcgroup_res ( $pcgroups ) ;
	
=cut
## START of SUB
sub add_min_max_for_pcgroup_res {
    ## Retrieve Values
    my $self = shift ;
    my ( $pcgroups, $delta ) = @_ ;

    my %temp = () ;
    
    if (!defined $delta) {
    	$delta = 0.01 ; ## default value
    }

	if ( defined $pcgroups) {
		%temp = %{$pcgroups} ;

		foreach my $pc (keys %temp) {
			my %mz_intervales = () ;
			if ( $temp{$pc}{'mzmed'} ) {
				my @temp = @{$temp{$pc}{'mzmed'}} ;
				foreach my $mz (@temp) {
					my ($min, $max) = lib::mapper::new->min_and_max_from_double_with_delta($mz, 'Da', $delta);
					$mz_intervales{$mz} = {'min' => $min, 'max' => $max } ;
				}
			}
			else {
				warn "No mzmed found for this pcgroup\n" ;
			}
			$temp{$pc}{'intervales'} = \%mz_intervales ;
			
		}
	} ## End IF
	else {
		warn "No pcgroup and threshold defined\n" ;
	}
    return (\%temp) ;
}
### END of SUB



=head2 METHOD min_and_max_from_double_with_delta

	## Description : returns the minimum and maximum double according to the delta
	## Input : \$double, \$delta_type, \$delta
	## Output : \$min, \$max
	## Usage : ($min, $max)= min_and_max_from_double_with_delta($double, $delta_type, $mz_delta) ;
	
=cut
## START of SUB
sub min_and_max_from_double_with_delta {
	## Retrieve Values
    my $self = shift ;
    my ( $double, $delta_type, $delta ) = @_ ;
    my ( $min, $max ) = ( undef, undef ) ;
    
	if ($delta_type eq 'ppm'){
		$min = $double - ($delta * 10**-6 * $double);
		$max = $double + ($delta * 10**-6 * $double) + 0.0000000001; ## it's to included the maximum value in the search
	}
	elsif ($delta_type eq 'Da'){
		$min = $double - $delta;
		$max = $double + $delta + 0.0000000001; ## it's to included the maximum value in the search
	}
	else {	croak "The double delta type '$delta_type' isn't a valid type !\n" ;	}
	
    return($min, $max) ;
}
## END of SUB


=head2 METHOD compute_ids_from_pcgroups_res

	## Description : get all ids returned by massbank with sent queries and keep only unique ones.
	## Input : $pcgroups
	## Output : $unique_ids
	## Usage : my ( $unique_ids ) = compute_ids_from_pcgroups_res ( $pcgroups ) ;
	
=cut
## START of SUB
sub compute_ids_from_pcgroups_res {
    ## Retrieve Values
    my $self = shift ;
    my ( $pcgroups ) = @_;
    my ( @ids, @unique ) = ( (), () ) ;
    
    if ( defined $pcgroups ) {
		
		foreach my $pc ( keys %{$pcgroups} ) {
			if ( $pcgroups->{$pc}{'massbank_ids'} ) {
				push (@ids , @{ $pcgroups->{$pc}{'massbank_ids'} } ) ;
			}
		}
		
		if ( ( scalar (@ids) ) > 0 ) {
#			print Dumper @ids ;
			@unique = do { my %seen; grep { !$seen{$_}++ if (defined $_) } @ids };
			@unique = sort { $a cmp $b } @unique;
		}
		else {
			@unique = () ;
		}
		
		
	}
    return (\@unique) ;
}
### END of SUB


=head2 METHOD get_massbank_records_by_chunk

	## Description : get massbank records from a complete list but send queries chunk by chunk.
	## Input : $ids, $chunk_size
	## Output : $records
	## Usage : my ( $records ) = get_massbank_records_by_chunk ( $ids, $chunk_size ) ;
	
=cut
## START of SUB
sub get_massbank_records_by_chunk {
    ## Retrieve Values
    my $self = shift ;
    my ( $server, $ids, $chunk_size ) = @_;
    my ( @records, @sent_ids ) = ( (), () ) ;
    
    my $current = 0 ;
    my $pos = 1 ; 
    my @temp_ids = () ;
    
    my $num_ids = scalar(@{$ids}) ;
#    print "The number of given massbank ids is: $num_ids\n" ;
    
    foreach my $id (@{$ids}) {
    	$current++ ;
#    	print "$id - - $current/$num_ids) - - $pos \n" ;
    	
    	if (  ($current == $num_ids) or ($pos == $chunk_size)  ) {
#    		print "Querying Massbank with...\n" ;
    		push (@temp_ids, $id) ;
    		## send query
    		my $omassbank = lib::massbank_api->new() ;
			my ($osoap) = $omassbank->selectMassBank($server) ;
			my ($records) = $omassbank->getRecordInfo($osoap, \@temp_ids) ;
			push (@records, @{$records}) ;
    		
    		@temp_ids = () ; 
    		$pos = 0 ;
    	}
    	elsif ($pos < $chunk_size) {
#    		print "store...\n";
    		push (@temp_ids, $id) ;
    		$pos ++ ;
    	}
    	else {
    		warn "Something goes wrong : out of range\n"
    	}
    	
    	
    }
    my $num_records = scalar(@records) ;
#    print "The number of received massbank records is: $num_records\n" ;
    return (\@records) ;
}
### END of SUB

=head2 METHOD set_massbank_matrix_object

	## Description : build the massbank_row under its ref form
	## Input : $header, $init_pcs, $init_mzs, $pcgroups, $records
	## Output : $massbank_matrix
	## Usage : my ( $massbank_matrix ) = set_lm_matrix_object( $header, $init_pcs, $init_mzs, $pcgroups, $records ) ;
	
=cut
## START of SUB
sub set_massbank_matrix_object {
	## Retrieve Values
    my $self = shift ;
    my ( $header, $init_pcs, $init_mzs, $pcgroups, $records ) = @_ ;
    my @massbank_matrix = () ;
	
	my $current_pos = 0 ;
    
    ## format massbank(score::name::mz::formula::adduct::id)
    if ( defined $header ) {
    	$header .= '(score::name::mz::formula::adduct::id)' ;
    	my @headers = () ;
    	push @headers, $header ;
    	push @massbank_matrix, \@headers ;
    }
			
	## foreach mz of the input file
	foreach my $mz (@{$init_mzs}) {
		
		my $nb_ids = 0 ;
		my @ids = () ;
		
		my $pc = $init_pcs->[$current_pos] ; ## get the rigth pcgroup with maz postion in list
#		print "---> Current PCGROUP is $pc\n" ;
		if ( $pcgroups->{$pc}{'enrich_annotation'}{$mz} ) {
			## get record_ids 
			my @massbank_ids = @{ $pcgroups->{$pc}{'enrich_annotation'}{$mz} } ; ## get validated ids relative to one mz
			$nb_ids = scalar (@massbank_ids) ;
#			print "- - - NB RECORDS FOR MZ $mz = $nb_ids - - STATUS => \t" ;
			my $massbank_ids_string = undef ;
			## manage empty array
			if (!defined $nb_ids) { carp "The number of massbank ids is not defined\n" ; }
			elsif ( $nb_ids > 0 ) {
				## get data from records and init_annotation
				my $index_entries = 0 ;
				foreach my $record_id (@massbank_ids) {
					my $massbank_name = $records->{$record_id}{names}[0] ;
		    		my $massbank_id = $record_id  ;
		    		my $massbank_formula = $records->{$record_id}{formula} ;
		    		my $massbank_cpd_mz = $records->{$record_id}{exact_mz}  ;
		    		my $massbank_adduct = $records->{$record_id}{precursor_type}  ;
		    		my $massbank_score = 0 ;
		    		
		    		## getting the score 
		    		my @filtered_records= @{ $pcgroups->{$pc}{'annotation'}{res} } ;
		    		foreach my $record (@filtered_records) {
		    			if ($record->{id} eq $massbank_id ) {
		    				$massbank_score = $record->{score} ;
		    				last ;
		    			}
		    			else {
		    				next ;
		    			}
		    		}
			    	
			    	## METLIN data display model
			   		## entry1= ENTRY_DELTA::ENTRY_ENTRY_NAME::ENTRY_CPD_MZ::ENTRY_FORMULA::ENTRY_ADDUCT::ENTRY_ENTRY_ID | entry2=VAR1::VAR2::VAR3::VAR4|...
			   		my $massbank_id_string = $massbank_score.'::['."$massbank_name".']::'.$massbank_cpd_mz.'::'.$massbank_formula.'::['.$massbank_adduct.']::'.$massbank_id ;
			   		
			   		# manage final pipe
			   		if ($index_entries < $nb_ids-1 ) { 	$massbank_ids_string .= $massbank_id_string.' | ' ; }
			   		else { 						   			$massbank_ids_string .= $massbank_id_string ; 	}
			   		$index_entries++;
				}
			}
			elsif ( $nb_ids == 0 ) { $massbank_ids_string = 'NONE' ; }
			else {
				$massbank_ids_string = 'NONE' ;
			}
#			print "$massbank_ids_string\n" ;
			push (@ids, $massbank_ids_string) ;
		} ## End if
		else {
			next;
		}
		$current_pos++ ;
		
		push (@massbank_matrix, \@ids) ;
	} ## End foreach mz
#	print "* * * * Start of the MATRIX: * * * *\n" ;
#    print Dumper @massbank_matrix ;
#    print "* * * * END of the MATRIX * * * *\n" ;
    return(\@massbank_matrix) ;
}
## END of SUB

=head2 METHOD add_massbank_matrix_to_input_matrix

	## Description : build a full matrix (input + lm column)
	## Input : $input_matrix_object, $massbank_matrix_object
	## Output : $output_matrix_object
	## Usage : my ( $output_matrix_object ) = add_massbank_matrix_to_input_matrix( $input_matrix_object, $massbank_matrix_object ) ;
	
=cut
## START of SUB
sub add_massbank_matrix_to_input_matrix {
	## Retrieve Values
    my $self = shift ;
    my ( $input_matrix_object, $massbank_matrix_object ) = @_ ;
    
    my @output_matrix_object = () ;
    my $index_row = 0 ;
    
    foreach my $row ( @{$input_matrix_object} ) {
    	my @init_row = @{$row} ;
    	
    	if ( $massbank_matrix_object->[$index_row] ) {
    		my $dim = scalar(@{$massbank_matrix_object->[$index_row]}) ;
    		
    		if ($dim > 1) { warn "the add method can't manage more than one column\n" ;}
    		my $lm_col =  $massbank_matrix_object->[$index_row][$dim-1] ;

   		 	push (@init_row, $lm_col) ;
	    	$index_row++ ;
    	}
    	push (@output_matrix_object, \@init_row) ;
    }
    return(\@output_matrix_object) ;
}
## END of SUB

=head2 METHOD map_res_to_generic_json

	## Description : build json structure with all massbank results
	## Input : $mzs, $pcs, $pcgroups_results
	## Output : $json_scalar
	## Usage : my ( $json_scalar ) = add_massbank_matrix_to_input_matrix( $mzs, $pcs, $pcgroups_results ) ;
	
=cut
## START of SUB
sub map_pc_to_generic_json {
    my $self = shift;
    my ( $pcs, $pcgroups, $records ) = @_ ;
    
#    print Dumper $pcgroups ;
#    print Dumper $records ;

    ## JSON DESIGN
	my %JSON = (
		QUERY => {},
		PARAM => {},
		TYPE => {}	
	) ;
	
	my %oEntry = (
		mzmed => undef,
		into  => undef,
		mzmin => undef,
		mzmax => undef,
		pcgroup => undef,
		num_res => undef,
		RECORDS => undef,
	) ;
	
	
	my %oRecord = (
		id  => undef,
		exact_mz  => undef,
		score  => undef,
		formula => undef,
		inchi => undef,
		ms_type => undef,
		precursor_type => undef,
		instrument_type => undef,
		name => undef,
		peaks => undef,
	) ;
	
	
    
    foreach my $pc (@{$pcs}) {
		
    	my $pc_res = {} ;
    	my $num_res = undef ;
    	
    	if ($pcgroups->{$pc}) {
			my $pos = 0 ;
			## foreach mz of the pcgroup
			foreach my $mz (@{ $pcgroups->{$pc}{mzmed} } ) {
				
				my %entry = %oEntry ;
				##
				if ( defined $mz ) 								{	$entry{mzmed} = $mz ; }
				if ( $pcgroups->{$pc}{intervales}{$mz} ) 		{	$entry{mzmin} = $pcgroups->{$pc}{intervales}{$mz}{min} ; }
				if ( $pcgroups->{$pc}{intervales}{$mz} ) 		{	$entry{mzmax} = $pcgroups->{$pc}{intervales}{$mz}{max} ; }
				if ( $pcgroups->{$pc}{into}[$pos] )	 			{	$entry{into} = $pcgroups->{$pc}{into}[$pos] ; }
				if ( defined $pc ) 								{	$entry{pcgroup} = $pc ; }
				## get RECORDS
				if ( $pcgroups->{$pc}{enrich_annotation}{$mz} ) 	{
					
					my @recs = @{ $pcgroups->{$pc}{enrich_annotation}{$mz} } ;
					$entry{num_res} = scalar(@recs) ;
					
					foreach my $recId (@recs) {
						
						my %record = %oRecord ;
						if ( $records->{$recId} ) 					{ 	$record{id} = $recId ; }
						if ( $records->{$recId}{exact_mz} ) 		{ 	$record{exact_mz} = $records->{$recId}{exact_mz} ; }
						if ( $records->{$recId}{formula} ) 			{ 	$record{formula} = $records->{$recId}{formula} ; }
						if ( $records->{$recId}{ms_type} ) 			{ 	$record{ms_type} = $records->{$recId}{ms_type} ; }
						if ( $records->{$recId}{precursor_type} ) 	{ 	$record{precursor_type} = $records->{$recId}{precursor_type} ; }
						if ( $records->{$recId}{instrument_type} ) 	{ 	$record{instrument_type} = $records->{$recId}{instrument_type} ; }
						if ( $records->{$recId}{names} ) 			{ 	$record{name} = $records->{$recId}{names}[0] ; }
						if ( $records->{$recId}{inchi} ) 			{ 	$record{inchi} = $records->{$recId}{inchi} ; }
						## peaks TODO...
						
						## Score / BIG SHIT /
						foreach my $record (@{ $pcgroups->{$pc}{'annotation'}{res} }) {
			    			if ($record->{id} eq $recId ) {
			    				$record{score} = $record->{score} ;
			    				last ;
			    			}
			    			else {
			    				next ;
			    			}
			    		} ## foreach record - - - for score
						$entry{RECORDS}{$recId} = \%record ;
					} ## foreach recId
				} ## end IF
				
				$JSON{QUERY}{$mz} = \%entry ;
				$pos ++ ;
			} ## End FOREACH MZ
    	}
    	else {
    		warn "The pc group $pc doesn't exist in results !" ;
    	}    	
    }
#    print Dumper %JSON ;
	return(\%JSON) ;
}
## END of SUB


=head2 METHOD mapGroupsWithRecords

	## Description : map records with pcgroups mz to adjust massbank id annotations
	## Input : $pcgroups, $records
	## Output : $pcgroups
	## Usage : my ( $var4 ) = mapGroupsWithRecords ( $$pcgroups, $records ) ;
	
=cut
## START of SUB
sub mapGroupsWithRecords {
    ## Retrieve Values
    my $self = shift ;
    my ( $pcgroups, $records ) = @_;

    my %temp = () ;
    my (%intervales, @annotation_ids) = ( (), ()  ) ;
    
    if ( ( defined $pcgroups ) and ( defined $records )  ) {
    	
		%temp = %{$pcgroups} ;
		my %unik_real_ids = () ;
		my @real_ids = () ;
		
		foreach my $pc (keys %temp) {

			if ( $temp{$pc}{'intervales'} ) { %intervales = %{$temp{$pc}{'intervales'}} ; }
			else { warn "Cant't find any intervale values\n" ; }
			if ( $temp{$pc}{'massbank_ids'} ) { @annotation_ids = @{$temp{$pc}{'massbank_ids'}} ; }
			else { warn "Cant't find any massbank id values\n" ; }
			
#			print Dumper %intervales;
#			print Dumper @annotation_ids ;
			
			## map with intervales 
			foreach my $mz (keys %intervales) {
				my @filteredIds = () ;
				my ( $min, $max ) = ( $intervales{$mz}{'min'}, $intervales{$mz}{'max'} ) ;
				
				foreach my $id (@annotation_ids) {
#					print "Analyse mzs of id: $id...\n" ;
					if ( (defined $id) and ( $records->{$id}) ) {
						
						my %currentRecord = %{$records->{$id}} ;
						
						if (scalar @{$currentRecord{'peaks'} } > 0 ) {
							## 
							foreach my $peak_mz (@{ $currentRecord{'peaks'} } ) {
								if ($peak_mz) {
									my $record_mz = $peak_mz->{'mz'} ;
									if ( ($record_mz > $min ) and ($record_mz < $max) ){
										
										if (!exists $unik_real_ids{$id}) {
											$unik_real_ids{$id} = 1 ;
											push (@filteredIds, $id) ;
		#									print "$mz - - $id\n" ;
										}
										
										
									}
									else {
										next ;
									}
								}
								else {
									warn "The mz field is not defined\n" ;
								}
							} ## foreach
						}
						else {
							warn "The record ($id) has no peak\n" ;
						}
					}
					else {
						if (defined $id) {
							warn "The id $id seems to be not present in getting records\n" ;
						}
						else {
							warn "This catched id seems to be undef in getting records\n" ;
						}
						
						next ;
					}
				} ## end foreach
				## to avoid multiple ids
#				foreach my $id (keys %unik_real_ids) {
#					push(@real_ids, $id) ;
#				}
				%unik_real_ids = () ;
#				my @temp = @real_ids ;
				my @temp = @filteredIds ;
				$temp{$pc}{'enrich_annotation'}{$mz} = \@temp ;
				@real_ids = () ;
				@filteredIds = () ;
			} ## End foreach mz
			@annotation_ids = () ;
		} ## End foreach pc
    }
    else {
    	warn"Can't find record or pcgroup data\n" ;
    }
    
    return (\%temp) ;
}
### END of SUB

=head2 METHOD set_html_tbody_object

	## Description : initializes and build the tbody object (perl array) needed to html template
	## Input : $nb_pages, $nb_items_per_page
	## Output : $tbody_object
	## Usage : my ( $tbody_object ) = set_html_tbody_object($nb_pages, $nb_items_per_page) ;
	
=cut
## START of SUB
sub set_html_tbody_object {
	my $self = shift ;
    my ( $nb_pages, $nb_items_per_page ) = @_ ;

	my ( @tbody_object ) = ( ) ;
	
	for ( my $i = 1 ; $i <= $nb_pages ; $i++ ) {
	    
	    my %pages = ( 
	    	# tbody feature
	    	PAGE_NB => $i,
	    	MASSES => [], ## end MASSES
	    ) ; ## end TBODY N
	    push (@tbody_object, \%pages) ;
	}
    return(\@tbody_object) ;
}
## END of SUB

=head2 METHOD add_mz_to_tbody_object

	## Description : initializes and build the mz object (perl array) needed to html template
	## Input : $tbody_object, $nb_items_per_page, $mz_list
	## Output : $tbody_object
	## Usage : my ( $tbody_object ) = add_mz_to_tbody_object( $tbody_object, $nb_items_per_page, $mz_list ) ;
	
=cut
## START of SUB
sub add_mz_to_tbody_object {
	my $self = shift ;
    my ( $tbody_object, $nb_items_per_page, $mz_list, $json ) = @_ ;

	my ( $current_page, $mz_index ) = ( 0, 0 ) ;
	
	foreach my $page ( @{$tbody_object} ) {
		
		my @colors = ('white', 'green') ;
		my ( $current_index, , $icolor ) = ( 0, 0 ) ;
		
		for ( my $i = 1 ; $i <= $nb_items_per_page ; $i++ ) {
			# 
			if ( $current_index > $nb_items_per_page ) { ## manage exact mz per html page
				$current_index = 0 ; 
				last ; ##
			}
			else {
				$current_index++ ;
				if ( $icolor > 1 ) { $icolor = 0 ; }
				
				if ( exists $mz_list->[$mz_index]  ) {
					
					my %mz = (
						# mass feature
						MASSES_ID_QUERY => "mz_0".sprintf("%04s", $mz_index+1 ) ,
						MASSES_MZ_QUERY => $mz_list->[$mz_index],
						MASSES_PCGROUP_QUERY => $json->{QUERY}{ $mz_list->[$mz_index] }{pcgroup} ,
						MZ_COLOR => $colors[$icolor],
						MASSES_NB => $mz_index+1,
						ENTRIES => [] ,
					) ;
					push ( @{ $tbody_object->[$current_page]{MASSES} }, \%mz ) ;
					# Html attr for mass
					$icolor++ ;
				}
			}
			$mz_index++ ;
		} ## foreach mz

		$current_page++ ;
	}
    return($tbody_object) ;
}
## END of SUB

=head2 METHOD add_entries_to_tbody_object

	## Description : initializes and build the entries object (perl array) needed to html template
	## Input : $tbody_object, $nb_items_per_page, $mz_list, $entries
	## Output : $tbody_object
	## Usage : my ( $tbody_object ) = add_entries_to_tbody_object( $tbody_object, $nb_items_per_page, $mz_list, $entries ) ;
	
=cut
## START of SUB
sub add_entries_to_tbody_object {
	## Retrieve Values
    my $self = shift ;
    my ( $tbody_object, $nb_items_per_page, $mz_list, $JSON ) = @_ ;
    
    my $index_page = 0 ;
    my $index_mz_continous = 0 ;
    
    foreach my $page (@{$tbody_object}) {
    	
    	my $index_mz = 0 ;
    	
    	foreach my $mz (@{ $tbody_object->[$index_page]{MASSES} }) {
    		my $index_entry = 0 ;
    		my $check_noentry = 0 ;
    		my @toSort = () ;
    		
    		foreach my $record (keys %{ $JSON->{QUERY}{$mz->{MASSES_MZ_QUERY}}{RECORDS} }) {
    			$check_noentry ++ ;

    				my %entry = (
		    			ENTRY_COLOR => $tbody_object->[$index_page]{MASSES}[$index_mz]{MZ_COLOR},
		    			ENTRY_ENTRY_NAME => $JSON->{QUERY}{$mz->{MASSES_MZ_QUERY}}{RECORDS}{$record}{name}, 	
		   				ENTRY_ENTRY_ID => $JSON->{QUERY}{$mz->{MASSES_MZ_QUERY}}{RECORDS}{$record}{id}, 
		   				ENTRY_ENTRY_ID2 => $JSON->{QUERY}{$mz->{MASSES_MZ_QUERY}}{RECORDS}{$record}{id},
						ENTRY_FORMULA => $JSON->{QUERY}{$mz->{MASSES_MZ_QUERY}}{RECORDS}{$record}{formula},
						ENTRY_CPD_MZ => $JSON->{QUERY}{$mz->{MASSES_MZ_QUERY}}{RECORDS}{$record}{exact_mz},
						ENTRY_MS_TYPE => $JSON->{QUERY}{$mz->{MASSES_MZ_QUERY}}{RECORDS}{$record}{ms_type},
						ENTRY_PRECURSOR_TYPE => $JSON->{QUERY}{$mz->{MASSES_MZ_QUERY}}{RECORDS}{$record}{precursor_type},
						ENTRY_INSTRUMENT_TYPE => $JSON->{QUERY}{$mz->{MASSES_MZ_QUERY}}{RECORDS}{$record}{instrument_type},
						ENTRY_SCORE => $JSON->{QUERY}{$mz->{MASSES_MZ_QUERY}}{RECORDS}{$record}{score},
						ENTRY_ENTRY_INCHI => $JSON->{QUERY}{$mz->{MASSES_MZ_QUERY}}{RECORDS}{$record}{inchi},			
		    		) ;
	    			push ( @{ $tbody_object->[$index_page]{MASSES}[$index_mz]{ENTRIES} }, \%entry) ;

    			$index_entry++ ;	
    		} ## end foreach record
    		if ($check_noentry == 0 ) {
    			my %entry = (
		    			ENTRY_COLOR => $tbody_object->[$index_page]{MASSES}[$index_mz]{MZ_COLOR},
		    			ENTRY_ENTRY_NAME  => 'UNKNOWN',
		   				ENTRY_ENTRY_ID => 'NONE',
		   				ENTRY_ENTRY_ID2 => '',
						ENTRY_FORMULA => 'n/a',
						ENTRY_CPD_MZ => 'n/a',
						ENTRY_MS_TYPE => 'n/a',
						ENTRY_PRECURSOR_TYPE => 'n/a',
						ENTRY_INSTRUMENT_TYPE => 'n/a',
						ENTRY_SCORE => 0,
						ENTRY_ENTRY_INCHI => 'n/a',
		    		) ;
		    		push ( @{ $tbody_object->[$index_page]{MASSES}[$index_mz]{ENTRIES} }, \%entry) ;
    		}
    		
    		## sorted by score
    		my @sorted = () ;
    		my @temp = @{ $tbody_object->[$index_page]{MASSES}[$index_mz]{ENTRIES} } ;
    		if (scalar (@temp) > 1 ) { ## for mz without record (only one entry with NA or 0 values)
    			@sorted = sort {  $b->{ENTRY_SCORE} <=> $a->{ENTRY_SCORE} } @temp;
    		}
    		else {
    			@sorted = @temp;
    		}
    		
    		$tbody_object->[$index_page]{MASSES}[$index_mz]{ENTRIES} = \@sorted ;
    		
    		$index_mz ++ ;
    		$index_mz_continous ++ ;

    	} ## End foreach mz
    	$index_page++ ;
    	
    } ## End foreach page
#    print Dumper $tbody_object ;
    return($tbody_object) ;
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