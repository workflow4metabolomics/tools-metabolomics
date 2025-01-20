package lib::writter ;

use strict;
use warnings ;
use Exporter ;
use Carp ;

use Data::Dumper ;
use JSON ;
use HTML::Template ;

use vars qw($VERSION @ISA @EXPORT %EXPORT_TAGS);

our $VERSION = "1.0";
our @ISA = qw(Exporter);
our @EXPORT = qw( write_csv_skel write_xls_skel write_json_skel write_html_skel );
our %EXPORT_TAGS = ( ALL => [qw( write_csv_skel write_xls_skel write_json_skel write_html_skel )] );

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

=head2 METHOD write_csv_skel

	## Description : prepare and write csv output file
	## Input : $csv_file, $rows
	## Output : $csv_file
	## Usage : my ( $csv_file ) = write_csv_skel( $csv_file, $rows ) ;
	
=cut
## START of SUB
sub write_csv_skel {
	## Retrieve Values
    my $self = shift ;
    my ( $csv_file, $rows ) = @_ ;
    
    my $ocsv = lib::csv::new( {is_binary => 1 , quote_binary => 0, quote_char => undef }) ;
	my $csv = $ocsv->get_csv_object("\t") ;
	$ocsv->write_csv_from_arrays($csv, $$csv_file, $rows) ;
    
    return($csv_file) ;
}
## END of SUB

=head2 METHOD write_xls_skel

	## Description : prepare and write xls output file
	## Input : $xls_file, $rows
	## Output : $xls_file
	## Usage : my ( $xls_file ) = write_xls_skel( $xls_file, $rows ) ;
	
=cut
## START of SUB
sub write_xls_skel {
	## Retrieve Values
    my $self = shift ;
    my ( $out_xls, $mzs, $pcs, $pcgroups, $records ) = @_ ;
    
    my $results = undef ;
    my $i = 0 ;
    
    open(XLS, '>:utf8', "$$out_xls") or die "Cant' create the file $$out_xls\n" ;
    print XLS "ID\tPCGROUP\tQuery(Da)\tScore\tMetabolite_name\tCpd_Mw(Da)\tFormula\tAdduct\tMASSBANK_ID\tInstrument\tMS_level\n" ;

    $results = ['ID','PCGROUP','Query(Da)','Score','Metabolite_name','Cpd_Mw(Da)','Formula','Adduct','MASSBANK_ID','Instrument','MS_level'] ;
    
    foreach my $pc (@{$pcs}) {
    	
    	if ($pcgroups->{$pc}) {
#    		print "------>$pc - $pcgroups->{$pc}{annotation}{num_res}\n" ;
    		
    		if ( $pcgroups->{$pc}{'annotation'} ) {
    			my $result = undef ;
    			my $well_id = "mz_0".sprintf("%04s", $i+1 ) ;
    			
    			if ($pcgroups->{$pc}{'annotation'}{'num_res'} > 0) {

    				my @entries = @{$pcgroups->{$pc}{'annotation'}{'res'} } ;
    				my $status = undef ;
    				foreach my $entry (@entries) {
    					my $match = undef ;
    					## manage if the queried mz is really in the mzs spectrum list...
    					
    					if ( $pcgroups->{$pc}{'enrich_annotation'}{$mzs->[$i]} ) {
    						
    						my @matching_ids = @{$pcgroups->{$pc}{'enrich_annotation'}{$mzs->[$i]}} ;
    						
    						## 
    						if ( scalar @matching_ids == 0 ) {
    							$result .= $well_id."\t".$pc."\t".$mzs->[$i]."\t".'0'."\t".'UNKNOWN'."\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\n" ;
	    						print XLS "$well_id\t$pc\t$mzs->[$i]\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n" ;
	    						last ;
    						}
    						else {
    							# search the massbank matched id
    							foreach (@matching_ids) {
	    							if ($_ eq $entry->{'id'} ) {
	    								$match = 'TRUE' ;
	    								last ;
	    							}
	    						}
	    						
	    						if ( ( defined $match ) and ($match eq 'TRUE') ) {
	    							## sort by ['ID','PCGROUP','Query(Da)','Score','Metabolite_name','Cpd_Mw(Da)','Formula','Adduct','MASSBANK_ID','Instrument','MS_level']
	    							
			    					## print mz_id
				    				if ($mzs->[$i]) { 	print XLS "$well_id\t" ; $result .= $well_id."\t" ; 	}
				    				else {				print XLS "NA\t" ; }
				    				## print submitted pcgroup
				    				if ($pc ) { 	 print XLS "$pc\t" ; $result .= $pc."\t" ; 	} ## pb de clean de la derniere ligne !!!!!!
				    				else {		print XLS "NA\t" ; }
				    				## print Query(Da)
				    				if ($mzs->[$i]) { 	print XLS "$mzs->[$i]\t" ; $result .= $mzs->[$i]."\t" ; 	}
				    				else {				print XLS "NA\t" ; }
				    				
				    				## print Score
				    				if ($entry->{'score'}) { 	print XLS "$entry->{'score'}\t" ; $result .= $entry->{'score'}."\t" ; 	}
				    				else {				print XLS "NA\n" ; }
				    				## print Met_name
				    				if ($entry->{'id'}) { 	print XLS "$records->{$entry->{'id'}}{names}[0]\t" ; $result .= $records->{$entry->{'id'}}{names}[0]."\t" ; 	}
				    				else {				print XLS "NA\t" ; }
				    				## print Cpd_mw
				    				if ($entry->{'exactMass'}) { 	print XLS "$entry->{'exactMass'}\t" ; $result .= $entry->{'exactMass'}."\t" ; 	}
				    				else {				print XLS "NA\t" ; }
				    				## print Formula
				    				if ($entry->{'formula'}) { 	print XLS "$entry->{'formula'}\t" ; $result .= $entry->{'formula'}."\t" ; 	}
				    				else {				print XLS "NA\t" ; }
				    				## print Adduct (precursor type)
				    				if ($entry->{'id'}) { 	print XLS "$records->{$entry->{'id'}}{precursor_type}\t" ; $result .= $records->{$entry->{'id'}}{precursor_type}."\t" ; 	}
				    				else {				print XLS "NA\t" ; }
				    				## print Massbank ID
				    				if ($entry->{'id'}) { 	print XLS "$entry->{'id'}\t" ; $result .= $entry->{'id'}."\t" ; 	}
				    				else {				print XLS "NA\t" ; }
				    				## print Instrument
				    				if ($entry->{'id'}) { 	print XLS "$records->{$entry->{'id'}}{instrument_type}\t" ; $result .= $records->{$entry->{'id'}}{instrument_type}."\t" ; 	}
				    				else {				print XLS "NA\t" ; }
				    				## print MS_Level
									if ($entry->{'id'}) { 	print XLS "$records->{$entry->{'id'}}{ms_type}\n" ; $result .= $records->{$entry->{'id'}}{ms_type}."\n" ; 	}
				    				else {				print XLS "NA\n" ; }
	
	    						}
	    						## else match is not TRUE
	    						else {
	    							next ;
	    						}
    						}
    					}
    				} ## End foreach entries
    			}
    			else {
    				$result .= $well_id."\t".$pc."\t".$mzs->[$i]."\t".'0'."\t".'UNKNOWN'."\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\n" ;
    				print XLS "$well_id\t$pc\t$mzs->[$i]\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n" ;
    			}
    		}
    		else{
    			warn "Not possible to get number of found ids on MassBank\n" ;
    		}
    	}
    	else {
    		croak "No such pc group exists in your pcgroups object - No xls written\n" ;
    	}
    	$i++ ;
    	
    } ## End foreach pcs
    
	close(XLS) ;
    return($results) ;
}
## END of SUB

=head2 METHOD write_json_skel

	## Description : prepare and write json output file
	## Input : $json_file, $scalar
	## Output : $json_file
	## Usage : my ( $json_file ) = write_json_skel( $csv_file, $scalar ) ;
	
=cut
## START of SUB
sub write_json_skel {
	## Retrieve Values
    my $self = shift ;
    my ( $json_file, $scalar ) = @_ ;
    
    my $utf8_encoded_json_text = encode_json $scalar ;
    open(JSON, '>:utf8', "$$json_file") or die "Cant' create the file $$json_file\n" ;
    print JSON $utf8_encoded_json_text ;
    close(JSON) ;
    
    return($json_file) ;
}
## END of SUB

=head2 METHOD write_html_skel

	## Description : prepare and write the html output file
	## Input : $html_file_name, $html_object, $html_template
	## Output : $html_file_name
	## Usage : my ( $html_file_name ) = write_html_skel( $html_file_name, $html_object ) ;
	
=cut
## START of SUB
sub write_html_skel {
	## Retrieve Values
    my $self = shift ;
    my ( $html_file_name,  $html_object, $pages , $search_condition, $html_template, $js_path, $css_path ) = @_ ;
    
    my $html_file = $$html_file_name ;
    
    if ( defined $html_file ) {
		open ( HTML, ">$html_file" ) or die "Can't create the output file $html_file " ;
		
		if (-e $html_template) {
			my $ohtml = HTML::Template->new(filename => $html_template);
			$ohtml->param(  JS_GALAXY_PATH => $js_path, CSS_GALAXY_PATH => $css_path  ) ;
			$ohtml->param(  CONDITIONS => $search_condition  ) ;
			$ohtml->param(  PAGES_NB => $pages  ) ;
			$ohtml->param(  PAGES => $html_object  ) ;
			print HTML $ohtml->output ;
		}
		else {
			croak "Can't fill any html output : No template available ($html_template)\n" ;
		}
		
		close (HTML) ;
    }
    else {
    	croak "No output file name available to write HTML file\n" ;
    }
    return(\$html_file) ;
}
## END of SUB

1 ;


__END__

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

 perldoc writter.pm

=head1 Exports

=over 4

=item :ALL is ...

=back

=head1 AUTHOR

Franck Giacomoni E<lt>franck.giacomoni@clermont.inra.frE<gt>

=head1 LICENSE

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=head1 VERSION

version 1 : 14 / 08 / 2015

version 2 : ??

=cut