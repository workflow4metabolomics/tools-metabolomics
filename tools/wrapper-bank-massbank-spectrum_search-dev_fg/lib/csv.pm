package lib::csv ;

use strict;
use warnings ;
use Exporter ;
use Carp ;

use Text::CSV ;

use Data::Dumper ;

use vars qw($VERSION @ISA @EXPORT %EXPORT_TAGS);

our $VERSION = "1.0";
our @ISA = qw(Exporter);
our @EXPORT = qw( get_csv_object get_value_from_csv get_value_from_csv_multi_header );
our %EXPORT_TAGS = ( ALL => [qw( get_csv_object get_value_from_csv get_value_from_csv_multi_header )] );

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

=head2 METHOD get_csv_object

	## Description : builds a csv object and etablishes format
	## Input : $separator
	## Output : $csv
	## Usage : my ( $csv ) = get_csv_object( $separator ) ;
	
=cut
## START of SUB
sub get_csv_object {
	## Retrieve Values
    my $self = shift ;
    my ( $separator ) = @_ ;
    
#    my $csv = Text::CSV->new({'sep_char' => "$separator"});
    my $csv = Text::CSV->new ( {'sep_char' => "$separator", binary => 1,  } )  # should set binary attribute.
    	or die "Cannot use CSV: ".Text::CSV->error_diag ();

    return($csv) ;
}
## END of SUB

=head2 METHOD get_value_from_csv

	## Description : extract a targeted column in a csv file 
	## Input : $csv, $file, $column, $is_header
	## Output : $value
	## Usage : my ( $value ) = get_value_from_csv( $csv, $file, $column, $is_header ) ;
	
=cut
## START of SUB
sub get_value_from_csv {
	## Retrieve Values
    my $self = shift ;
    my ( $csv, $file, $column, $is_header ) = @_ ;
    
    my @value = () ;
    
    ## Adapte the number of the colunm : (nb of column to position in array)
	$column = $column - 1 ;
    
    open (CSV, "<", $file) or die $! ;
	
	my $line = 0 ;
	
	while (<CSV>) {
		$line++ ;
	    chomp $_ ;
		# file has a header
		if ( defined $is_header ) { if ($line == 1) { next ; } }
		# parsing the targeted column
	    if ( $csv->parse($_) ) {
	        my @columns = $csv->fields();
	        push ( @value, $columns[$column] ) ;
	    }
	    else {
	        my $err = $csv->error_input;
	        die "Failed to parse line: $err";
	    }
	}
	close CSV;
    return(\@value) ;
}
## END of SUB

=head2 METHOD get_value_from_csv_multi_header

	## Description : extract a targeted column in a csv file 
	## Input : $csv, $file, $column, $is_header, $nb_header
	## Output : $value
	## Usage : my ( $value ) = get_value_from_csv_multi_header( $csv, $file, $column, $is_header, $nb_header ) ;
	
=cut
## START of SUB
sub get_value_from_csv_multi_header {
	## Retrieve Values
    my $self = shift ;
    my ( $csv, $file, $column, $is_header, $nb_header ) = @_ ;
    
    my @value = () ;
    
    ## Adapte the number of the colunm : (nb of column to position in array)
	$column = $column - 1 ;
    
    open (CSV, "<", $file) or die $! ;
	
	my $line = 0 ;
	
	while (<CSV>) {
		$line++ ;
	    chomp $_ ;
		# file has a header
		if ( defined $is_header and $is_header eq 'yes') { if ($line <= $nb_header) { next ; } }
		# parsing the targeted column
	    if ( $csv->parse($_) ) {
	        my @columns = $csv->fields();
	        my $value = $columns[$column] ;
	        $value =~s/\r|\n// ;
	        push ( @value,  $value ) ;
	    }
	    else {
	        my $err = $csv->error_input;
	        die "Failed to parse line: $err";
	    }
	}
	close CSV;
    return(\@value) ;
}
## END of SUB

=head2 METHOD parse_csv_object

	## Description : parse_all csv object and return a array of rows
	## Input : $csv, $file
	## Output : $csv_matrix
	## Usage : my ( $csv_matrix ) = parse_csv_object( $csv, $file ) ;
	
=cut
## START of SUB
sub parse_csv_object {
	## Retrieve Values
    my $self = shift ;
    my ( $csv, $file ) = @_ ;
    
    my @csv_matrix = () ;
    
	open my $fh, "<:encoding(utf8)", $$file or die "Can't open csv file $$file: $!";
	
	while ( my $row = $csv->getline( $fh ) ) {
	    push @csv_matrix, $row;
	}
	$csv->eof or $csv->error_diag();
	close $fh;
    
    return(\@csv_matrix) ;
}
## END of SUB

=head2 METHOD parse_allcsv_object

	## Description : parse_all csv object and return a array of rows with or without header
	## Input : $csv, $file, $keep_header
	## Output : $csv_matrix
	## Usage : my ( $csv_matrix ) = parse_csv_object( $csv, $file, $keep_header ) ;
	
=cut
## START of SUB
sub parse_allcsv_object {
	## Retrieve Values
    my $self = shift ;
    my ( $csv, $file, $keep_header ) = @_ ;
    
    my @csv_matrix = () ;
    my $line = 1 ;
    
	open my $fh, "<:encoding(utf8)", $$file or die "Can't open csv file $$file: $!";
	
	while ( my $row = $csv->getline( $fh ) ) {
		if ( ( $keep_header eq 'n' )  and ($line == 1) ) {  }
		else { push @csv_matrix, $row; }
	    $line ++ ;
	}
	my $status = $csv->eof or $csv->error_diag();
	close $fh;
    
    return(\@csv_matrix, $status) ;
}
## END of SUB


=head2 METHOD write_csv_from_arrays

	## Description : write a csv file from list of rows
	## Input : $csv, $file_name, $rows
	## Output : $csv_file
	## Usage : my ( $csv_file ) = write_csv_from_arrays( $csv, $file_name, $rows ) ;
	
=cut
## START of SUB
sub write_csv_from_arrays {
	## Retrieve Values
    my $self = shift ;
    my ( $csv, $file_name, $rows ) = @_ ;
    
    my $fh = undef ;
    $csv->eol ("\n"); ##  end-of-line string to add to rows
    $csv->quote_char(undef) ;
    open $fh, ">:encoding(utf8)", "$file_name" or die "$file_name: $!";
    
	my $status = $csv->print ($fh, $_) for @{$rows};
	close $fh or die "$file_name: $!";
    
    return(\$file_name) ;
}
## END of SUB

1 ;


__END__

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

 perldoc csv.pm

=head1 Exports

=over 4

=item :ALL is get_csv_object, get_value_from_csv

=back

=head1 AUTHOR

Franck Giacomoni E<lt>franck.giacomoni@clermont.inra.frE<gt>

=head1 LICENSE

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=head1 VERSION

version 1 : 23 / 10 / 2013

version 2 : ??

=cut