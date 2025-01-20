package lib::conf ;

use strict;
use warnings ;
use Exporter ;
use Carp ;
use Data::Dumper ;

use vars qw($VERSION @ISA @EXPORT %EXPORT_TAGS);

our $VERSION = "1.0" ;
our @ISA = qw(Exporter) ;
our @EXPORT = qw( as_conf get_value_from_conf check_path_and_file ) ;
our %EXPORT_TAGS = ( ALL => [qw( as_conf get_value_from_conf )] ) ;

=head1 NAME

conf - A module for manage pfem conf file

=head1 SYNOPSIS

    use conf ;
    my $object = conf->new() ;

=head1 DESCRIPTION

This module does manage conf file (extract all or selected fields)

=head1 METHODS

Methods are :

=head2 METHOD new

	## Description : new
	## Input : $self
	## Ouput : bless $self ;
	## Usage : new() ;

=cut
## START of SUB
sub new {
    ## Variables
    my $self={};
    bless($self) ;
    return $self ;
}
### END of SUB

=head2 METHOD as_conf

	## Description : permet de creer l'object conf a partir d'un fichier de conf de type KEY=VALUE
	## Input : $file
	## Ouput : $oConf (a hash)
	## Usage : my ( $oConf ) = as_conf( $file ) ;
	
=cut
## START of SUB
sub as_conf {
	## Retrieve Values
    my $self = shift ;
    my ( $file, $separator ) = @_ ;
    
#    if (!defined $separator) { $separator = ';' } ## set separator to ;
    
    if ( !defined $file )  {  croak "Can't create object with an none defined file\n" ; }
    
    my %Conf = () ; ## Hash devant contenir l'ensemble des parametres locaux
	
	if (-e $file) {
		open (CFG, "<$file") or die "Can't open $file\n" ;
		while (<CFG>) {
			chomp $_ ;
			if ( $_ =~ /^#(.*)/)  {	next ; }
			elsif ($_ =~/^(\w+?)=(.*)/) { ## ALPHANUMERIC OR UNDERSCORE ONLY FOR THE KEY AND ANYTHING ELSE FOR VALUE
				
				my ($key, $value) = ($1, $2) ;
				
				if (defined $separator) {
					if ( $value=~/$separator/ ) { ## is a list to split
						my @tmp = split(/$separator/ , $value) ;
						$Conf{$key} = \@tmp ;
					}
				}
				else {
					$Conf{$key} = $value ;
				}
			}
		}
		close(CFG) ;
	}
	else { 
		croak "Can't create object with an none existing file\n" ;
	}
	
    return ( \%Conf ) ;
}
## END of SUB

=head2 METHOD as_conf_list

	## Description : permet de charger une liste txt en array
	## Input : $file
	## Output : elements
	## Usage : my ( elements ) = as_conf_list( $conf_file ) ;
	
=cut
## START of SUB
sub as_conf_list {
	## Retrieve Values
    my $self = shift ;
    my ( $file ) = @_ ;
    
    my @elements = () ;
    if ( !defined $file )  {  croak "Can't create object with an none defined file\n" ; }
    
    if (-e $file) {
		open (CFG, "<$file") or die "Can't open $file\n" ;
		while (<CFG>) {
			chomp $_ ;
			if ( $_ =~ /^#(.*)/)  {	next ; }
			elsif ($_ =~/^(.*)/) { if (defined $1) { push (@elements, $1) ; } 	}
		}
    }
    else {
		croak "Can't create object with an none existing file\n" ;
	}
    return(\@elements) ;
}
## END of SUB

=head2 METHOD get_value_from_conf

	## Description : permet de retourner une valeur du hash de conf a partir d'une key
	## Input : $oConf, $Key
	## Ouput : $Value
	## Usage : my ( $Value ) = get_value_from_conf( $oConf, $Key ) ;
	
=cut
## START of SUB
sub get_value_from_conf {
	## Retrieve Values
    my $self = shift ;
    my ( $oConf, $Key ) = @_ ;
    
    my $Value = undef ;
    
    if ( defined $oConf ) {
    	if ( defined $oConf->{$Key} ) {
    		$Value = $oConf->{$Key} ;
    	}
    }
    else {
    	croak "Can't manage value with undefined object\n" ;
    }
    return($Value) ;
}
## END of SUB

=head2 METHOD get_value_from_conf

	## Description : permet de retourner une valeur du hash de conf a partir d'une key
	## Input : $oConf, $Key
	## Ouput : $Value
	## Usage : my ( $Value ) = get_value_from_conf( $oConf, $Key ) ;
	
=cut
## START of SUB
sub split_value_from_conf {
	## Retrieve Values
    my $self = shift ;
    my ( $oConf, $Key, $sep ) = @_ ;
    
    my $value = undef ;
    my @values = () ;
    
    if ( defined $oConf ) {
    	if ( defined $oConf->{$Key} ) {
    		$value = $oConf->{$Key} ;
    		@values = split ( /$sep/, $value) ;
    	}
    }
    else {
    	croak "Can't manage value with undefined object\n" ;
    }
    return(\@values) ;
}
## END of SUB


=head2 METHOD check_path_and_file

	## Description : permet de verifier les path et la presence des exe decrits dans le file conf. Bloque le script en cas de probleme
	## Input : $oConfs
	## Ouput : NA
	## Usage : &get_value_from_conf( $oConf ) ;
	
=cut
## START of SUB
sub check_path_and_file {
	
	my $self = shift ;
	my ( $oConfs ) = @_ ;
	
	foreach my $conf ( keys %{ $oConfs } ) {
		if ( $conf =~ /^FILE/ ) {
			if ( -e $oConfs->{$conf} ) {
				if ( -s $oConfs->{$conf} ) { next ; }
				else { carp "[Warning] : The size of file $oConfs->{$conf} is null\n" ; }
			}
			else {
					carp "[Warning] : The file $oConfs->{$conf} doesn't exist\n" ;
			}
		}
		elsif ( $conf =~ /^PATH/ ) {
			if ( -d $oConfs->{$conf} ) { next ; }
			else { carp "[Warning] :  The dir $oConfs->{$conf} doesn't exist\n" ;	}
		}
		else { 	next ; 	}
	}
	return ;
}
## END of SUB

1 ;


__END__

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

 perldoc conf.pm


=head1 Exports

=over 4

=item :ALL is as_conf get_value_from_conf

=back

=head1 AUTHOR

Franck Giacomoni E<lt>franck.giacomoni@clermont.inra.frE<gt>

=head1 LICENSE

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=head1 VERSION

version 1 : 10 / 02 / 2013

version 2 : ??

=cut