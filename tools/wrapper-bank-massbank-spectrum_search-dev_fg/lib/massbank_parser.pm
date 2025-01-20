package lib::massbank_parser ;

use strict;
use warnings ;
use Exporter ;
use Carp ;

use File::Basename;

use Data::Dumper ;

use vars qw($VERSION @ISA @EXPORT %EXPORT_TAGS);

our $VERSION = "1.0" ;
our @ISA = qw(Exporter) ;
our @EXPORT = qw( getChemNamesFromString getPeaksFromString ) ;
our %EXPORT_TAGS = ( ALL => [qw( getChemNamesFromString getPeaksFromString )] ) ;

=head1 NAME

parser::chem::massbank - An example module

=head1 SYNOPSIS

    use parser::chem::massbank ;
    my $object = parser::chem::massbank->new();
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

=head2 METHOD get_list_of_analysis_intrument_names

	## Description : permt de retourner la liste des nom uniques des instruments utilises
	## Input : $dir, $ms_files (a list of files)
	## Output : $names
	## Usage : my ( $names ) = get_list_of_analysis_intrument_names( $ms_files ) ;
	
=cut
## START of SUB
sub get_list_of_analysis_intrument_names {
	## Retrieve Values
    my $self = shift ;
    my ( $dir, $ms_files ) = @_ ;
    my (%tmp_names, @names) = ( (), () ) ;
    foreach my $ms_file (@{$ms_files}) {
    	my $file = $dir.'\\'.$ms_file ;
    	if ( ( defined $file ) and ( -e $file )) {
	    	open(MS, "<$file") or die "Cant' read the file $file\n" ;
	    	while ( my $field = <MS> ){
	    		chomp $field ;
	        	if ($field =~/AC\$INSTRUMENT:(.*)/) {
	        		if ( $tmp_names{$1} ) { last ; }
	        		else { $tmp_names{$1} = 1 ; push (@names, $1) ; }
	        	}
	    	}
	    	close(MS) ;
	    }
	    else { 
	    	croak "Can't work with a undef / none existing massbank file\n" ;
	    }
    }
    return(\@names) ;
}
## END of SUB

=head2 METHOD get_analysis_instruments_data

	## Description : permet de recuperer tous les champs d'un object massbank
	## Input : $ms_file
	## Output : $features
	## Usage : my ( $features ) = get_analysis_instruments_data( $ms_file ) ;
	
=cut
## START of SUB
sub get_analysis_instruments_data {
	## Retrieve Values
    my $self = shift ;
    my ( $ms_file ) = @_ ;
    
    my $control = 0 ;
    my %features = (
    	'name' => undef,
    	'type' => undef,
    ) ;
    if ( ( defined $ms_file ) and ( -e $ms_file )) {
    	open(MS, "<$ms_file") or die "Cant' read the file $ms_file\n" ;
    	while ( my $field = <MS> ){
    		chomp $field ;    		
    		if ($field =~/AC\$INSTRUMENT: (.*)/) { $features{'name'} = $1 ; $control++; }
	    	elsif ($field =~/AC\$INSTRUMENT_TYPE: (.*)/) { $features{'type'} = $1 ; $control++; }
	    	else { next ; }
    	}
    	close(MS) ;
    }
    else {
    	croak "Can't work with a undef / none existing massbank file\n" ;
    }
    if ($control == 0) { %features = () ;  }
    return(\%features) ;
}
## END of SUB
     
=head2 METHOD get_ms_methods_data

	## Description : permet de recuperer tous les champs d'un object massbank
	## Input : $ms_file
	## Output : $features
	## Usage : my ( $features ) = get_ms_methods_data( $ms_file ) ;
	
=cut
## START of SUB
sub get_ms_methods_data {
	## Retrieve Values
    my $self = shift ;
    my ( $ms_file ) = @_ ;
    
    my $control = 0 ;
    my %features = (
    	'ion_mode' => undef,
    	'ms_type' => undef,
    	'collision_energy' => undef,
    	'collision_gas' => undef,
    	'desolvation_gas_flow' => undef,
    	'desolvation_temperature' => undef,
    	'ionization_energy' => undef,
    	'laser' => undef,
    	'matrix' => undef,
    	'mass_accuracy' => undef,
    	'reagent_gas' => undef,
    	'scanning' => undef
    ) ;
    if ( ( defined $ms_file ) and ( -e $ms_file )) {
    	open(MS, "<$ms_file") or die "Cant' read the file $ms_file\n" ;
    	while ( my $field = <MS> ){
    		chomp $field ;    		
    		if ($field =~/AC\$MASS_SPECTROMETRY: ION_MODE:(.*)/) { $features{'ion_mode'} = $1 ; $control++; } # mandatory
	    	elsif ($field =~/AC\$MASS_SPECTROMETRY: MS_TYPE:(.*)/) { $features{'ms_type'} = $1 ; $control++; } # mandatory
	    	elsif ($field =~/AC\$MASS_SPECTROMETRY: COLLISION_ENERGY(.*)/) { $features{'collision_energy'} = $1 ; $control++; } # optionnal
	    	elsif ($field =~/AC\$MASS_SPECTROMETRY: COLLISION_GAS(.*)/) { $features{'collision_gas'} = $1 ; $control++; } # optionnal
	    	elsif ($field =~/AC\$MASS_SPECTROMETRY: DESOLVATION_GAS_FLOW(.*)/) { $features{'desolvation_gas_flow'} = $1 ; $control++;  } # optionnal
	    	elsif ($field =~/AC\$MASS_SPECTROMETRY: DESOLVATION_TEMPERATURE(.*)/) { $features{'desolvation_temperature'} = $1 ; $control++; } # optionnal
	    	elsif ($field =~/AC\$MASS_SPECTROMETRY: IONIZATION_ENERGY(.*)/) { $features{'ionization_energy'} = $1 ; $control++;  } # optionnal
	    	elsif ($field =~/AC\$MASS_SPECTROMETRY: LASER(.*)/) { $features{'laser'} = $1 ; $control++; } # optionnal
	    	elsif ($field =~/AC\$MASS_SPECTROMETRY: MATRIX(.*)/) { $features{'matrix'} = $1 ; $control++; } # optionnal
	    	elsif ($field =~/AC\$MASS_SPECTROMETRY: MASS_ACCURACY(.*)/) { $features{'mass_accuracy'} = $1 ; $control++; } # optionnal
	    	elsif ($field =~/AC\$MASS_SPECTROMETRY: REAGENT_GAS(.*)/) { $features{'reagent_gas'} = $1 ; $control++; } # optionnal
	    	elsif ($field =~/AC\$MASS_SPECTROMETRY: SCANNING(.*)/) { $features{'scanning'} = $1 ; $control++; } # optionnal
	    	else { next ; }
    	}
    	close(MS) ;
    }
    else {
    	croak "Can't work with a undef / none existing massbank file\n" ;
    }
    ## vide l'object si undef
    if ($control == 0) { %features = () ;  }
    return(\%features) ;
}
## END of SUB

=head2 METHOD get_solvents_data

	## Description : permet de recuperer tous les champs d'un object massbank
	## Input : $ms_file
	## Output : $features
	## Usage : my ( $features ) = get_solvents_data( $ms_file ) ;
	
=cut
## START of SUB
sub get_solvents_data {
	## Retrieve Values
    my $self = shift ;
    my ( $ms_file ) = @_ ;
    
    my @features = () ;
    if ( ( defined $ms_file ) and ( -e $ms_file )) {
    	open(MS, "<$ms_file") or die "Cant' read the file $ms_file\n" ;
    	while ( my $field = <MS> ){
    		chomp $field ;    		
    		if ($field =~/AC\$CHROMATOGRAPHY: SOLVENT(.*)/) { push(@features, 'Solvent '.$1 ) ;  }
	    	else { next ; }
    	}
    	close(MS) ;
    }
    else {
    	croak "Can't work with a undef / none existing massbank file\n" ;
    }
    return(\@features) ;
}
## END of SUB

=head2 METHOD get_sample_data

	## Description : permet de recuperer tous les champs d'un object massbank
	## Input : $ms_file
	## Output : $features
	## Usage : my ( $features ) = get_sample_data( $ms_file ) ;
	
=cut
## START of SUB
sub get_sample_data {
	## Retrieve Values
    my $self = shift ;
    my ( $ms_file ) = @_ ;
    
    my $control = 0;
    my %features = (
    	'sample_type' => undef,
    ) ;
    if ( ( defined $ms_file ) and ( -e $ms_file )) {
    	open(MS, "<$ms_file") or die "Cant' read the file $ms_file\n" ;
    	while ( my $field = <MS> ){
    		chomp $field ;    		
    		if ($field =~/SP\$SAMPLE(.*)/) { $features{'sample_type'} = $1 ; $control++ ; }
	    	else { next ; }
    	}
    	close(MS) ;
    }
    else {
    	croak "Can't work with a undef / none existing massbank file\n" ;
    }
    if ($control == 0) { %features = () ;  }
    return(\%features) ;
}
## END of SUB

=head2 METHOD get_chromato_methods_data

	## Description : permet de recuperer tous les champs d'un object massbank
	## Input : $ms_file
	## Output : $features
	## Usage : my ( $features ) = get_chromato_methods_data( $ms_file ) ;
	
=cut
## START of SUB
sub get_chromato_methods_data {
	## Retrieve Values
    my $self = shift ;
    my ( $ms_file ) = @_ ;
    
    my $control = 0 ;
    my %features = (
    	'capillary_voltage' => undef,
    	'column_name' => undef,
    	'column_temperature' => undef,
    	'flow_gradient' => undef,
    	'flow_rate' => undef,
    	'retention_time' => undef,
    ) ;
    if ( ( defined $ms_file ) and ( -e $ms_file )) {
    	open(MS, "<$ms_file") or die "Cant' read the file $ms_file\n" ;
    	while ( my $field = <MS> ){
    		chomp $field ;    		
    		if ($field =~/AC\$CHROMATOGRAPHY: CAPILLARY_VOLTAGE (.*)/) { $features{'capillary_voltage'} = $1 ; $control++ ; }
	    	elsif ($field =~/AC\$CHROMATOGRAPHY: COLUMN_NAME (.*)/) { $features{'column_name'} = $1 ; $control++ ; }
	    	elsif ($field =~/AC\$CHROMATOGRAPHY: COLUMN_TEMPERATURE (.*)/) { $features{'column_temperature'} = $1 ; $control++ ; }
	    	elsif ($field =~/AC\$CHROMATOGRAPHY: FLOW_GRADIENT (.*)/) { $features{'flow_gradient'} = $1 ; $control++ ; }
	    	elsif ($field =~/AC\$CHROMATOGRAPHY: FLOW_RATE (.*)/) { $features{'flow_rate'} = $1 ; $control++ ; }
	    	elsif ($field =~/AC\$CHROMATOGRAPHY: RETENTION_TIME (.*)/) { $features{'retention_time'} = $1 ; $control++ ; }
	    	else { next ; }
    	}
    	close(MS) ;
    	# for db field
    }
    else {
    	croak "Can't work with a undef / none existing massbank file\n" ;
    }
    if ($control == 0) { %features = () ;  }
    return(\%features) ;
}
## END of SUB

=head2 METHOD get_analytical_conditions_data

	## Description : permet de recuperer tous les champs d'un object massbank .. for massbank version < 2.0
	## Input : $ms_file
	## Output : $features
	## Usage : my ( $features ) = get_analytical_conditions_data( $ms_file ) ;
	
=cut
## START of SUB
sub get_analytical_conditions_data {
	## Retrieve Values
    my $self = shift ;
    my ( $ms_file ) = @_ ;
    my $control_ms = 0 ;
    my %features_ms = (
    	'ion_mode' => undef,
    	'ms_type' => undef,
    	'collision_energy' => undef,
    	'collision_gas' => undef,
    	'desolvation_gas_flow' => undef,
    	'desolvation_temperature' => undef,
    	'ionization_energy' => undef,
    	'laser' => undef,
    	'matrix' => undef,
    	'mass_accuracy' => undef,
    	'reagent_gas' => undef,
    	'scanning' => undef    	
    ) ;
    my $control_chrom = 0 ;
    my %features_chrom = (
    	'capillary_voltage' => undef,
    	'column_name' => undef,
    	'column_temperature' => undef,
    	'flow_gradient' => undef,
    	'flow_rate' => undef,
    	'retention_time' => undef   	
    ) ;
    if ( ( defined $ms_file ) and ( -e $ms_file )) {
    	open(MS, "<$ms_file") or die "Cant' read the file $ms_file\n" ;
    	while ( my $field = <MS> ){
    		chomp $field ;
    		## new = chromato_method	
    		if ($field =~/AC\$ANALYTICAL_CONDITION: CAPILLARY_VOLTAGE (.*)/) { $features_chrom{'capillary_voltage'} = $1 ; $control_chrom++ ; }
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: COLUMN_NAME (.*)/) { $features_chrom{'column_name'} = $1 ; $control_chrom++ ; }
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: COLUMN_TEMPERATURE( .*)/) { $features_chrom{'column_temperature'} = $1 ; $control_chrom++ ; }
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: FLOW_GRADIENT (.*)/) { $features_chrom{'flow_gradient'} = $1 ; $control_chrom++ ;  }
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: FLOW_RATE (.*)/) { $features_chrom{'flow_rate'} = $1 ; $control_chrom++ ; }
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: RETENTION_TIME (.*)/) { $features_chrom{'retention_time'} = $1 ; $control_chrom++ ; }
	    	## new = ms_method
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: ION_MODE (.*)/) { $features_ms{'ion_mode'} = $1 ; $control_ms++ ; } # mandatory
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: MS_TYPE (.*)/) { $features_ms{'ms_type'} = $1 ; $control_ms++ ; } # mandatory
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: COLLISION_ENERGY (.*)/) { $features_ms{'collision_energy'} = $1 ; $control_ms++ ; } # optionnal
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: COLLISION_GAS (.*)/) { $features_ms{'collision_gas'} = $1 ; $control_ms++ ; } # optionnal
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: DESOLVATION_GAS_FLOW (.*)/) { $features_ms{'desolvation_gas_flow'} = $1 ; $control_ms++ ; } # optionnal
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: DESOLVATION_TEMPERATURE (.*)/) { $features_ms{'desolvation_temperature'} = $1 ; $control_ms++ ; } # optionnal
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: IONIZATION_ENERGY (.*)/) { $features_ms{'ionization_energy'} = $1 ; $control_ms++ ; } # optionnal
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: LASER (.*)/) { $features_ms{'laser'} = $1 ; $control_ms++ ; } # optionnal
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: MATRIX (.*)/) { $features_ms{'matrix'} = $1 ; $control_ms++ ; } # optionnal
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: MASS_ACCURACY (.*)/) { $features_ms{'mass_accuracy'} = $1 ; $control_ms++ ; } # optionnal
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: REAGENT_GAS (.*)/) { $features_ms{'reagent_gas'} = $1 ; $control_ms++ ; } # optionnal
	    	elsif ($field =~/AC\$ANALYTICAL_CONDITION: SCANNING (.*)/) { $features_ms{'scanning'} = $1 ; $control_ms++ ; } # optionnal
	    	else { next ; }
    	}
    	close(MS) ;
    	# for db field
    }
    else {
    	croak "Can't work with a undef / none existing massbank file\n" ;
    }
    if ($control_ms == 0) { %features_ms = () ;  }
    if ($control_chrom == 0) { %features_chrom = () ;  }
    return(\%features_chrom, \%features_ms) ;
}
## END of SUB

=head2 METHOD get_spectrums_data

	## Description : permet de recuperer tous les champs d'un object massbank
	## Input : $ms_file
	## Output : $features
	## Usage : my ( $features ) = get_spectrums_data( $ms_file ) ;
	
=cut
## START of SUB
sub get_spectrums_data {
	## Retrieve Values
    my $self = shift ;
    my ( $ms_file ) = @_ ;
    my $control = 0 ;
    my %features = (
    	'ion_type' => undef,
    	'precursor_mz' => undef,
    	'precursor_type' => undef,
    	'num_peaks' => undef,
    ) ;
    if ( ( defined $ms_file ) and ( -e $ms_file )) {
    	open(MS, "<$ms_file") or die "Cant' read the file $ms_file\n" ;
    	while ( my $field = <MS> ){
    		chomp $field ;    		
    		if ($field =~/MS\$FOCUSED_ION: ION_TYPE(.*)/) { $features{'ion_type'} = $1 ; $control++ ; }
	    	elsif ($field =~/MS\$FOCUSED_ION: PRECURSOR_M\/Z(.*)/) { $features{'precursor_mz'} = $1 ; $control++ ; }
	    	elsif ($field =~/MS\$FOCUSED_ION: PRECURSOR_TYPE(.*)/) { $features{'precursor_type'} = $1 ; $control++ ; }
	    	elsif ($field =~/PK\$NUM_PEAK: (.*)/) { $features{'num_peaks'} = $1 ; $control++ ; }
	    	else { next ; }
    	}
    	close(MS) ;
    	# for db field
    }
    else {
    	croak "Can't work with a undef / none existing massbank file\n" ;
    }
    if ($control == 0) { %features = () ;  }
    return(\%features) ;
}
## END of SUB

=head2 METHOD get_peaks_data

	## Description : permet de recuperer tous les champs d'un object massbank
	## Input : $ms_file
	## Output : $features
	## Usage : my ( $features ) = get_peaks_data( $ms_file ) ;
	
=cut
## START of SUB
sub get_peaks_data {
	## Retrieve Values
    my $self = shift ;
    my ( $ms_file ) = @_ ;
    
    my @features = () ;
    my $peaks = 0 ;
    if ( ( defined $ms_file ) and ( -e $ms_file )) {
    	open(MS, "<$ms_file") or die "Cant' read the file $ms_file\n" ;
    	while ( my $field = <MS> ){
    		chomp $field ;
    		if ($field =~/PK\$PEAK: m\/z int\. rel\.int\./) { $peaks = 1 ; }
    		elsif ( $peaks == 1 ) { ## detected peak area
    			if ($field =~/\s+(\d+)\s+(\d+)\s+(\d+)/) {
    				my %tmp = ( 'mz' => $1, 'intensity' => $2, 'relative_intensity' => $3 ) ;
    				push (@features, \%tmp) ;
    			}
    			## for int = xx.xxx and mz = xxx.xxx
    			elsif ($field =~/\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+)/) {
    				my %tmp = ( 'mz' => $1, 'intensity' => $2, 'relative_intensity' => $3 ) ;
    				push (@features, \%tmp) ;
    			}
    			## for int = xx and mz = xxx.xxx
    			elsif ($field =~/\s+(\d+\.\d+)\s+(\d+)\s+(\d+)/) {
    				my %tmp = ( 'mz' => $1, 'intensity' => $2, 'relative_intensity' => $3 ) ;
    				push (@features, \%tmp) ;
    			}
    			## for int = xxxxx.xxx and mz = xxx
    			elsif ($field =~/\s+(\d+)\s+(\d+\.\d+)\s+(\d+)/) {
    				my %tmp = ( 'mz' => $1, 'intensity' => $2, 'relative_intensity' => $3 ) ;
    				push (@features, \%tmp) ;
    			}
    		}
	    	else { next ; }
    	}
    	close(MS) ;
    	# for db field
    }
    else {
    	croak "Can't work with a undef / none existing massbank file\n" ;
    }
    return(\@features) ;
}
## END of SUB

=head2 METHOD getPeaksFromString

	## Description : permet de recuperer la data peaks d'un record handler massbank
	## Input : $record
	## Output : $features
	## Usage : my ( $features ) = getPeaksFromString( $record ) ;
	
=cut
## START of SUB
sub getPeaksFromString {
	## Retrieve Values
    my $self = shift ;
    my ( $record ) = @_ ;
    
    my @features = () ;
    my $peaks = 0 ;
    if ( defined $record ) {
    	my @tmp = split(/\n/, $record) ;
    	foreach my $field (@tmp) {
    		if ($field =~/PK\$PEAK: m\/z int\. rel\.int\./) { $peaks = 1 ; }
    		elsif ( $peaks == 1 ) { ## detected peak area
    			if ($field =~/\s+(\d+)\s+(\d+)\s+(\d+)/) {
    				my %tmp = ( 'mz' => $1, 'intensity' => $2, 'relative_intensity' => $3 ) ;
    				push (@features, \%tmp) ;
    			}
    			## for int = xx.xxx and mz = xxx.xxx
    			elsif ($field =~/\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+)/) {
    				my %tmp = ( 'mz' => $1, 'intensity' => $2, 'relative_intensity' => $3 ) ;
    				push (@features, \%tmp) ;
    			}
    			## for int = xx and mz = xxx.xxx
    			elsif ($field =~/\s+(\d+\.\d+)\s+(\d+)\s+(\d+)/) {
    				my %tmp = ( 'mz' => $1, 'intensity' => $2, 'relative_intensity' => $3 ) ;
    				push (@features, \%tmp) ;
    			}
    			## for int = xxxxx.xxx and mz = xxx
    			elsif ($field =~/\s+(\d+)\s+(\d+\.\d+)\s+(\d+)/) {
    				my %tmp = ( 'mz' => $1, 'intensity' => $2, 'relative_intensity' => $3 ) ;
    				push (@features, \%tmp) ;
    			}
    			## for int = x.xxxex and m/z = xxx.xxx (int with exposant)
    			elsif ($field =~/\s+(\d+\.\d+)\s+(\d+\.\d+)e(\d)\s+(\d+)/) {
    				my %tmp = ( 'mz' => $1, 'intensity' => ($2*(10*$3)), 'relative_intensity' => $4 ) ;
    				push (@features, \%tmp) ;
    			}
    		}
	    	else { next ; }
    	}
    	# for db field
    }
    else {
    	croak "Can't work with a undef / none existing massbank handler\n" ;
    }
    return(\@features) ;
}
## END of SUB

=head2 METHOD getIdFromString

	## Description : get the accesion id of massbank record
	## Input : $record
	## Output : $id
	## Usage : my ( $id ) = getIdFromString ( $record ) ;
	
=cut
## START of SUB
sub getIdFromString {
    ## Retrieve Values
    my $self = shift ;
    my ( $record ) = @_;
    my ( $id ) = ( undef ) ;
    
    if ( defined $record ) {
    	my @tmp = split(/\n/, $record) ;
    	foreach my $field (@tmp) {
    		if ($field =~/ACCESSION:\s+(.+)/) { 
				$id = $1;
    		}
    	}
    	# for db field
    }
    else {
    	croak "Can't work with a undef / none existing massbank handler\n" ;
    }
    
    return ($id) ;
}
### END of SUB



=head2 METHOD getInstrumentTypeFromString

	## Description : get the instrument type of massbank record
	## Input : $record
	## Output : $instrumentType
	## Usage : my ( $instrumentType ) = getInstrumentTypeFromString ( $record ) ;
	
=cut
## START of SUB
sub getInstrumentTypeFromString {
    ## Retrieve Values
    my $self = shift ;
    my ( $record ) = @_;
    my ( $instrumentType ) = ( undef ) ;
    
    if ( defined $record ) {
    	my @tmp = split(/\n/, $record) ;
    	foreach my $field (@tmp) {
    		if ($field =~/INSTRUMENT_TYPE:\s+(.+)/) { 
				$instrumentType = $1;
    		}
    	}
    	# for db field
    }
    else {
    	croak "Can't work with a undef / none existing massbank handler\n" ;
    }
    
    return ($instrumentType) ;
}
### END of SUB

=head2 METHOD getFormulaFromString

	## Description : get the elementar formula of massbank record
	## Input : $record
	## Output : $formula
	## Usage : my ( $formula ) = getFormulaFromString ( $record ) ;
	
=cut
## START of SUB
sub getFormulaFromString {
    ## Retrieve Values
    my $self = shift ;
    my ( $record ) = @_;
    my ( $formula ) = ( undef ) ;
    
    if ( defined $record ) {
    	my @tmp = split(/\n/, $record) ;
    	foreach my $field (@tmp) {
    		if ($field =~/CH\$FORMULA:\s+(.+)/) { 
				$formula = $1;
    		}
    	}
    	# for db field
    }
    else {
    	croak "Can't work with a undef / none existing massbank handler\n" ;
    }
    
    return ($formula) ;
}
### END of SUB

=head2 METHOD getInchiFromString

	## Description : get the IUPAC InCHi of massbank record
	## Input : $record
	## Output : $inchi
	## Usage : my ( $inchi ) = getInchiFromString ( $record ) ;
	
=cut
## START of SUB
sub getInchiFromString {
    ## Retrieve Values
    my $self = shift ;
    my ( $record ) = @_;
    my ( $inchi ) = ( undef ) ;
    
    if ( defined $record ) {
    	my @tmp = split(/\n/, $record) ;
    	foreach my $field (@tmp) {
    		if ($field =~/CH\$IUPAC:\s+(.+)/) {
				$inchi = $1;
    		}
    	}
    	# for db field
    }
    else {
    	croak "Can't work with a undef / none existing massbank handler\n" ;
    }
    
    return ($inchi) ;
}
### END of SUB

=head2 METHOD getExactMzFromString

	## Description : get the exact mass of massbank record
	## Input : $record
	## Output : $exactMass
	## Usage : my ( $exactMass ) = getExactMzFromString ( $record ) ;
	
=cut
## START of SUB
sub getExactMzFromString {
    ## Retrieve Values
    my $self = shift ;
    my ( $record ) = @_;
    my ( $exactMass ) = ( undef ) ;
    
    if ( defined $record ) {
    	my @tmp = split(/\n/, $record) ;
    	foreach my $field (@tmp) {
    		if ($field =~/CH\$EXACT_MASS:\s+(.+)/) { 
				$exactMass = $1;
    		}
    	}
    	# for db field
    }
    else {
    	croak "Can't work with a undef / none existing massbank handler\n" ;
    }
    
    return ($exactMass) ;
}
### END of SUB


=head2 METHOD getPrecursorTypeFromString

	## Description : get the precursor type of massbank record
	## Input : $record
	## Output : $precursorType
	## Usage : my ( $precursorType ) = getPrecursorTypeFromString ( $record ) ;
	
=cut
## START of SUB
sub getPrecursorTypeFromString {
    ## Retrieve Values
    my $self = shift ;
    my ( $record ) = @_;
    my $id = undef ;
    my $precursorType = undef ;
    my $precursorType_first  = undef ;
    my $ionType_first  = undef ;
    my $precursorType_optionnal = undef ;
    
    if ( defined $record ) {
    	my @tmp = split(/\n/, $record) ;
    	foreach my $field (@tmp) {
    		if ($field =~/ACCESSION:\s+(.+)/) { 
				$id = $1;
    		}
    		if ($field =~/RECORD_TITLE:\s+(.+)/) { 
				my @title = split(/;/, $1) ;
				$precursorType_optionnal = $title[-1] ;
				$precursorType_optionnal =~ s/\s//g ;
    		}
    		if ($field =~/PRECURSOR_TYPE(.+)/) {
				$precursorType_first = $1;
				last;
    		}
    		if ($field =~/ION_TYPE(.+)/) {
				$ionType_first = $1;
				last;
    		}
    	}
    	# for db field
    }
    else {
    	croak "Can't work with a undef / none existing massbank handler\n" ;
    }
    
    ## manage undef precursor/ion type field 
#    print "ID:$id-//-$precursorType_first-//-$ionType_first-//-$precursorType_optionnal\n" ;
    if (defined $precursorType_first) {
    	$precursorType = $precursorType_first ;
    }
    elsif ( (!defined $precursorType_first) and (defined $ionType_first) ) {
    	$precursorType = $ionType_first ;
    }
    elsif ( (!defined $precursorType_first) and (!defined $ionType_first) and (defined $precursorType_optionnal) ) {
    	$precursorType = $precursorType_optionnal ;
    }
    else {
    	$precursorType = 'NA' ;
    }
    
    return ($precursorType) ;
}
### END of SUB

=head2 METHOD getMsTypeFromString

	## Description : get the MS type of massbank record
	## Input : $record
	## Output : $msType
	## Usage : my ( $msType ) = getMsTypeFromString ( $record ) ;
	
=cut
## START of SUB
sub getMsTypeFromString {
    ## Retrieve Values
    my $self = shift ;
    my ( $record ) = @_;
    my ( $msType ) = ( undef ) ;
    
    if ( defined $record ) {
    	my @tmp = split(/\n/, $record) ;
    	foreach my $field (@tmp) {
    		if ($field =~/AC\$MASS_SPECTROMETRY:\s+MS_TYPE\s+(.+)/) { 
				$msType = $1;
    		}
    	}
    	# for db field
    }
    else {
    	croak "Can't work with a undef / none existing massbank handler\n" ;
    }
    
    return ($msType) ;
}
### END of SUB

=head2 METHOD getChemNamesFromString

	## Description : get lits of names of a massbank record
	## Input : $record
	## Output : $names
	## Usage : my ( $names ) = getChemNamesFromString( $record ) ;
	
=cut
## START of SUB
sub getChemNamesFromString {
	## Retrieve Values
    my $self = shift ;
    my ( $record ) = @_ ;
    
    my @names = () ;
    if ( defined $record ) {
    	my @tmp = split(/\n/, $record) ;
    	foreach my $field (@tmp) {   		
    		if ($field =~/CH\$NAME: (.*)/) { 
    			push(@names, $1 ) ;  }
	    	else { next ; }
    	}
    }
    else {
    	croak "Can't work with a undef / none existing massbank record (string)\n" ;
    }
    return(\@names) ;
}
## END of SUB





=head2 METHOD getMassBankHandler

	## Description : get a massbank handler from a file
	## Input : $record
	## Output : $massbankHandler
	## Usage : my ( $massbankHandler ) = getMassBankHandler ( $record ) ;
	
=cut
## START of SUB
sub getMassBankHandler {
    ## Retrieve Values
    my $self = shift ;
    my ( $record ) = @_;
    my ( $massbankHandler ) = ( undef ) ;
    
    ## TODO...
    
    return ($massbankHandler) ;
}
### END of SUB

=head2 METHOD get_annotations_data

	## Description : permet de recuperer tous les champs d'un object massbank
	## Input : $ms_file
	## Output : $features
	## Usage : my ( $features ) = get_annotations_data( $ms_file ) ;
	
=cut
## START of SUB
sub get_annotations_data {
	## Retrieve Values
    my $self = shift ;
    my ( $ms_file ) = @_ ;
    
    my @features = () ;
    if ( ( defined $ms_file ) and ( -e $ms_file )) {
    	open(MS, "<$ms_file") or die "Cant' read the file $ms_file\n" ;
    	while ( my $field = <MS> ){
    		chomp $field ;    		
    		if ($field =~/PK\$ANNOTATION:(.*)/) { push( @features, $1) ;  }
	    	else { next ; }
    	}
    	close(MS) ;
    	# for db field
    }
    else {
    	croak "Can't work with a undef / none existing massbank file\n" ;
    }
    return(\@features) ;
}
## END of SUB

=head2 METHOD get_links_data

	## Description : permet de recuperer tous les champs d'un object massbank
	## Input : $ms_file
	## Output : $features
	## Usage : my ( $features ) = get_annotations_data( $ms_file ) ;
	
=cut
## START of SUB
sub get_links_data {
	## Retrieve Values
    my $self = shift ;
    my ( $ms_file ) = @_ ;
    
    my %features = () ;
    my $control = 0 ;
    
    my ( @CAS, @KEGG, @PUBCHEM ) = ((), (), ()) ;
    
    if ( ( defined $ms_file ) and ( -e $ms_file )) {
    	open(MS, "<$ms_file") or die "Cant' read the file $ms_file\n" ;
    	while ( my $field = <MS> ){
    		chomp $field ;    		
    		if ($field =~/CH\$LINK: CAS (.*)/) { push (@CAS, $1) ; $control++; }
    		elsif ($field =~/CH\$LINK: KEGG (.*)/) { push (@KEGG, $1) ; $control++; }
    		elsif ($field =~/CH\$LINK: PUBCHEM CID (.*)/) { push (@PUBCHEM, $1) ; $control++; }
    		## others !!?
    		
	    	else { next ; }
    	}
    	close(MS) ;
    	# for db field
    }
    else {
    	croak "Can't work with a undef / none existing massbank file\n" ;
    }
    
    $features{'CAS'} = \@CAS ;
    $features{'KEGG'} = \@KEGG ;
    $features{'PUBCHEM'} = \@PUBCHEM ;
    
    return(\%features) ;
}
## END of SUB

=head2 METHOD get_ms_record_links_data

	## Description : permet de recuperer tous les champs d'un object massbank
	## Input : $ms_file
	## Output : $features
	## Usage : my ( $features ) = get_ms_record_links_data( $ms_file ) ;
	
=cut
## START of SUB
sub get_ms_record_links_data {
	## Retrieve Values
    my $self = shift ;
    my ( $ms_file ) = @_ ;
    
    ## Internal reference for MASSBANK and RESPECT
    
    my @massbank_id = ( 'TUE', 'GLS', 'AU', 'MSJ', 'ML','FIO', 'UF', 'CO', 'UO', 'TT', 'OUF', 'MCH', 'NU', 'KNA', 'MT', 'CE', 'KO', 'KZ', 'JEL', 'JP', 'PR', 'BML', 'CA', 'TY', 'PB', 'FU', 'EA', 'UT', 'BSU', 'WA' ) ;
    my @respect_id = ( 'PS', 'PT', 'PM' ) ;
    
    my $dabase_used = undef ;
    my %db = ( 'accession' => undef, 'name' => undef ) ;
    my $control = 0 ;
    
    if ( $ms_file ) {
    	my $filename = basename("$ms_file",  ".txt");
    	
    	if ( $filename =~ /(\w+)$/ ) { # keep only record id (0001-PS0002 => PS0002 or BJ0045 => BJ0045) 
    		$db{'accession'} = $1 ;
    		$control++ ;
    		if ( ( defined $db{'accession'} ) and ( $db{'accession'} =~ /(\D+)(\d+)/) ) {
    			my ($key, $eval) = ($1, 0) ;
    			foreach (@respect_id) { if ($_ eq $key) { $db{'name'} = 'RESPECT' ; $eval = 1 ; last ; } }
    			foreach (@massbank_id) { if ($_ eq $key) { $db{'name'} = 'MASSBANK' ; $eval = 1 ; last ; } }
    			if ( $eval == 0 ){ 	carp "The following key ($key) for $db{'accession'} has an unknown reference (not a Massbank or ReSpect source)\n" ; }
    		}
    	}
    }
    if ($control == 0) { %db = () ;  }
	return(\%db) ;
}
## END of SUB


1 ;


__END__

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

 perldoc parser::chem::massbank.pm

=head1 Exports

=over 4

=item :ALL is ...

=back

=head1 AUTHOR

Franck Giacomoni E<lt>franck.giacomoni@clermont.inra.frE<gt>

=head1 LICENSE

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=head1 VERSION

version 1 : 25 / 06 / 2013

version 2 : ??

=cut