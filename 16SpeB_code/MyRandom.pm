#!/usr/bin/perl -w
use strict;

package MyRandom;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(rand srand);

#use Math::Random::MT::Perl;

#EXAMPLE:
#
#use MyRandom;
#
#my $count = 100;
#my $srand = 10;
#
#my $rng = new MyRandom;
#$rng->srand($srand);
#
#for(my $i = 0; $i < $count; $i++) {
#    print $rng->rand() . "\n";
#}

sub new 
{
	my ($class, $rndseed, $generator) ;
	if(scalar(@_) == 1) {
		($class) = @_;
		$rndseed = undef;
		$generator = undef;
	}
	elsif(scalar(@_) == 2) {
		($class, $rndseed) = @_;
		#$generator = Math::Random::MT::Perl->new($rndseed);
	}
	else {
		die "Invalid number of parameters";
	}

	my $self = {
		'RAND_SEED' => $rndseed,
		'GENERATOR' => $generator,
	};
	bless($self, $class);
	return $self;
}

sub rand
{
	my $self = shift;
	#return $self->{GENERATOR}->rand();
	return rand();
}

sub srand
{
	my ($self, $seed) = @_;

	$self->{RAND_SEED} = $seed;
	#$self->{GENERATOR} = Math::Random::MT::Perl->new($seed);
	srand($seed);
}


1;

__END__
