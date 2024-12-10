#!/usr/bin/perl -w
use strict;

die "Usage: GammaLogToPaup.txt FastTreeLogFiles > siteloglks.txt\n"
    . "   Creates a file in PAUP-like format for CONSEL's makermt --paup\n"
if @ARGV==0;

print join("\t","Tree","-lnL","Site","-lnL")."\n";

my $ntree = 1;
my $npos = undef;
foreach my $file (@ARGV) {
    open(FILE, "<", $file) || die "Cannot read $file";
    while(<FILE>) {
	chomp;
	if (m/^Gamma20LogLk/) {
	    my @F = split /\t/, $_;
	    print join("\t", $ntree, -$F[1])."\n";
	    $ntree++;
	} elsif (m/^Gamma20/) {
	    my @F = split /\t/, $_;
	    print join("\t", "", "", $F[1]+1, -$F[2])."\n" unless $F[1] eq "Site";
	}
    }
    close(FILE) || die "Error reading $file";
}
