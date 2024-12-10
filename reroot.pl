#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib $RealBin;
use MOTree;

my $usage = <<END
Usage: reroot -midpoint < tree > outtree
       reroot split1file split2file < tree > outtree
           where each split file is a list of identifiers, one per line, and
           it tries to reroot so that the two sets of identifiers are on
           different sides.
END
  ;

{
  my ($debug, $midpoint);
  die $usage
    unless GetOptions('debug' => \$debug,
                      'midpoint' => \$midpoint);
  my %split1 = ();
  my %split2 = ();
  if (defined $midpoint) {
    die $usage unless @ARGV == 0;
  } else {
    die $usage unless @ARGV == 2;
    my ($file1,$file2) = @ARGV;
    open(IN,"<",$file1) || die "Cannot read $file1\n";
    while(<IN>) {
      chomp;
      $split1{$_} = 1;
    }
    close(IN) || die "Error reading $file1\n";
    open(IN,"<",$file2) || die "Cannot read $file2\n";
    while(<IN>) {
      chomp;
      $split2{$_} = 1;
    }
    close(IN) || die "Error reading $file2\n";
  }
  my $nTrees = 0;
  while(my $tree = MOTree::new(fh => \*STDIN)) {
    $nTrees++;
    if ($midpoint) {
      $tree->rerootMidpoint();
    } else {
      my %leaves = map {$tree->id($_) => $_} $tree->get_leaf_nodes();
      my @split1 = ();
      my @split2 = ();
      foreach my $name (keys %split1) {
        die "name '$name' is in both splits!" if exists $split2{$name};
        push @split1, $leaves{$name} if exists $leaves{$name};
      }
      foreach my $name (keys %split2) {
        push @split2, $leaves{$name} if exists $leaves{$name};
      }
      my $lca1 = $tree->lca(@split1);
      my $lca2 = $tree->lca(@split2);
      die "Invalid split -- both have the same least common ancestor"
        if ($lca1 == $lca2);
      my $rootat = $lca1;
      $rootat = $lca2 if $lca1 == $tree->get_root_node;
      print STDERR "rootat $rootat\n" if $debug;
      $tree->reroot($rootat);
    }
    print $tree->toNewick(),"\n";
  }
  die "No tree in STDIN" unless $nTrees > 0;
}
