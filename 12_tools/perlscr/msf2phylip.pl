#! /usr/bin/perl -w
use IO::Handle;         #autoflush
# FH -> autoflush(1);
$length = 0;
$no_names = 0;
$first_name = "";
$first = 1;
$all = "";
$i =0;
%used = ();
@names = ();
while ( <> ) {
    last if ( /\/\// );
    if ( /MSF:/ ) {
	@aux = split;
	$almt_length = $aux[1];
    } elsif ( /Name:/ ) {
	@aux = split;
	push @names, $aux[1];
   }
}

#main problem with phylip: it want the names no longer than 10 characters
foreach $name ( @names ) {
    $l = length($name);
    if ( $l == 10 ) {
	$short_name{$name} = $name;
    } elsif ( $l > 10 ) {
	$short_name{$name} = substr $name, 0, 10;
    } else {
	$short_name{$name} = $name;
	for ($i=$l; $i<10; $i++) {
	    $short_name{$name} .= " ";
	}
    }
    $ctr = 2;
    while ( defined $used{$short_name{$name}} && $ctr < 10)  {
	substr ($short_name{$name}, $l-1, 1) = $ctr;
	$ctr ++;
    }
    if ( $ctr == 10 ) {
	die "In msf2phylip.pl: error resolving name clashes ($name).\n";
    }
}
	
while ( <> ) {
    if ( /\S/ ) {
	chomp;
	$line = $_;
	$line =~ s/\./\-/g;
	@aux = split ' ', $line;
	$name = $aux[0];
	$name =~ s/\./\-/g;
	
	if ( ! defined $found{$name} ) {
	    $found{$name } = 1;
	    $all .= $short_name{$name}." ";
	    $no_names ++;
	    if ( $first ) { 
		$first = 0; 
		$first_name = $name;
	    } 
	} elsif ( $name =~ $first_name ) {
	    $all .= "\n";
	}
	foreach $piece ( @aux[1..$#aux] ) {
	    $all .=  "$piece ";
	    $length += length $piece;
	}
	$all .= "\n";
    
    }
}

print "     $no_names     $almt_length     \n";
print $all;
