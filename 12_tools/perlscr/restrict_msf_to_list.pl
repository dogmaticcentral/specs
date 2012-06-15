#! /usr/bin/perl -w
use IO::Handle;         #autoflush
# FH -> autoflush(1);

 
defined $ARGV[1]  ||
    die "Usage: $0  <msf_file_name>  <list of sequences (file name)>.\n"; 


($msf, $file_name) =  @ARGV;


@protected_names = ();
open (IF, "<$file_name") || die "Cno $file_name: $!\n";

while ( <IF>) {
    chomp;
    $name =  $_;
    $name =~ s/\s//g;
    push @protected_names, $name;
}
close IF;


open ( MSF, "<$msf") ||
    die "Cno $msf: $!\n";

# read in the msf file:
while ( <MSF> ) {
    last if ( /\/\//);
}

%sequence = ();
@names = ();
do {
    if ( /\w/ ) {
	@aux = split;
	$name = $aux[0];
	$aux_str = join ('', @aux[1 .. $#aux] );
	if ( defined $sequence{$name} ) {
	    $sequence{$name} .= $aux_str;
	} else {
	    $sequence{$name}  = $aux_str;
	    push @names, $name;
	}
		
    } 
} while ( <MSF>);


# turn the msf into a table (first index= sequence, 2nd index= position
$seq = 0;
@protected_seqs = ();
foreach $name ( @names ) {
    @aux = split '', $sequence{$name};
    foreach $pos ( 0 .. $#aux ) {
	$array[$seq][$pos] = $aux[$pos];
    }

    foreach $query_name (@protected_names) {
	if ( $name =~ $query_name && $query_name =~ $name) {
	    push @protected_seqs, $seq;
	    last;
	}
    }
    $seq++;
    
}

$no_seqs = $seq;   # number of seqs
$max_seq = $seq-1; # max index a seq can have
$max_pos = $#aux;  # max index a position can have

# sanity check:
$no_seqs || die "Error in sift.pl: no seqs found.\n"; 


foreach $pos ( 0 .. $max_pos ) {

    $delete[$pos] = 1;
    foreach $query_seq (@protected_seqs)  {
	if ( $array[$query_seq][$pos] !~ '[\.\-]' ) {
	    $delete[$pos] = 0;
	    last;
	}
    }
}

foreach $seq ( 0 .. $max_seq ) {
    $seq_new[$seq] = "";
}

foreach $seq ( 0 .. $max_seq ) {
    foreach $pos ( 0 .. $max_pos ) {
	if ( ! $delete[$pos] ) {
	    $seq_new[$seq] .= $array[$seq][$pos];
	}
    }
}

$deleted = 0;
foreach $pos ( 0 .. $max_pos ) {
    $deleted += $delete[$pos];
}

$new_length = length $seq_new[0];

@aux = split '\.', $msf;


print  "PileUp\n\n";
print  "            GapWeight: 30\n";
print  "            GapLengthWeight: 1\n\n\n";
printf   ("  MSF: %d  Type: P    Check:  9554   .. \n\n",$new_length) ;
foreach $name ( @names ) {
    printf  (" Name: %-40s   Len: %5d   Check: 9554   Weight: 1.00\n", $name, $new_length);
}
print  "\n//\n\n\n\n";

for ($j=0; $j  < $new_length; $j += 50) {
    $seq = 0;
    foreach $name ( @names ) {
	printf  "%-40s", $name;
	for ( $k = 0; $k < 5; $k++ ) {
	    if ( $j+$k*10+10 >= $new_length ) {
		printf  "%-10s ",   substr ($seq_new[$seq], $j+$k*10 );
		last;
	    } else {
		printf  "%-10s ",   substr ($seq_new[$seq], $j+$k*10, 10);
	    }
	}
	print  "\n";
	$seq++;
    } 
    print  "\n";
}

$max_pos ++;
#printf "removed $deleted columns   (out of $max_pos)\n";


