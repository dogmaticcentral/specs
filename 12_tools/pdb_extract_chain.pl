#! /usr/bin/perl -w
use IO::Handle;         #autoflush
# FH -> autoflush(1);

defined ( $ARGV[0]  ) ||
    die "Usage: pdb_chain_rename.pl   <pdb_file>   [<chain_name>].\n";

$pdbfile = $ARGV[0];
if ( defined $ARGV[1] ) {
    $query_chain_name =$ARGV[1] ;
} else {
    $query_chain_name ="" ;
}

open ( IF, "<$pdbfile") ||
    die "Cno $pdbfile: $!.\n";


    
while ( <IF> ) {

    last if ( /^ENDMDL/);
    next if ( ! /^ATOM/  );
    $line = $_;
    $chain_name = substr ( $line,  21, 1) ;
    if ( ! $query_chain_name ||   $chain_name eq " " ||
	$chain_name eq $query_chain_name ) {
	print $line;
    }
}

close IF;
