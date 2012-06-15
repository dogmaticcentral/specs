#!/usr/bin/perl -w 
# shorten the names staring with gi (like in entres\z return)


while (<STDIN>) {
    if (0) {
    } elsif ( />\s*gi\|(\d+)\|/ ) { # for gi names
	print ">$1\n";
    } elsif ( />.*\|(\w+?)\s/ ) { # for uniprot names
	print "> $1\n";
    } else {
	print;
    }
 }


