#!/usr/bin/perl -w

 

(defined $ARGV[3]) ||
    die "Usage:  cbcvg.pl  <method [rvet|majf|entr]>  <specs score file>  <pdb_file_full_path>  <output name> [<chain> and/or -r and/or -b] \n"; 
($method, $ranks_file, $pdb_file, $output_file) = @ARGV;

##################################################
#set the pallette:
$COLOR_RANGE = 20;
$green = $blue = $red = 0;


$N = 5;
$C1 = $COLOR_RANGE-1;

$red = 1.00;
$green =  0.83;
$blue =    0.17;
$color[0] = "[$red, $green, $blue]"; 
$color_name[0] = "c0";

$bin_size = $C1/$N;
for ( $ctr=1; $ctr <= int ($COLOR_RANGE/$N); $ctr++ ) {

    $ratio =  ( int ( 100*($bin_size- $ctr+1)/$bin_size) ) /100;
    $red   = $ratio;
    $green = $blue = 0;
		 
    $color[$ctr] = "[$red, $green, $blue]"; 
    $color_name[$ctr] = "c$ctr";

}

for ( $ctr= int ($COLOR_RANGE/$N)+1 ; $ctr <= $COLOR_RANGE; $ctr++ ) {

    $ratio =  ( $ctr -  $COLOR_RANGE/$N)/ ($COLOR_RANGE*($N-1)/$N);
    $red = $ratio;
    $green = $blue = $red;
		 
    $color[$ctr] = "[$red, $green, $blue]"; 
    $color_name[$ctr] = "c$ctr";
}


##################################################
# input
$chain = "";
$reverse = 0;
$backbone = 0;

for  $argctr ( 4 .. 5 ) {
    if ( defined  $ARGV[ $argctr ] ){
	if ( $ARGV[$argctr ] eq "-r" ) {
	    $reverse  = 1;
	} elsif ( $ARGV[ $argctr ] eq "-b" ) {
	    $backbone  = 1;
	} else {
	    $chain =  $ARGV[4];
	}
    }
}



open (RANKS_FILE, "<$ranks_file") || 
    die "cno $ranks_file\n";
    

$method_column = -1;

while ( <RANKS_FILE> ) {
    next if ( !/\S/ );
    if ( /\%/ ){
	@aux = split;
	shift @aux;
	for ($ctr=0; $ctr< $#aux; $ctr++) {
	    if ($aux[$ctr] eq $method ) {
		$method_column = $ctr;
		last;
	    }
	}
    } elsif ($method_column > -1) {
	chomp;
	@aux = split;
	$pdb_id = $aux[1];
	next if ($pdb_id =~ '-' );
	$cvg{$pdb_id} = $aux[$method_column];
	if ( $reverse ) {
	    $cvg{$pdb_id} = 1 - $cvg{$pdb_id};
	}
    }
}


close RANKS_FILE;

##################################################
# output
format FPTR = 
load @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< , struct_name 
     $pdb_file
zoom complete=1
bg_color white
color white, struct_name
show spheres, struct_name
hide lines, struct_name
#spacefill
.

# open the output file
if  ($reverse ) {
    $filename = $output_file.".rev";
} else {
    $filename = $output_file;
}

open (FPTR, ">$filename") || die "cno $filename\n";

write FPTR;

for $ctr ( 0 .. $#color ) {
    print  FPTR "set_color $color_name[$ctr] = $color[$ctr]\n";
}

foreach $pos ( keys %cvg ) {

    $color_index = int ($cvg{$pos}*$COLOR_RANGE );
    
    print FPTR "color $color_name[$color_index], resid $pos ";
    if ( $chain ){
	print  FPTR  "and chain $chain";
    }
    print  FPTR "\n";
     
}

=pod
# load in the  structure as a "movie"
# color each "state" according to the changing level of conservation

fetch 1nmr
set all_states
 
# the object has 20 states, so we can set separate line colors
# for each state as follows:
for a in range(1,21): cmd.set("line_color","auto","1nmr",a)
# in my case it will probably be something like

for a in range(1,21): cmd.set("color",  $color_name[$color_index], "resid $pos",a)


# creating text on the screen
pseudoatom foo
hide everything, foo
set label_font_id ,12
set label_size, -4  # angstroms
label foo, "Once upon a time..."

=cut



close FPTR; 




