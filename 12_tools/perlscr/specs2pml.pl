#!/usr/bin/perl -w

 

(defined $ARGV[3]) ||
    die "Usage:  $0  <method [rvet|majf|entr]>  <specs score file>  <pdb_file_full_path>  <output name> [<chain> and/or -r and/or -b] \n"; 
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
$pdb_id_column = -1;

while ( <RANKS_FILE> ) {
    next if ( !/\S/ );
    if ( /^\%/ ){
	@aux = split;
	shift @aux;
	for ($ctr=0; $ctr< $#aux; $ctr++) {
	    if ($aux[$ctr] eq $method ) {
		$method_column = $ctr;
	    } elsif ($aux[$ctr] eq "pdb_id") { 
		$pdb_id_column = $ctr;
	    }
	}
    } elsif ($method_column > -1 && $pdb_id_column > -1) {
	chomp;
	@aux = split;
	$pdb_id = $aux[$pdb_id_column];
	next if ($pdb_id =~ '-' || $pdb_id =~ '\.'  );
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
hide everything
.


# open the output file
if  ($reverse ) {
    $filename = $output_file.".rev";
} else {
    $filename = $output_file;
}

open (FPTR, ">$filename") || die "cno $filename\n";

write FPTR;

if (!$chain) {
    print  FPTR "color white, struct_name\n";
    print  FPTR "show spheres, struct_name\n";
} else {
    print  FPTR "select chain$chain, struct_name and chain $chain and polymer \n";
    print  FPTR "color white, chain$chain \n";
    print  FPTR "show spheres, chain$chain \n";
}



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


if ( $chain ) {

    print  FPTR "select heteroatom, hetatm and not solvent \n";
    print  FPTR "select other_chains, not chain $chain \n";
    print  FPTR "select struct_water, solvent and chain $chain \n";
    print  FPTR "cartoon putty \n";
    print  FPTR "show  cartoon, other_chains \n";
    print  FPTR "show  sticks, heteroatom \n";
    print  FPTR "show  spheres,  struct_water \n";
    print  FPTR "color palecyan, struct_water \n";
    print  FPTR "color lightteal, other_chains or heteroatom \n";
    print  FPTR "zoom chain$chain\n";
}


print  FPTR "deselect \n";

close FPTR; 




