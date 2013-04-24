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
$bad_cvg       = 0.5;
$min_cvg       = 1.0;

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
	($cvg{$pdb_id}< $min_cvg) && ($min_cvg=$cvg{$pdb_id});
	if ( $reverse ) {
	    $cvg{$pdb_id} = 1 - $cvg{$pdb_id};
	}
    }
}
close RANKS_FILE;


($bad_cvg < $min_cvg) &&  ($bad_cvg = $min_cvg);

##################################################
# output


# open the output file
if  ($reverse ) {
    $filename = $output_file.".rev";
} else {
    $filename = $output_file;
}

open (FPTR, ">$filename") || die "cno $filename\n";
print  FPTR "load $pdb_file, the_whole_thing\n";
print  FPTR "zoom complete=1\n";
print  FPTR "bg_color white\n";
print  FPTR "hide everything\n";

if (!$chain) {
    print  FPTR "color white,  the_whole_thing\n";
    print  FPTR "show cartoon, the_whole_thing\n";
    print  FPTR "show spheres, the_whole_thing\n";
} else {
    print  FPTR "select chain$chain, the_whole_thing and chain $chain and polymer \n";
    print  FPTR "color white,  chain$chain \n";
    print  FPTR "show cartoon, chain$chain \n";
    print  FPTR "show spheres, chain$chain \n";
}



for $ctr ( 0 .. $#color ) {
    print  FPTR "set_color $color_name[$ctr] = $color[$ctr]\n";
}

@poorly_scoring = ();

foreach $pos ( keys %cvg ) {

    $color_index = int ($cvg{$pos}*$COLOR_RANGE );
    
    print FPTR "color $color_name[$color_index], resid $pos ";
    if ( $chain ){
	print  FPTR  "and chain $chain";
    }
    print  FPTR "\n";
    if ( $cvg{$pos} > $bad_cvg ) {
         push @poorly_scoring,$pos; 
    }      
}

$from = 0;
$slc   = "";
while ($from < int (@poorly_scoring) ) {
    $to = $from + 10;
    ($to > $#poorly_scoring)  &&  ($to=$#poorly_scoring);
    $reslist = join "+", @poorly_scoring[$from .. $to];
    if ($slc) {
       $slc .= "select poorly_scoring, poorly_scoring or resid $reslist\n"
    } else {
       $slc  = "select poorly_scoring, resid $reslist\n"
    }
    $from = $to+1;

} 


if ( $chain ) {

    if ($slc) {
      $slc  .= "select poorly_scoring, poorly_scoring and chain$chain\n";
        
    }
 
    print  FPTR "select heteroatoms, hetatm and not solvent \n";
    print  FPTR "select other_chains, not chain $chain \n";
    print  FPTR "select struct_water, solvent and chain $chain \n";
    print  FPTR "cartoon putty \n";
    print  FPTR "show  cartoon,  other_chains \n";
    print  FPTR "show  sticks,   heteroatoms \n";
    print  FPTR "show  spheres,  struct_water \n";
    print  FPTR "color palecyan, struct_water \n";
    print  FPTR "color lightteal, other_chains or heteroatoms \n";
    print  FPTR "zoom  chain$chain\n";
}

#$slc .= "hide spheres, poorly_scoring\n";

print  FPTR $slc;


print  FPTR "deselect \n";

close FPTR; 




