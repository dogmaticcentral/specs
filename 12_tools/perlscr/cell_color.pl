#!/usr/bin/perl -w

 

(@ARGV ==3 ) ||
    die "Usage:  cbcvg.pl  <method [rvet|majf|entr]>  <specs score file>   <output name>  \n"; 
($method, $ranks_file,  $output_file) = @ARGV;

##################################################
#set the pallette:
$COLOR_RANGE = 20;
$green = $blue = $red = 0;


$N = 5;
$C1 = $COLOR_RANGE-1;

$red = 1.00;
$green =  0.83;
$blue =    0.17;
	
$red = sprintf "%3d", $red*255;
$blue = sprintf "%3d", $blue*255;
$green = sprintf "%3d", $green*255;

$color[0] = "$red, $green, $blue"; 

$bin_size = $C1/$N;
for ( $ctr=1; $ctr <= int ($COLOR_RANGE/$N); $ctr++ ) {

    $ratio =  ( int ( 100*($bin_size- $ctr+1)/$bin_size) ) /100;
    $red   = $ratio;
    $green = $blue = 0;
	
    $red = sprintf "%3d", $red*255;
    $blue = sprintf "%3d", $blue*255;
    $green = sprintf "%3d", $green*255;

    $color[$ctr] = "$red, $green, $blue"; 
 
}

for ( $ctr= int ($COLOR_RANGE/$N)+1 ; $ctr <= $COLOR_RANGE; $ctr++ ) {

    $ratio =  ( $ctr -  $COLOR_RANGE/$N)/ ($COLOR_RANGE*($N-1)/$N);
    $red = $ratio;
    $green = $blue = $red;

    $red = sprintf "%3d", $red*255;
    $blue = sprintf "%3d", $blue*255;
    $green = sprintf "%3d", $green*255;

		 
    $color[$ctr] = "$red, $green, $blue"; 
}


##################################################
# input


open (RANKS_FILE, "<$ranks_file") || 
    die "cno $ranks_file\n";
    

$method_column = -1;
@cvg = ();

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
	push @cvg, $aux[$method_column];
	$color_index =int ($aux[$method_column]*$COLOR_RANGE );
    }
}

close RANKS_FILE;

##################################################
# output

# open the output file


$ctr = 0;
foreach  ( @cvg) {

    if ( ! ( $ctr % 200 )  ) {
	$ctr && close FPTR;
	$filename = $output_file.".$ctr";
	open (FPTR, ">$filename") || die "cno $filename\n";
	
    }

    $color_index = int ($_*$COLOR_RANGE );
    $ctr++;

    print FPTR "Cells ($ctr, 2).Select\n";
    print FPTR "With Selection.Interior\n";
    print FPTR "    .Color = RGB( $color[$color_index] )\n";
    print FPTR "    .Pattern = xlSolid\n";
    print FPTR "End With\n";
    
     
}



close FPTR; 





