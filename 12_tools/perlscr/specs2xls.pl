#! /usr/bin/perl -w
 
# to install perl Spreadsheet module
# perl -MCPAN -e "install 'Spreadsheet::WriteExcel'"
# doc: 
# http://search.cpan.org/dist/Spreadsheet-WriteExcel/lib/Spreadsheet/WriteExcel.pm

use Spreadsheet::WriteExcel;                            

(@ARGV >= 2) ||
    die "Usage:  cbcvg.pl   <specs score file>  <output name>  \n"; 


($ranks_file, $outname) = @ARGV;

@method_names = ('rvet', 'entr', 'ivet', 'majf', 'valdar', 'rvn', 'pheno', 'carb', 
		 'entr_s', 'rvs', 'majf_s', 'ivet_s', 'phen_h', 'phen_n');
$server_version = 1;

my $workbook = Spreadsheet::WriteExcel->new("$outname.xls");
$worksheet   = $workbook->add_worksheet();                  
$worksheet->keep_leading_zeros();

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
$xl_color[0] = $workbook->set_custom_color(8, $red*255, $green*255, $blue*255 );	

$bin_size = $C1/$N;
for ( $ctr=1; $ctr <= int ($COLOR_RANGE/$N); $ctr++ ) {

    $ratio =  ( int ( 100*($bin_size- $ctr+1)/$bin_size) ) /100;
    $red   = $ratio;
    $green = $blue = 0;
		 
    $color[$ctr] = "[$red, $green, $blue]"; 
    $color_name[$ctr] = "c$ctr";

    $xl_color[$ctr] = $workbook->set_custom_color($ctr+8, $red*255, $green*255, $blue*255 );	
    #printf "  %d   %d    %d , $ctr, $xl_color[$ctr] \n",$red*255, $green*255, $blue*255;
}

for ( $ctr= int ($COLOR_RANGE/$N)+1 ; $ctr <= $COLOR_RANGE; $ctr++ ) {

    $ratio =  ( $ctr -  $COLOR_RANGE/$N)/ ($COLOR_RANGE*($N-1)/$N);
    $red = $ratio;
    $green = $blue = $red;

    $xl_color[$ctr] = $workbook->set_custom_color($ctr+8, $red*255, $green*255, $blue*255 );	
    #printf "  %d   %d    %d , $ctr, $xl_color[$ctr] \n",$red*255, $green*255, $blue*255;

  
    $color[$ctr] = "[$red, $green, $blue]"; 
    $color_name[$ctr] = "c$ctr";
}



##################################################
# input/ouput

open (RANKS_FILE, "<$ranks_file") || 
    die "cno $ranks_file\n";
    

$format_centered = $workbook->add_format();
$format_centered->set_align('center');
$format_right = $workbook->add_format();
$format_right->set_align('right');
$format_float = $workbook->add_format();
$format_float->set_num_format('0.00');
#
$format_hdr = $workbook->add_format(); # Add a format
$format_hdr->set_bold();
$format_hdr->set_border();
$format_hdr->set_align('center');

$line = 1;
$last_column = 0;
@header = ();

while ( <RANKS_FILE> ) {
    next if ( !/\S/ );
    if ( /\%/ ){ # find the column that the specified method is in
	@aux = split;
	shift @aux;
	@header = @aux;

	if ( $server_version) {
	    @aux2 = ();
	    for $hdr_field (@aux) {
		if ($hdr_field eq 'surf') {
		    push @aux2, "surface";
		} elsif ( grep ($_ =~ $hdr_field, @method_names ) ) {
		    push @aux2, " conservation score ($hdr_field) ";
		} else {
		    push @aux2, $hdr_field;
		}
	    }
	    @header = @aux2;
	}

	# header
 	$last_column = @header;
	$column =  chr (ord('A') + $last_column-1);
        $worksheet->set_column("$column:$column", 20);
	$worksheet->write_row(0, 0, \@header, $format_hdr);

   } else {
        # print;
	chomp;
	@aux = split;
	$worksheet->write("A$line", '',  $format);

	foreach $i (0 .. $#aux) {
	    $column =  chr (ord('A') + $i);
	    if ( grep ($_ =~$header[$i], @method_names ) ) {

		$cvg =  $aux[$i];
		$color_index = int ($cvg*$COLOR_RANGE );
		# the first column is the color strip
		$format = $workbook->add_format(
		    bg_color => $xl_color[$color_index],
		    pattern  => 1,
		    border   => 1
		    );
		$format ->set_num_format('0.00');
		$worksheet->write("$column$line", $aux[$i], $format);

	    } elsif ($header[$i] eq 'gaps') {
		$worksheet->write("$column$line", $aux[$i], $format_float);   
 
	    } elsif ($header[$i] eq 'surface') {
		if ( $aux[$i] == 1) {
		    $worksheet->write("$column$line", "surface", $format_centered);  
		} 
 
	    } else {

		# generic behavior
		if ( $aux[$i] =~ /\d/ ) {
		    $worksheet->write("$column$line", $aux[$i], $format_right);  
		} else {  
		    $worksheet->write("$column$line", $aux[$i], $format_centered);  
		}  

	    }
	}

    }
    $line++;
}


close RANKS_FILE;



##################################################
##################################################
##################################################
##################################################
# colorbar

$column =  $last_column + 2;
$row = 2;
$worksheet->write($row, $column, "conserved", $format_right);   
$row = 2+$COLOR_RANGE;
$worksheet->write($row, $column, "variable", $format_right);   



$column =  $last_column + 3;
$row = 0;

$worksheet->write($row, $column, "conservation colorbar", $format_hdr);   

for ($ctr = 0; $ctr <= $COLOR_RANGE; $ctr++) {
    $format = $workbook->add_format(
	bg_color => $xl_color[$ctr],
	pattern  => 1,
	border   => 1
	);
    $row = $ctr+2;
    $worksheet->write($row, $column, '',  $format);
}
