#!/usr/bin/perl -w
use strict;

my @filecontents;
my $line;
my @atom;
my $atom;
my @het;
my $het;
my $dist;
my $name;
my $distance;
print "Enter the location of a .pdb file:  ";
my $pdbFile = <STDIN>; # grab pdb file from user
$pdbFile =~ s/\s+$//; # trim right side of name, basically chomp, but I think it has been more consistent
my $filename = '/Users/Hubble_Space_Telescope/Desktop/somePymolStuff/testFolder/testing.pml'; # file to write to
my %close;


#
# Example PDB file format for atom type (PERL numbering start w/ 0)
#
#          1         2         3         4         5         6         7         8
#012345678901234567890123456789012345678901234567890123456789012345678901234567890
#
#ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92      A1   N
#ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85      A1   C
#HETATM15786  U   IUM A1102      18.227  57.101 113.666  0.30 28.22           U
#
#Record Format
#
#COLUMNS      DATA TYPE        FIELD      DEFINITION
#------------------------------------------------------
# 1 -  6      Record name      "ATOM    "
# 7 - 11      Integer          serial     Atom serial number.
#13 - 16      Atom             name       Atom name.
#17           Character        altLoc     Alternate location indicator.
#18 - 20      Residue name     resName    Residue name.
#22           Character        chainID    Chain identifier.
#23 - 26      Integer          resSeq     Residue sequence number.
#27           AChar            iCode      Code for insertion of residues.
#31 - 38      Real(8.3)        x          Orthogonal coordinates for X in Angstroms
#39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in Angstroms
#47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in Angstroms
#55 - 60      Real(6.2)        occupancy  Occupancy.
#61 - 66      Real(6.2)        tempFactor Temperature factor.
#77 - 78      LString(2)       element    Element symbol, right-justified.
#79 - 80      LString(2)       charge     Charge on the atom.

open (INPUT, "< $pdbFile") || die "Could not open file\n";
chomp(@filecontents = <INPUT>);
#my $moleculeName = substr $pdbFile, -8, 4; #method 1, grabs 8 char from end, until 4 from end (.pdb)
my $moleculeName;
if ($pdbFile =~ /\/([\d\w]{4})\./){ #method 2, grabs 4 char between '/' and '.', aka the /(2hck).pdb bit
	$moleculeName = $1;
} #I think the other way works fine if you don't chomp filename, but I'm not sure if there's a reason not to do this way
printf "$moleculeName\n"; #just prints the molecule name in command line; not essential

foreach $line (@filecontents) {
        if (($line =~ /^ATOM/) && ($line =~ /\bCA/)) {
       # if ($line =~ /^ATOM/)  {
                push(@atom, $line);
        } elsif (($line =~ /^HETATM/) && !($line =~ /HOH/)) { #excluding water seems to have improved results, less clutter
                push(@het, $line);
        }
}

open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
print $fh "fetch $moleculeName, async=0\n hi ev\n sh cartoon\n bg white \n"; #begin pymol script
#A few lines of the beginning, hides everything but cartoon; white bg is nice for screenshots
my $count = 0; #used for labeling selections
foreach $atom (@atom) {
        my $x = substr($atom, 31, 8);
        my $y = substr($atom, 39, 8);
        my $z = substr($atom, 47, 8);
        my $resname = substr($atom, 17, 3);
        my $chainid = substr($atom, 21, 1);
        my $resnum  = substr($atom, 22, 4);

        foreach $het (@het) {
                my $a = substr($het, 31, 8);
                my $b = substr($het, 39, 8);
                my $c = substr($het, 47, 8);
        			my $hetresname = substr($het, 17, 3);
       			my $hetchainid = substr($het, 21, 1);
     			   	my $hetresnum  = substr($het, 22, 4);
        	
                my $dist = sqrt( ($x - $a)**2 + ($y - $b)**2 + ($z - $c)**2 );
                if ($dist <= 5) {
        				  $close{"$hetresname:$hetchainid:$resname:$chainid:$resnum"}  = $dist;
        	   	
        	   	      	  $resnum =~ s/^\s+//;
        	            print $fh "sele sele_$count, $chainid/$resname`$resnum/\n";
        	         	   print $fh "sele sele_$count, sele_$count expand 5\n\n"; # 5 was $dist
        	           	   # select the predicted active site (and the ligand for now, recolored later)
        	         	   $count++;
        	   
              	}
                			
        }
}
my $i = 0;
print $fh "set cartoon_side_chain_helper, on\n";
print $fh "sele binding, ("; 
while ($i <  $count) {
			print $fh "sele_$i, ";
			$i++;
}
print $fh ")\n";

$i = 0;
print $fh "delete ("; 
while ($i <  $count) {
			print $fh "sele_$i, ";
			$i++;
}
print $fh ")\n";
print $fh "sele nonbind, !binding and $moleculeName\n";
print $fh "color gray40, nonbind\n"; # recolor things not in the selection
print $fh "desele\n";
print $fh "sele ligand, binding and hetatm\n"; #selects the hetatms from the selection colored in foreach loop
print $fh "sele ligand, ligand expand 2\n";
print $fh "sh sticks, binding\n";
print $fh "color red, ligand\n"; #red seems to be a good color to make it stick out well
print $fh "desele\n";
print $fh "dist distance, ligand, binding, mode=2\n";
print $fh "select closeh2o, resn HOH w. 5 of ligand\n";
print $fh "sh spheres, closeh2o\n";
print $fh "set sphere_scale, 0.25\n";
print $fh "color aqua, closeh2o\n";
print $fh "color orange, distance\n";
print $fh "desele\n";
print $fh "orient binding\n";
print $fh "zoom binding, 15 \n"; #trying to come up with a good starting viewpoint; this seems okay. maybe room to improve
system('open /Users/Hubble_Space_Telescope/Desktop/somePymolStuff/testFolder/testing.pml');	 #runs pymol script.
close $fh;

	
while (($name, $distance) = each %close) {
printf "%s %5.3f\n",$name, $distance;

					}

					
					
					

