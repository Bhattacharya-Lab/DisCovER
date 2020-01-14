#!/usr/bin/perl
# Credit: 
# 	1. Zheng,W., Wuyun,Q., et al. (2019) Detecting distant-homology protein structures by aligning deep neural-network based contact maps. PLOS Computational Biology 
#	2. Zhang,C. et al. DeepMSA: constructing deep multiple sequence alignment to improve contact prediction and fold-recognition for distant-homology pro- teins. Bioinformatics.
# Modified by: Sutanu

use Math::Trig;
use File::Basename;
use Cwd 'abs_path';


######## set your variables here
$DisCovER_path="/home/project/conThreader/DisCovER/DisCovER-master/";
$a3m2psiblast=$DisCovER_path."preprocessing/Target/a3m2psiblast_v2.pl";
#printf "$a3m2psiblast\n";
$cpunum=4;


#################################
##### report node -------->
`hostname`=~/(\S+)/;
$node=$1;
printf "hostname: $node\n";
$time=`date`;
printf "starting time: $time";
$pwd=`pwd`;
printf "pwd: $pwd";
#^^^^^^^^^^^^^^^^^^^^^^^^^^

%ts=(
     'GLY'=>'G',
     'ALA'=>'A',
     'VAL'=>'V',
     'LEU'=>'L',
     'ILE'=>'I',
     'SER'=>'S',
     'THR'=>'T',
     'CYS'=>'C',
     'MET'=>'M',
     'PRO'=>'P',
     'ASP'=>'D',
     'ASN'=>'N',
     'GLU'=>'E',
     'GLN'=>'Q',
     'LYS'=>'K',
     'ARG'=>'R',
     'HIS'=>'H',
     'PHE'=>'F',
     'TYR'=>'Y',
     'TRP'=>'W',

     'ASX'=>'B',
     'GLX'=>'Z',
     'UNK'=>'X',

     'G'=>'GLY',
     'A'=>'ALA',
     'V'=>'VAL',
     'L'=>'LEU',
     'I'=>'ILE',
     'S'=>'SER',
     'T'=>'THR',
     'C'=>'CYS',
     'M'=>'MET',
     'P'=>'PRO',
     'D'=>'ASP',
     'N'=>'ASN',
     'E'=>'GLU',
     'Q'=>'GLN',
     'K'=>'LYS',
     'R'=>'ARG',
     'H'=>'HIS',
     'F'=>'PHE',
     'Y'=>'TYR',
     'W'=>'TRP',

     'a'=>'CYS',
     'b'=>'CYS',
     'c'=>'CYS',
     'd'=>'CYS',
     'e'=>'CYS',
     'f'=>'CYS',
     'g'=>'CYS',
     'h'=>'CYS',
     'i'=>'CYS',
     'j'=>'CYS',
     'k'=>'CYS',
     'l'=>'CYS',
     'm'=>'CYS',
     'n'=>'CYS',
     'o'=>'CYS',
     'p'=>'CYS',
     'q'=>'CYS',
     'r'=>'CYS',
     's'=>'CYS',
     't'=>'CYS',
     'u'=>'CYS',
     'v'=>'CYS',
     'w'=>'CYS',
     'x'=>'CYS',
     'y'=>'CYS',
     'z'=>'CYS',

     'B'=>'ASX',
     'Z'=>'GLX',
     'X'=>'CYS',
    );

@AA=qw(
       C
       M
       F
       I
       L
       V
       W
       Y
       A
       G
       T
       S
       Q
       N
       E
       D
       H
       R
       K
       P
       );


open(seqtxt,"seq.fasta");
$sequence="";
while($line=<seqtxt>)
{
    goto pos1_1 if($line=~/^>/);
    if($line=~/(\S+)/)
    {
        $sequence .=$1;
    }
  pos1_1:;
}
close(seqtxt);
$Lch=length $sequence;
for($i=1;$i<=$Lch;$i++)
{
    $a=substr($sequence,$i-1,1);
    $seqQ{$i}=$a; 
}


`cp seq.aln deepmsa_protein.aln`;
`cp seq.a3m deepmsa_protein.a3m`;

#printf "$a3m2psiblast. deepmsa_protein.a3m psitmp -a3m -neff 7";
`$a3m2psiblast deepmsa_protein.a3m psitmp -a3m -neff 7`;   ## for ss


`cp psitmp.mtx  protein.mtx`;
`cp psitmp.pssm protein.pssm`;

open(msa,"deepmsa_protein.aln");
$it=0;
my $queryseq=<msa>;
$L=length $queryseq;
#printf "length::::: $L\n";
$it++;
while($line=<msa>)
{
    for($i=1;$i<=$L;$i++)
    {
        $am{$it,$i}=substr($line,$i-1,1);
    }
    $it++;
}
#printf "it::::: $it\n";

####### include query #######
for($i=1;$i<=$L;$i++)
{
    $am{$it,$i}=substr($queryseq,$i-1,1); #query residue
    #printf "$am{$it,$i}\n";
}

####### Henikoff weight $wei{i_seq} ----------->
##### nA{A,i_pos}: number of times A appear at the position:
undef %nA; 
for($i=1;$i<=$Lch;$i++){
    for($j=1;$j<=$it;$j++){
	$nA{$am{$j,$i},$i}++;
	#printf "$nA{$am{$j,$i},$i}\n";
    }
    #printf "$nA{$am{$j,$i},$i}\n";
}
##### henikoff weight w(i)=sum of 1/rs:
for($i=1;$i<=$it;$i++){
    for($j=1;$j<=$Lch;$j++){
	####### r: number of different residues in j'th position:
	$r=0;
	foreach $A(@AA){
	    $r++ if($nA{$A,$j}>0);
	}
	$A=$am{$i,$j};
	$s2=$nA{$A,$j};
	$w1=1.0/($r*$s2);
	$w{$i}+=$w1;
    }
    $w_all+=$w{$i};
}
#printf "w_all::::: $w_all\n";

#### normalization of w(i):
for($i=1;$i<=$it;$i++){
    $w{$i}/=$w_all;
}
#^^^^^ Henikoff weight finished ^^^^^^^^^^^^^^^^^

########### weighted frequence #################
undef %log;
for($i=1;$i<=$it;$i++){
    for($j=1;$j<=$Lch;$j++){
	$A=$am{$i,$j};
	$log{$j,$A}+=$w{$i};
    }
}
#^^^^^^^^^ Henikoff frequence finished ^^^^^^^^^^^^^

open(freq,">seq.prf");
printf freq "$Lch\n";
for($i=1;$i<=$Lch;$i++){
    printf freq "%3d $seqQ{$i} %3d",$i,$i;
    $norm=0;
    foreach $A(@AA){
        $norm+=$log{$i,$A};
    }
    foreach $A(@AA){
        printf freq "%10.7f",$log{$i,$A}/$norm;
    }
    printf freq "\n";
}
close(freq);



sub getRandomString
{
    my $len = shift;
    my @W = ('0' .. '9', 'a' .. 'z', 'A' .. 'Z');
    my $str;
    my $i = 0;
    while ($i++ < $len) {
        $str .= $W[rand(@W)];
    }
    return $str;
}

