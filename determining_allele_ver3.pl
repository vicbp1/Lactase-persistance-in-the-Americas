#!usr/bin/perl
use Getopt::Long;
#This script creates an vcf file using inputs and outputs of Rfmix and also introduce NA for ancestries non required (masking)
#We consider that the first ancestry class of the rfmix classes file was the Natives: 1: Natives, 2: African 3:European

use Data::Dumper;
GetOptions (
	"vcf=s"=>\$vcf,
	"variantD=s"=>\$variantD,
	"variantL=s"=>\$variantL,
	"identity=s"=>\$identity,
	"forback=s"=>\$forback,
	"output=s"=>\$output,
	"help!"=>\$help,
)or die(showhelp());

if($vcf eq "" or $identity eq ""  or $variantD eq "" or $variantL eq "" or $forback eq "" or $output eq ""){	#Se o usuario pedir ajuda ou algum dos arquivos nao forem setados
	showhelp();							#apresente ajuda
}

open (IF,"$vcf") or die ("We cannot find the vcf file\n");

$homozygotesancestral=0;
$hap1=0;
$hap2=0;
$hap3=0;
while (<IF>){
	chomp;
	if (m/^#CHROM/){
		@line = split /\s+/;
		foreach $word (0..$#line){
			if (@line[$word]=~ m/$identity/){  ### delicate, depends on pattern of ind names
				push (@inds, @line[$word]);
				push (@positions, $word);
					
			}
		}
		$nind = scalar @inds;
	}

	if (m/$variantD/){
		@genotype = split /\s+/;
		foreach $i(0..$#positions){
			if (@genotype[@positions[$i]]=~ m/0\|0/){  ### to identify hoimozygotes for the ancestral
				$homozygotesancestral++;
					
			}elsif (@genotype[@positions[$i]]=~ m/1\|0/){
				$hash_genotype_hap1{@inds[$i]}="hap1";
				$hap1++;
	
			}elsif (@genotype[@positions[$i]]=~ m/0\|1/){
				$hash_genotype_hap2{@inds[$i]}="hap2";
				$hap2++;

			}elsif (@genotype[@positions[$i]]=~ m/1\|1/){		
				$hash_genotype_twohaps{@inds[$i]}="hap1hap2";				
				$hap3++;
			}
		}

	}
}

$derived=$hap1+$hap2+($hap3*2);
print "$homozygotesancestral of $nind are Homozygotes for Ancestral allele and $derived derived alleles\n";


$europeanderived=0;
$noneuropeanderived=0;
open (IF2,"$forback") or die ("We cannot find the forward backward file\n");

while (<IF2>){
	chomp;

	if (m/^chromosome/){
		@fbheader = split /\s+/;
		foreach $j (3..$#fbheader){			### Considering that the first three columns are chr-physical-genetic position info	
			if (@fbheader[$j]=~ m/$identity/){  ### delicate, depends on pattern of ind names
	
				@header = split (/:::/,@fbheader[$j]);

				
				if(exists($hash_genotype_hap1{@header[0]})){				
					
					push (@positionsinds_1, $j);		#### just if it exist in the hash with derive allele
					push (@indslocal_1, @header[0]);
					push (@haplocal_1, @header[1]);
					push (@ancestrylocal_1, @header[2]);
				
				}

				if(exists($hash_genotype_hap2{@header[0]})){				
					
					push (@positionsinds_2, $j);		#### just if it exist in the hash with derive allele
					push (@indslocal_2, @header[0]);
					push (@haplocal_2, @header[1]);
					push (@ancestrylocal_2, @header[2]);
				
				}

				if(exists($hash_genotype_twohaps{@header[0]})){				
					
					push (@positionsinds_3, $j);		#### just if it exist in the hash with derive allele
					push (@indslocal_3, @header[0]);
					push (@haplocal_3, @header[1]);
					push (@ancestrylocal_3, @header[2]);
				
				}
			}
		}
	}







	if (m/$variantL/){
		@fbprobs = split /\s+/;
		foreach $k (0..$#positionsinds_1){			### Array created with the information of the header of derived allele
			$place1 = @positionsinds_1[$k];
			if (@fbprobs[$place1] > 0.9499){
				if (@haplocal_1[$k] eq "hap1"){
					if (@ancestrylocal_1[$k] eq "EUR"){
										$europeanderived++;
					}else{
										$noneuropeanderived++;
					}
				}
			}
		}


		foreach $l (0..$#positionsinds_2){			### Array created with the information of the header of derived allele
			$place2 = @positionsinds_2[$l];
			if (@fbprobs[$place2] > 0.9499){
				if (@haplocal_2[$l] eq "hap2"){
					if (@ancestrylocal_2[$l] eq "EUR"){
										$europeanderived++;
					}else{
										$noneuropeanderived++;
					}
				}
			}
		}

		foreach $l (0..$#positionsinds_3){			### Array created with the information of the header of derived allele
			$place3 = @positionsinds_3[$l];
			if (@fbprobs[$place3] > 0.9499){
					if (@ancestrylocal_3[$l] eq "EUR"){
										$europeanderived++;
					}else{
										$noneuropeanderived++;
			
					}
			}
		}


	}


}


print "$europeanderived\n$noneuropeanderived\n";

