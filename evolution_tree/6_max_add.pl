## as the blackground of class simialrity
open(IN,"Total_dup_species.Cor.ann_subcluster.sort_1.filter.txt");
open(IN1,"Total_dup_species.subclass_1.Cor.txt");
<IN1>;
while(<IN1>){
	chomp;
	my @row=split/\t/;
	if($row[2]>0.8){
		$flag{$row[0]}{$row[1]}=$row[2];
		$ff{$row[0]}=$_;
		}
	}
close(IN1);	

open(IN1,"Total_dup_species.subclass_1.Cor.txt");
<IN1>;
while(<IN1>){
	chomp;
	my @row=split/\t/;
	#if($row[2]<0.6 || defined $ff{$row[0]}){next;}
	#else{
	if(defined $flag2{$row[0]}){
			my @temp=split/\t/,$flag2{$row[0]};
			
			if($temp[2]>=$row[2]){next;}else{$flag2{$row[0]}=$_;}
		}else{
				$flag2{$row[0]}=$_;	
		}
	
	#}
}

close(IN1);
my $cnt = 1;
open(OUU,">TT_08_1-1.out");
foreach my $k(keys %flag){
	foreach my $kk(keys %{$flag{$k}}){
	print OUU "$k\t$kk\t$flag{$k}{$kk}\n";
	$cnt++;
	}
}
print $cnt."\n";


my $cnt = 1;
open(OUTT,">TT_08-1.out");
foreach my $k(keys %flag2){
	my @temp=split/\t/,$flag2{$k};
	$flag3{$temp[0]}{$temp[1]}=$temp[2];
	print OUTT $flag2{$k}."\n";
	$cnt++;
	}
print $cnt."\n";




my $head=<IN>;	
my $cnt = 1;
while(<IN>){
	chomp;
	my @row=split/\t/;
	# max or >0.8
	if(defined $flag{$row[2]}{$row[5]} || defined $flag3{$row[2]}{$row[5]}){
	$cnt++;
	#if( defined $flag2{$row[2]}{$row[5]}){
		# save max
		if(defined $hash{$row[0]}){
			my @temp=split/\t/,$hash{$row[0]};
			if($temp[6]>=$row[6]){next;}
			else{$hash{$row[0]}=$_;}
		}else{
				$hash{$row[0]}=$_;	
		}
	}
}
open(OUT,">Total_dup_species.Cor.ann.sort.max_8_subclass_1.txt");
print OUT $head;
print $cnt."\n";



foreach my $k(keys %hash){
	
		print OUT $hash{$k}."\n";
		
	
	}
	
	
	
print OUT "Cnidocytean.n	Nematostella	Cnidocytean.n	root	Root	root	1
Dig_filamentsan.n	Nematostella	DigFilaments.n	root	Root	root	1
Epitheliuman.n	Nematostella	Epithelial.n	root	Root	root	1
Gastrodermisan.n	Nematostella	Gastrodermis.n	root	Root	root	1
Gland_secretoryan.n	Nematostella	Secretory.n	root	Root	root	1
Musclean.n	Nematostella	Muscle.n	root	Root	root	1
Neuronan.n	Nematostella	Neuron.n	root	Root	root	1
Precursorsan.n	Nematostella	Proliferating.n	root	Root	root	1"

