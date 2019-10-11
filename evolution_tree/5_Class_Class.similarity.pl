open(IN,"Total_dup_species.Cor.ann_subcluster.sort_1.txt");
my $head=<IN>;
while(<IN>){
	chomp;
	my @row=split/\t/;
	$result{$_}++;
	}

open(OUT,">Total_dup_species.Cor.ann_subcluster.sort_1.filter.txt");
print OUT $head;
foreach my $k(keys %result){print OUT "$k\n";}
close(IN);
close(OUT);



my %hash;
my %count;
open(IN,"Total_dup_species.Cor.ann_subcluster.sort_1.filter.txt");
while(<IN>){
	chomp;
	my @row=split/\t/;
	
	
	#if($row[6]>0.1){
		if(defined $hash{$row[2]}{$row[5]}){
			$hash{$row[2]}{$row[5]}+=$row[6];
			$count{$row[2]}{$row[5]}++;
		}else{
				$hash{$row[2]}{$row[5]}=$row[6];
				$count{$row[2]}{$row[5]}++;	
	}
	#}

}

open(OUT,">Total_dup_species.subclass_1.Cor.txt");
print OUT "Cluster1\tCluster2\tAUROC\n";
foreach my $k(keys %hash){
	foreach my $kk(keys %{$hash{$k}}){
		if(defined $count{$k}{$kk}){
			$oo=$hash{$k}{$kk}/$count{$k}{$kk};
			print OUT "$k\t$kk\t$oo\n";
			
			}
		
		}
	
	}
close(IN);
close(OUT);
