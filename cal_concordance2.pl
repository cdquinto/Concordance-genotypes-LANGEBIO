#!/usr/bin/perl
use warnings;
use File::Slurp;

@array=read_file("/home/cdquinto/test_langebio/analysis_entrenamiento_clusterGLOBAL/merged_data_hapmap_3_ok.ped");

open(OUT,">results.txt");
open(OUT2,">resultsR.txt");

#@matrix = ([0,6],[1,4],[2,5],[3,7] ); 
@matrix=([0,7],[1,5],[2,8],[3,6],[4,9]);

#for($i=0;$i<4;$i++){
for($i=0;$i<5;$i++){
	chomp $array[$matrix[$i][0]];
	chomp $array[$matrix[$i][1]];
	@split1=split(/\-9/,$array[$matrix[$i][0]]);
	$split1[1] =~ s/ //g;
	@string1=split(//,$split1[1]);

	@split2=split(/\-9/,$array[$matrix[$i][1]]);
	@name=split(/ /,$split2[0]);
	$split2[1] =~ s/ //g;
	@string2=split(//,$split2[1]);

	$cont=0;
	
	foreach $i (scalar($string1)) {
  		if($string1[$i] eq $string2[$i]){
  			$cont++;
  		}
  		elsif ($string1[$i] eq 'G' && $string2[$i] eq 'C'){
  			$cont++;
  		}
  		elsif ($string1[$i] eq 'C' && $string2[$i] eq 'G'){
  			$cont++;
  		}
  		elsif ($string1[$i] eq 'T' && $string2[$i] eq 'A'){
  			$cont++;
  		}
  		elsif ($string1[$i] eq 'A' && $string2[$i] eq 'T'){
  			$cont++;
  		}
	}
	
	$long=length($split1[1]);
	$ratio=$cont/$long;
	$minus_ratio=1-$ratio;
	#print OUT "$name[0]\t$matches\t$long\t$ratio\t$minus_ratio\n";
	#print OUT2 "$name[0]\t$ratio\t$minus_ratio\n";
	print OUT "$name[1]\t$matches\t$long\t$ratio\t$minus_ratio\n";
	print OUT2 "$name[1]\t$ratio\t$minus_ratio\n";
}

die;

open(RCMD,">commandsR.txt");
print RCMD "data <- read.delim(\"resultsR.txt\", header=FALSE, row.names=1)\n";
#print RCMD "pdf(\"Concordance_1KG_MEGA_clusterGLOBAL.pdf\")\n";
print RCMD "pdf(\"Concordance_HAPMAP_MEGA_clusterGLOBAL.pdf\")\n";
print RCMD "at_tick <- seq(0.1,0.9,0.2)\n";
#print RCMD "barplot(t(as.matrix(data)), ylim=c(0,1), col=c(\"chocolate1\",\"cadetblue1\"), main=\"Concordance between 1KG and MEGA\")\n";
print RCMD "barplot(t(as.matrix(data)), ylim=c(0,1), col=c(\"chocolate1\",\"cadetblue1\"), main=\"Concordance between HAPMAP3 and MEGA\")\n";
print RCMD "abline(h = mean(data\$V2), col = \"black\", lwd=5)\n";
print RCMD "axis(2, at = at_tick,tcl = -0.25)\n";
print RCMD "dev.off()\n";
close(RCMD);

system ( "R --slave < commandsR.txt" );
