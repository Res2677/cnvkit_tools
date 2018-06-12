use Getopt::Long;
GetOptions(
	"g=s" => \$gene,
	"rf=s" => \$rf,
	"cnr=s"  => \$cnr,
	"pcolor=s" => \$pcolor,
	"o=s" => \$outdir,
	"h|?" => \$help
);

if (!defined $gene || !defined $rf || !defined $cnr || !defined $pcolor || !defined $outdir || defined $help) {
        die << "USAGE";
        description: CNV plot
        usage: perl $0 [options]
options:
        -g	gene symbol
	-rf	file like ALL.gene.exon_intron.txt
	-cnr	file of cnvkit output,like ZF001-C_C.cnr
	-pcolor	color of CNV point
	-o	outdir
	-help|?	      information for help
e.g.:
perl $0 -g EGFR -rf /public/home/hangjf/data/database/hg38/ALL.gene.exon_intron.txt -cnr /public/home/hangjf/cancer-changhai/cnvkit1/data/ZF001_C/ZF001-C_C.cnr -pcolor "#E41A1C" -o /public/home/hangjf/cancer-changhai/cnvkit1/data/ZF001_C
USAGE
}

$trans = `grep '$gene'\$'\\t' $rf`;
$scatter = `grep -A10 -B10 '$gene' $cnr`;
open OUT1,">$outdir/$gene.trans.txt",or die $!;
print OUT1 "$trans";
close OUT1;
open OUT2,">$outdir/$gene.scatter.txt",or die $!;
print OUT2 "$scatter";
close OUT2;

open my $fh_rcode,">","$outdir/plot_CNV.R" or die $!;
print $fh_rcode <<CMD;

library(ggplot2)
library(RColorBrewer)
cols<-brewer.pal(n=4,name="Set1")

pdf("$outdir/$gene.CNV.pdf",8,8)
tt <- data.frame(read.table('$outdir/$gene.trans.txt'))
scatter <- data.frame(read.table('$outdir/$gene.scatter.txt'))
interval <- max(scatter\$V2) -min(scatter\$V2)
leftlab <- min(scatter\$V2)+(50000 - (min(scatter\$V2)%%50000))
rightlab <- max(scatter\$V2) - (max(scatter\$V2)%%50000)
xlabat <- seq(leftlab,rightlab,50000)
xlab <- paste(seq(leftlab/1000000,rightlab/1000000,0.05),"M")
ylab <- seq(-10,10,5)
xmin <- (max(scatter\$V2) + min(scatter\$V2))/2

plot(scatter\$V2,scatter\$V5,ylab="",xlab="",cex.lab=1.5,yaxt="n",xaxt="n",axes=F,main="",col = alpha("$pcolor", 0.7), cex =1.3,pch=19,ylim = c(-14,13),xlim = c(min(scatter\$V2)-interval/10,max(scatter\$V2)+(interval/10)*4))
k <- seq(1,length(tt[tt\$V4=="intron",]\$V2),1)
for (i in k)
{
	segments(tt[tt\$V4=="intron",]\$V2[i],-11,tt[tt\$V4=="intron",]\$V3[i],-11,lwd=2,col="lightblue3")
}
k <- seq(1,length(tt[tt\$V4=="exon",]\$V2),1)
for (i in k)
{
	segments(tt[tt\$V4=="exon",]\$V2[i],-11,tt[tt\$V4=="exon",]\$V3[i],-11,lwd=5,col="steelblue")
}

segments(min(scatter\$V2)-interval/10,-13,min(scatter\$V2)-interval/10,10,lwd=2,col="grey40")
segments(max(scatter\$V2)+interval/10,-13,max(scatter\$V2)+interval/10,10,lwd=2,col="grey40")
segments(min(scatter\$V2)-interval/10,-10,max(scatter\$V2)+interval/10,-10,lwd=2,col="grey40")
segments(min(scatter\$V2)-interval/10,-13,max(scatter\$V2)+interval/10,-13,lwd=2,col="grey40")
segments(min(scatter\$V2)-interval/10,0,max(scatter\$V2)+interval/10,0,lwd=1.1,col="grey70",lty=2)
axis(1,col="grey40",col.axis="grey40",pos=-13,lwd=2,at=xlabat,labels=xlab)
axis(2,col="grey40",col.axis="grey40",pos=min(scatter\$V2)-interval/10,lwd=2,at=ylab,labels=ylab)
text(xmin,-12,labels="$gene")
par(new=T)
par(mar=c(6,6,6,10),xpd=TRUE)
points(max(scatter\$V2)+(interval/10)*2,0,col = alpha(cols[4], 0.7), cex =1.3,pch=19)
segments(max(scatter\$V2)+(interval/10)*2,-11,max(scatter\$V2)+(interval/10)*2.5,-11,lwd=2,col="lightblue3")
segments(max(scatter\$V2)+(interval/10)*2,-12,max(scatter\$V2)+(interval/10)*2.5,-12,lwd=5,col="steelblue")
text(max(scatter\$V2)+(interval/10)*3.5,0,labels="log2(CNV)")
text(max(scatter\$V2)+(interval/10)*3.5,-11,labels="intron")
text(max(scatter\$V2)+(interval/10)*3.5,-12,labels="exon")
mtext(side=2,"values of log2(CNV)",line=3,cex=1.3)
mtext(side=3,paste("CNV of $gene (",scatter[scatter\$V4=="$gene",][1,1],":",min(tt\$V2),"-",max(tt\$V3),")",sep=""),line=1,cex=1.3)
mtext(side=1,"Position in chromosome",line=2,cex=1.3)
dev.off()
CMD

#`Rscript $outdir/plot_CNV.R`;
system("Rscript $outdir/plot_CNV.R 2>/dev/null");
