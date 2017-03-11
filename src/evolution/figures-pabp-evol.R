# Scripts used to generate Fig. 4 and Fig. S4 in Riback/Katanski et al. Cell 2017
# http://drummondlab.org/papers/paper/adaptive-phase-separation
# Citation:
# 	Riback JA, Katanski CD, Kear-Scott JL, Pilipenko EV, Rojek AE, Sosnick TR, Drummond DA. Stress-triggered phase separation is an adaptive, evolutionarily tuned response. Cell 168(6):1028–1040 (2017).

# Relative paths used
# Set working directory to the pab1-phase-2017/src/evolution/ directory
# E.g. setwd("~/research/pab1-phase-2017/src/evolution/")


library(Hmisc)
library(ggplot2)
library(reshape2)
library(cowplot)
library(data.table)

file.out = F
output.type = 'svg'

data.dir <- '../../data'
figures.dir <- '../../data'
# github.com/dad/base
base.dir <- '../../../base'

pdir <- function(dirname, filename) {
	paste(dirname, filename, sep='/')
}

source(pdir(base.dir, "stat-lib.R"))

load.data = T
fig.pab1.aa.composition = T
fig.pab1.aa.hydrophobicity = T
fig.pabp.length.distribution = T
fig.all = T

stderr.nona <- function(x) {
	res <- sderr(x)
	if (is.na(res)) {
		res <- 0
	}
	res
}

all.aas = charlist('ACDEFGHIKLMNPQRSTVWY')
aas = as.character(all.aas)
# Colors
cols <- c(P="#000000ff", dP="#781c86ff", disprot="#56B4E9ff", all='#E69F00ff', dP.all="#781c86ff",P.all="#000000ff")


if (load.data) {
	xdis <- read.table(pdir(data.dir, "disprot-aa-freqs.txt"), header=T)
	xall <- read.table(pdir(data.dir, "scer-proteome-aa-freqs.txt"), header=T)
	xp.all <- read.table(pdir(data.dir, "pabp-pdomain-freqs.txt"), header=T)
	xdp.all <- read.table(pdir(data.dir, "pabp-non-pdomain-freqs.txt"), header=T)
	xp <- subset(xp.all, orf=='S.cerevisiae_S288c')
	xdp <- subset(xdp.all, orf=='S.cerevisiae_S288c')
	#xp <- subset(xp, 
	#xdp <- read.table(pdir(data.dir, "pab1dp-props.txt"), header=T)
}

# Fig. 4B in Riback/Katanski et al.
if (fig.pab1.aa.composition) {

	aas <- as.character(all.aas)
	f <- p.0("f",aas)
	d <- list("P"=xp[,f], "dP"=xdp[,f], "disprot"=xdis[,f], "all"=xall[,f])

	xs <- lapply(names(d), function(n){
		x <- as.data.frame(d[[n]]); 
		colnames(x) <- aas;
		# Compute summary stats
		x.summary <- data.frame(
			aa=aas,
			mean=apply(x, 2, mean),
			median=apply(x, 2, median),
			se=apply(x, 2, function(d){stderr.nona(d)*1.96}),
			sd=apply(x, 2, sd),
			supper=apply(x, 2, quantile, probs=0.975),
			slower=apply(x, 2, quantile, probs=0.025)
			)
		rownames(x.summary) <- aas;
		x.summary
		})
	names(xs) <- names(d)
	# Order by overall frequency first, then by enrichment
	ord1 <- order(xs[['all']][,'mean'])
	ord <- order(ratio <- xs[['P']][ord1,'mean']/xs[['dP']][ord1,'mean'], decreasing=T)
	# Now melt
	xmf <- lapply(names(xs), function(n){
		x.summary <- xs[[n]];
		x.summary$aa <- factor(x.summary$aa, levels=aas[ord1][ord])
		m <- melt(x.summary, id.vars=c("aa",'median','se','sd','supper','slower')); 
		m$name <- n;
		m})
	xm <- do.call(rbind, xmf)
	xm$name <- factor(xm$name, levels=c('dP','disprot','all','P'))

	lwd <- 0.5
	dodgewd <- 0.25
	p <- ggplot(data=subset(xm, name!='P'), aes(x=aa, y=value, group=name, colour=name, fill=name)) + 
		geom_line(size=lwd, position=position_dodge(width=dodgewd)) + 
		geom_pointrange(aes(ymin=value-se, ymax=value+se), size=lwd/3, position=position_dodge(width=dodgewd), alpha=1)  +
		coord_cartesian(ylim=c(0,0.2)) +
		scale_color_manual(values=cols) + xlab('Amino acid') + ylab('Proportion')
	p <- p + geom_pointrange(data=subset(xm, name=='P'), aes(ymin=value-se, ymax=value+se), size=lwd/3, 
		position=position_dodge(width=dodgewd)) + 
	geom_line(data=subset(xm, name=='P'), size=lwd, position=position_dodge(width=dodgewd))

	pg <- plot_grid(p)
	print(pg)
	if (file.out) save_plot(pdir(figures.dir,'fig_pab1_aa_proportions.svg'), pg, base_width=5, base_height=3)

	p2 <- ggplot(data=xm, aes(x=aa, y=value, group=name, colour=name, fill=name)) + 
		geom_line(size=lwd, position=position_dodge(width=dodgewd)) + 
		geom_pointrange(aes(ymin=slower, ymax=supper), size=lwd/3, position=position_dodge(width=dodgewd), alpha=1)  +
		coord_cartesian(ylim=c(0,0.2)) +
		scale_color_manual(values=cols) + xlab('Amino acid') + ylab('Proportion')
	p2g <- plot_grid(p2)
	print(p2g)
	if (file.out) save_plot(pdir(figures.dir,'fig_pab1_aa_proportions_variation.svg'), p2g, base_width=5, base_height=3)

	# Reviewer wants statistics.
	y <- data.frame(P=xs$P$mean,dP=xs$dP$mean,disprot=xs$disprot$mean,all=xs$all$mean)
	rc <- rcormat(y, meth='p')
	# These correlations are the ones reported in the text.

	sub.aas <- charlist("DEKRFYILVMA")
	f2 <- p.0('f', sub.aas)

	xs.s <- lapply(xs, function(x) {
		y <- x[sub.aas,]
		y$aa <- factor (y$aa, levels=sub.aas)
		y
	})
	xs.sm <- rbindlist(xs.s, idcol='name') #do.call(rbind, xs.s) #melt(xs.s, id=c("aa",'se','sd','supper','slower'))
	p3 <- ggplot(data=xs.sm, aes(x=aa, y=mean, group=name, colour=name, fill=name)) + 
		geom_line(size=lwd, position=position_dodge(width=dodgewd)) + 
		geom_pointrange(aes(ymin=mean-se, ymax=mean+se), size=lwd/3, position=position_dodge(width=dodgewd), alpha=1)  +
		coord_cartesian(ylim=c(0,0.12)) +
		scale_color_manual(values=cols) + xlab('Amino acid') + ylab('Proportion')
	p3g <- plot_grid(p3)
	print(p3g)
	if (file.out) save_plot(pdir(figures.dir,'fig_pab1_aa_proportions_subset.svg'), p3g, base_width=5, base_height=3)

}

# Fig. 4D in Riback/Katanski et al.
if (fig.pab1.aa.hydrophobicity) {
	# Filter out too-short sequences?

	hyd <- read.table(pdir(base.dir, "data/hydrophobicity-scales.txt"), header=T, sep='\t')
	rownames(hyd) <- hyd$aa
	hyd[,'Hopp.Z'] <- as.numeric(scale(-hyd[,'Hopp.Woods']))
	hyd[,'Hessa.phobic'] <- -hyd[,'Hessa']
	hyd.scale <- 'Hopp.Z'

	faa <- p.0("f",aas)
	f <- c(p.0("f",aas),p.0("n",aas))
	#d <- list("P"=xp[,f], "dP"=xdp[,f], "disprot"=xdis[,f], "all"=xall[,f], "P.all"=xp.all[,f], "dP.all"=xdp.all[,f])
	d <- list("P.all"=xp.all[,f], "dP.all"=xdp.all[,f], "disprot"=xdis[,f], "all"=xall[,f])
	target.aas <- charlist('ALIMV')
	xml <- lapply(names(d), function(n){
		x <- as.data.frame(d[[n]][,faa]); 
		colnames(x) <- aas;
		# Compute summary stats
		x.summary <- data.frame(
			aa=aas,
			hydrophobicity=as.numeric(hyd[aas,hyd.scale]),
			mean=apply(x, 2, mean),
			se=apply(x, 2, function(d){stderr.nona(d)*1.96}),
			sd=apply(x, 2, sd)
			)
		rownames(x.summary) <- aas;
		#print(x.summary)
		m <- melt(x.summary[target.aas,], id.vars=c("aa",'hydrophobicity','se','sd')); 
		m$name <- n; 
		m})
	xm <- do.call(rbind, xml)
	# Correlations between target amino acid frequency and hydrophobicity
	corf <- function(x,hyds, pval=FALSE) {
		pseudocount <- 1
		ct <- cor.test(hyds, log(x+pseudocount), method='p', exact=FALSE)
		if (pval) return(ct$p.value)
		return(ct$estimate)
	}
	xml.cor <- lapply(names(d), function(n){
		x <- as.data.frame(d[[n]])[,p.0('n',target.aas)];
		colnames(x) <- target.aas;
		cors <- apply(x, 1, corf, hyds=as.numeric(hyd[target.aas,hyd.scale]))
		cors})
	names(xml.cor) <- names(d)
	
	xml.cor.p <- lapply(names(d), function(n){
		x <- as.data.frame(d[[n]])[,p.0('n',target.aas)];
		colnames(x) <- target.aas;
		cors <- apply(x, 1, corf, pval=TRUE, hyds=as.numeric(hyd[target.aas,hyd.scale]))
		cors})
	names(xml.cor.p) <- names(d)
	print(sapply(xml.cor, median))
	print(sapply(xml.cor.p, median))
	#xml.lm <- lapply(names(d), function(n){
	#	x <- as.data.frame(d[[n]])[,p.0('f',target.aas)];
	#	colnames(x) <- target.aas;
	#	slopes <- apply(x, 1, function(x,y){coef(lm(log(y)~x))[2]}, y=as.numeric(hyd[target.aas,hyd.scale]))
	#	slopes})
	#names(xml.lm) <- names(d)
	lwd <- 0.5
	gd <- ggplot(melt(xml.cor), aes(x=value, colour=L1)) + stat_ecdf(size=lwd) + 
		scale_color_manual(values=cols) + 
		geom_vline(xintercept=0, linetype='dotted') + 
		xlab(bquote('Frequency-hydrophobicity correlation ('*R*')')) +
		ylab('Cumulative distribution')
	gd2 <- ggplot(melt(xml.cor), aes(x=value^2, colour=L1)) + stat_ecdf(size=lwd) + 
		scale_color_manual(values=cols) + 
		xlab(bquote('Variance in frequency explained by hydrophobicity ('*R^2*')')) +
		ylab('Cumulative distribution')
	varplot <- plot_grid(gd+no.legend, gd2+no.legend)
	#print(varplot)
	if (file.out) save_plot(pdir(figures.dir,'fig_pab1_aa_hydrophobicity_corrs.svg'), varplot, base_height=3, base_aspect_ratio=2.2)
	#print(wilcox.test(xml.cor$P.all, xml.cor$dP.all, alternative='less'))
	#print(wilcox.test(xml.cor$P.all, xml.cor$disprot, alternative='less'))
	#print(wilcox.test(xml.cor$P.all, xml.cor$all, alternative='less'))
	#xm <- do.call(rbind, xml)

	p <- ggplot(data=xm, aes(x=hydrophobicity, y=value, label=aa, group=name, colour=name, fill=name)) + 
		scale_y_log10() + annotation_logticks(sides="l") +
		geom_line(size=lwd) + geom_point(size=0.6) + geom_pointrange(aes(ymin=value-se, ymax=value+se), size=lwd/2) +
		scale_color_manual(values=cols) + xlab('Hydrophobicity') + ylab('Proportion')
	p <- p + geom_text(data=subset(xm, name=='P.all'), stat='identity', size=5, nudge_x=0.03, nudge_y=0.05)
	pg <- plot_grid(p+no.legend)
	#print(pg)
	if (file.out) save_plot(pdir(figures.dir,'fig_pab1_aa_hydrophobicity.svg'), pg, base_height=3, base_aspect_ratio=1)

	# Reviewer 4 wants statistics.
	# What are those statistics?
	# 1) 
	stat.res <- lapply(xml.cor, function(x){wilcox.test(x, alternative='less')$p.value})
	names(stat.res) <- names(xml.cor)
	#print(stat.res)
	#print(wilcox.test(xml.cor$P.all, xml.cor$dP.all, alternative='less'))

	cor.thresh <- -0.25
	for (n in names(xml.cor)) {
		y <- xml.cor[[n]]
		cat("# There are", length(y[y< cor.thresh])/length(y), "entries in", n, "less than", cor.thresh, "\n")
		cat("# Median", median(y), "\n")
	}
}

if (fig.pabp.length.distribution) {
	xm <- melt(d <- list(P.all=xp.all, dP.all=xdp.all, all=xall, disprot=xdis), id='orf', meas='f.P')

	lwd <- 0.5
	gL <- ggplot(xp.all, aes(x=length)) + stat_ecdf(size=lwd) +
		xlab('P domain length') + ylab('Cumulative distribution')
	gP <- ggplot(xm, aes(x=value, colour=L1)) + stat_ecdf(size=lwd) +
		scale_color_manual(values=cols) +
		xlab('Proline content (fraction)') + ylab('Cumulative distribution')
	
	g.length.dist <- plot_grid(gL, gP+no.legend)
	print(g.length.dist)
	if (file.out) save_plot(pdir(figures.dir,'fig_pab1p_length_distribution.svg'), g.length.dist, base_height=3, base_aspect_ratio=2.2)

	print(wilcox.test(d[['P.all']]$f.P, d[['all']]$f.P, alternative='greater'))
	print(wilcox.test(d[['P.all']]$f.P, d[['dP.all']]$f.P, alternative='greater'))
	print(wilcox.test(d[['P.all']]$f.P, d[['disprot']]$f.P, alternative='greater'))

	L.thresh <- 20
	P.thresh <- 0.6
	n.small.p <- nrow(nss <- subset(d[['P.all']], length < L.thresh & f.P>P.thresh))
	cat("# All but", n.small.p, "contain a proline-rich LCR of at least", L.thresh, "amino acids\n")

	y <- data.frame(P=d[['P.all']][,c('length','f.P')], dP=d[['dP.all']][,c('length','f.P')], diff=d[['P.all']]$f.P-d[['dP.all']]$f.P)
	cat("# There are", nrow(subset(y, diff<0 | P.length<L.thresh)), "PABP LCR regions shorter than", L.thresh, "or without proline enrichment relative to the rest of PABP\n")
}


# Fig. S4 in Riback/Katanski et al.
if (fig.all) {
	pg <- plot_grid(
		plot_grid(p3, labels='A'),
		plot_grid(gL+no.legend, gP+no.legend, labels=c("B","C")),
		plot_grid(gd+no.legend, gd2+no.legend, labels=c("D","E")),
		nrow=3
		)
	print(pg)
	if (file.out) save_plot(pdir(figures.dir,'fig_pab1_evol_suppl.svg'), pg, base_width=7.5, base_height=7.5*(11/8.5))

}