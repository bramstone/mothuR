#Function to load in distance matrices created in mothur
#x = The name (and file path if necessary) of the .dist file to read in
read.mothurdist<-function(x) {
  groups<-scan(x,what=integer(),nlines=1,quiet=T)
  y<-data.matrix(read.table(x,sep='\t',fill=T,row.names=1,skip=1,col.names=1:groups))
  y<-cbind(y,rep(NA,dim(y)[1]))
  colnames(y)<-rownames(y)
  y<-as.dist(y,diag=FALSE,upper=FALSE)
  return(y)
}

#reads in shared file (specified by filename and file path if necessary) as an integer matrix (more memory efficient than data frame)
read.shared<-function(x) {
  y<-scan(x,what=character(),skip=1,sep='\t',quiet=T)
  y.names<-scan(x,what=character(),nlines=1,sep='\t',quiet=T)
  y<-matrix(y,ncol=length(y.names),byrow=T)
  colnames(y)<-y.names;rm(y.names)
  rownames(y)<-y[,'Group']
  otu_label<-y[1,'label']
  drop<-c('label','Group','numOtus')
  y<-y[,!colnames(y) %in% drop]
  storage.mode(y)<-"integer"
  colnames(y)<-sub('Otu(0)*','OTU ',colnames(y))
  attr(y, 'label') <- otu_label
  attr(y, 'filename') <- gsub('(.*)/(.*).shared', '\\2.shared', x)
  return(y)
}

#writes out shared file, either with user-specified filename or by appending 'editR' to the long shared file name
write.shared <- function(x, filename='', file.append='editR') {
  numOtus <- ncol(x)
  Group <- rownames(x)
  label <- attr(x, 'label')
  otu_names <- as.numeric(sub('OTU ([0-9*])', '\\1', colnames(x)))
  max_otu <- nchar(max(otu_names))
  otu_names <- formatC(otu_names, width=max_otu, format='d', flag='0')
  otu_names <- paste0('Otu', otu_names)
  if(filename=='') {
    append <- paste0('.', file.append, '.shared')
    filename <- paste0(sub('(.*).shared', '\\1', attr(x, 'filename')), append)
  }
  y <- x
  colnames(y) <- otu_names
  y <- data.frame(label, Group, numOtus, y)
  write.table(y, file=filename,
              quote=F, 
              row.names=F,
              sep='\t') #may need to tweek 'eol' (\n or \r\n for newline)
}
	
#Function to create ordination plot, useful for putting in for loops
plot.ord<-function(axisscores,grouping.1,grouping.2,metric,legend=FALSE,legend.names=c(),xlimit=c(-1,1),ylimit=c(-1,1),axes=c(1,2)) {
  plot.points<-c(21,24,23,22,25)
  if(missing(grouping.1)) {plot.colors<-c('lightskyblue4'); group.1.proxy<-1} else {plot.colors<-rainbow(length(unique(grouping.1))); group.1.proxy<-grouping.1}
  if(missing(grouping.2)) {plot.points<-c(21); group.2.proxy<-1} else {plot.points<-plot.points[1:length(unique(grouping.2))]; group.2.proxy<-grouping.2}
  if(missing(metric)) metric<-''
  #plot
  plot(axisscores[,axes[1]]~axisscores[,axes[2]],data=axisscores,bty='n',pch=plot.points[group.2.proxy],
    bg=plot.colors[group.1.proxy],xlab=paste(metric,'axis',axes[1]),ylab=paste(metric,'axis',axes[2]),
    ylim=ylimit,xlim=xlimit,xaxt='n',yaxt='n',bty='o',cex=1.25)
  #axes
  axis(1,tick=FALSE,at=seq(-1,1,by=0.4),line=-1);axis(2,tick=FALSE,at=seq(-1,1,by=0.4),line=-1)
  #legend
  if(legend==TRUE) {
    if(missing(grouping.1)) {
      warning('No grouping supplied to create legend')
    }
    else if(missing(grouping.2)) {
      L<-legend('topleft',legend=unique(grouping.1),title=legend.names[1],pt.cex=1.5,
        pt.bg=plot.colors,bty='n',xjust=0,x.intersp=0.75,y.intersp=1,pch=plot.points[1],
        text.font=ifelse(legend.names[1]=='species'|legend.names[1]=='species'|legend.names[1]=='host',3,1))
    }
    else {	
    #1st legend
      L<-legend('topleft',legend=unique(grouping.1),title=legend.names[1],pt.cex=1.5,
        pt.bg=plot.colors,bty='n',xjust=0,x.intersp=0.75,y.intersp=1,pch=plot.points[1],
        text.font=ifelse(legend.names[1]=='species'|legend.names[1]=='species'|legend.names[1]=='host',3,1))
    #2nd legend
    legend('bottomleft',legend=unique(grouping.2),
      title=legend.names[2],pch=plot.points,bty='n',
      xjust=0,x.intersp=0.75,y.intersp=1,pt.bg='black',
      text.font=ifelse(legend.names[1]=='species'|legend.names[1]=='species'|legend.names[1]=='host',3,1))
    }
  }
}
	

#Function which puts mothur's taxonomy and cons.taxonomy into an R data frame
taxonomy_to_data_frame<-function(x,confidence=FALSE,keep.counts=TRUE) {
  y<-scan(x,what=character(),sep='\t',skip=1,quiet=T)
  y.names<-scan(y,what=character(),nlines=1,quiet=T)
  colnames(y)<-y.names
  #remove heirarchy tags
  z<-gsub(pattern='.__','',as.character(y$Taxonomy))
  y<-y[,-3]
  #remove boostrapped confidence values
  if(confidence==FALSE) {
    z<-gsub(pattern='\\([^\\)]+\\)','',z)
  }
  #taxonomy is now only separated by ';', split by this into a data frame
  d<-scan(text=z,what=character(),sep=';',quiet=T,fill=TRUE)
  d<-d[,-8]
  names(d)<-c('Kingdom','Phyla','Class','Order','Family','Genus','Species')
  if(keep.counts==FALSE) {
    d<-cbind('OTU'=y[,'OTU'],d)
  }
  else {
    d<-cbind('OTU'=x[,'OTU'],'Size'=integer(x[,'Size']),d)
  }
  return(d)
}

#Function which summarizes mothur's tax.summary file by phyla
#x = filename (and path, if necessary) specifying a tax.summary file to load
#rel.abund = whether to calculate relative or absolute abundances
tax_summary_by_phyla<-function(x,rel.abund=TRUE) {
  #extract only phyla info
  x<-x[x$taxlevel==2,]
  #drop unneccesary columns
  drops<-c('taxlevel','rankID','daughterlevels','total')
  x<-x[,!names(x) %in% drops]
  #generate phyla totals across rows
  names(x)[1]<-'taxon'
  x$total<-rowSums(x[,-1])
  #scale values to relative abundances
  if(rel.abund==TRUE) {
    sample.sums<-colSums(x[,-1])
    for(i in 2:length(x)) {
        x[,i]<-x[,i]/sample.sums[i-1]}
  }
  #sort
  x<-x[rev(order(x$total)),]
  x$taxon<-gsub(pattern='.__','',as.character(x$taxon))
  x$taxon<-as.factor(x$taxon)
  return(x)
}

#Groups the relative percent of all phyla to be grouped as 'Other' together 
#for plotting
truncate_phyla<-function(x,cutoff=5,keep.unclassified=TRUE) {
	if(keep.unclassified==TRUE) {
		sum.y<-apply(x[(cutoff+1):nrow(x),-1],2,sum)
		x<-x[-c((cutoff+1):nrow(x)),]
		y<-data.frame('taxon'='Other',t(sum.y))
		names(y)<-names(x)
		x<-rbind(x,y)
		x$taxon<-factor(x$taxon,levels=as.ordered(x$taxon))
		return(x)}
	else{
		unclass<-x[which(x$taxon=='unclassified'|x$taxon=='Bacteria_unclassified'),]
		sum.y<-apply(x[(cutoff+1):nrow(x),-1],2,sum)
		x<-x[-c((cutoff+1):nrow(x)),]
		y<-data.frame('taxon'='Other',t(sum.y))
		names(y)<-names(x)
		y[,-1]<-unclass[,-1]+y[,-1]
		x<-x[-which(x$taxon=='unclassified'|x$taxon=='Bacteria_unclassified'),]
		x<-rbind(x,y)
		x$taxon<-factor(x$taxon,levels=as.ordered(x$taxon))
		return(x)}
}

#function for plotting sequences of dominant phyla from a truncated phyla table
#groups should be supplied in a list, and can come from either the sample names, or from an outside factor (recommended)
#group spacing should be supplied as a a logical vector specifying if the graph should be spaced by any particular factor (TRUE/FALSE for each factor)
plot_phyla<-function(x,type=c('both','summary','sample'),metric=c('total','mean'),groups=list(),grouping.in.names=FALSE,group.spacing=c(),label.angle=90,label.interval=1) {
	require('ggplot2')
	require('reshape2')
	require('gridExtra')
	x.colors<-rainbow(nrow(x)-1)
	col.seq<-seq(from=1,to=length(x.colors),by=2)
	opp.seq<-seq(from=2,to=length(x.colors),by=2);opp.seq<-rev(opp.seq)
	col.order<-paste(col.seq,opp.seq)
	col.order<-as.numeric(unlist(strsplit(col.order,' ')))
	x.colors<-x.colors[col.order]
	x.colors[nrow(x)]<-gray(0.8)
	type<-match.arg(type)
	if(type=='both'|type=='sample') {
		#If plotting by sample, the counting data to be plotted needs to be all in the same column
		#so it must be melted using the reshape2 package
		tp<-data.frame(t(x[,-1]))
		names(tp)<-as.character(x$taxon)
		tp<-cbind(data.frame('sample'=rownames(tp)),tp)
		tp$sample<-as.factor(tp$sample)
		total.row<-tp[tp$sample=='total',]
		tp<-tp[-nrow(tp),]
		rownames(tp)<-1:nrow(tp)
		composite.grouping.factor<-rep('',nrow(tp))
		if(length(groups)>0) {
			#To make groupings with sample-based bargraph, the tp object is edited prior to melting
			factor.list<-data.frame('sample'=droplevels(tp$sample)) #droplevels removes the level "total" from the sample column
			if(grouping.in.names==TRUE) {
				sample.names<-colnames(x)[2:(length(colnames(x))-1)]
				factor.list$sample<-as.character(sample.names)
				#create a data frame of factors from sample symbols
				for(i in 1:length(groups)) {
					var<-substr(factor.list$sample,groups[[i]][1],groups[[i]][2])
					factor.list<-data.frame(cbind(factor.list,as.factor(var)))
					names(factor.list)[i+1]<-paste0('var',i)
				}
			}
			else { #create a data frame of factors from input variables
				for(i in 1:length(groups)) {
					names(groups)<-paste0('var',i)
				}
				for(i in 1:length(groups)) {
					factor.list<-data.frame(cbind(factor.list,groups[i]))
				}
			}
			factor.list<-factor.list[,-1,drop=F]
			#lastly, tp is ordered by each of the factors, starting with the first and breaking ties by further arguments
			plot.order<-do.call(order,as.list(factor.list))
			tp<-tp[plot.order,]
			factor.list<-factor.list[plot.order,,drop=F]
			if(length(group.spacing)>0) {
				group.spacing<-as.logical(group.spacing)
				if(anyNA(group.spacing)==TRUE) stop('Spacing arguments must be logical or able to be coerced to logical.')
				if(length(groups)>length(group.spacing)) {
					warning('Spacing arguments should match the number of groups supplied.\nFALSE Spacing codes will be added to match the number of groups supplied.',call.=FALSE)
					diff<-length(groups)-length(group.spacing)
					to.add<-rep(FALSE,diff)
					group.spacing<-c(group.spacing,to.add)
				}
				else if(length(groups)<length(group.spacing)) {
					warning('Spacing arguments should match the number of groups supplied.\nSpacing codes will be removed to match the number of groups supplied.',call.=FALSE)
					groups.spacing<-group.spacing[1:length(groups)]
				}
				composite.grouping.factor<-as.character(factor.list[,1]) 
				j<-2
				while(j<=length(group.spacing)) {
					composite.grouping.factor<-as.factor(paste(composite.grouping.factor,as.character(factor.list[,j]),sep='\n'))
					j<-j+1
				}
			}
		}
		rownames(total.row)<-'total'
		tp<-rbind(tp,total.row)
		melt.tp<-melt(tp,id.var='sample')
		names(melt.tp)[2:3]<-c('taxon','rel_abund')
		total<-which(melt.tp$sample=='total',)
		melt.x<-melt.tp[-total,]
		melt.x$grouping_factor<-as.factor(rep(composite.grouping.factor,nrow(x)))
		melt.x$taxon<-factor(melt.x$taxon,levels=as.ordered(unique(melt.x$taxon)))			
	#Construct plot with all samples..............................................
		sample.colors<-rep(x.colors,length(unique(melt.x$sample)))
	#Spacing between facets is assigned based on user input from the function arguments
		all.factor.combos<-nlevels(melt.x$grouping_factor)
		split.combos<-1
		for(i in 1:length(group.spacing)) {
			if(group.spacing[i]==TRUE) {
				split.combos<-split.combos*nlevels(factor.list[,i])
			}
		}
		spacing.interval<-seq.int(from=all.factor.combos/split.combos,to=(all.factor.combos-1),by=all.factor.combos/split.combos)
		split.factor<-rep(0,all.factor.combos-1)
		split.factor[spacing.interval]<-0.5 #edit here to change spacing in the graph
	#plot
		gg.sample<-ggplot(melt.x,aes(sample,rel_abund,fill=taxon))+
		geom_bar(position='fill',stat='identity',col='black',fill=sample.colors)+
		facet_grid(~grouping_factor,scales='free')+
		ylab('Relative Abundance')+
		theme_classic(base_size=10)+
		theme(legend.position='none',
			panel.spacing.x=unit(split.factor,units='lines'),
			strip.background=element_blank(),
			strip.text.x=element_blank(),
			axis.line=element_blank(),
			axis.title.x=element_blank(),
			axis.ticks.x=element_blank(),
			axis.text=element_text(colour='black'),
			axis.text.x=element_text(angle=label.angle,vjust=0.2,margin=margin(t=-5)))
	}
#Construct plot of summary by phyla.................................................
	if(type=='both'|type=='summary') {
		metric<-match.arg(metric)
		if(metric=='mean') {
			x$mean<-apply(as.matrix(x[,-c(1,(length(x)-1))]),1,mean)
			x$se<-apply(as.matrix(x[,-c(1,(length(x)-1),length(x))]),1,sd)
			x$se<-x$se/sqrt(unique(length(melt.x$sample)))
			#plot
			gg.phyla<-ggplot(x,aes(taxon,mean,fill=taxon))+
			geom_bar(stat='identity',col='black',fill=x.colors)+
			coord_flip()+
			ylab('Relative Abundance')+
			theme_classic(base_size=10)+
			theme(legend.position='none',
				axis.line=element_blank(),
				axis.title.y=element_blank(),
				axis.ticks.y=element_blank(),
				axis.text=element_text(colour='black'))
		}
		else {
			gg.phyla<-ggplot(x,aes(taxon,total,fill=taxon))+
			geom_bar(stat='identity',col='black',fill=x.colors)+
			coord_flip()+
			ylab('Relative Abundance')+
			theme_classic(base_size=10)+
			theme(legend.position='none',
				axis.line=element_blank(),
				axis.title.y=element_blank(),
				axis.ticks.y=element_blank(),
				axis.text=element_text(colour='black'))
		}
	}
#plot.....................................................................
	if(type=='summary') gg.phyla
	else if(type=='sample') gg.sample
	else if(type=='both') {
		grid.arrange(gg.phyla,gg.sample,ncol=2,widths=c(.5,1))
	}
}
