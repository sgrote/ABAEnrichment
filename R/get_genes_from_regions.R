
# input: 
	# 0/1-vector with genomic regions as names (chr:from-to)
	# gene_coords
# output: list with elements
	# test-regions (bed)
	# background_regions (bed)
	# genes-vector ("normal hyper-input" for genes from test-regions) 

get_genes_from_regions = function(genes, gene_pos, circ_chrom){
	
	# convert coordinates from 'genes'-names to bed-format 
	bed = do.call(rbind, strsplit(names(genes), "[:-]"))
#	bed[,1] = substring(bed[,1], 4)  ## if chrom is given as chr21 instead of 21 
	bed = as.data.frame(bed)
	bed[,2:3] = apply(bed[,2:3], 2, as.integer)
	
	# check that start < stop, ##TODO:maybe list the wrong ones
	if(any(bed[,2] > bed[,3])){
		stop("In chr:start-stop start must always be smaller than stop.")
	}	
	
	# TODO test that regions are non-overlapping (separately for candidate and background)
	
	# split in test and background
	test_reg = bed[genes==1,]
	bg_reg = bed[genes==0,]	
	
	# sort  (assuming regions do not overlap)
	test_reg = test_reg[order(test_reg[,1], test_reg[,2]),] 
	bg_reg = bg_reg[order(bg_reg[,1], bg_reg[,2]),] 
		
	# if rolling chrom: remove unused bg chroms and warn, check that all candidate chroms have bg
	if(circ_chrom == TRUE){				
		if(!(all(bg_reg[,1] %in% test_reg[,1]))){
			not_used = paste(unique(bg_reg[!(bg_reg[,1] %in% test_reg[,1]),1]), collapse=", ")
			# TODO: order not used chrom-numbers for message (mixedsort from gtools - add gtools to Dependencies) warum sind die nicht schon sortiert? - weil character ->  1 10 12 2 20 3
			# maybe a sentence like with circ_chrom=T only background of same chrom is used oder so
			warning(paste("Unused chromosomes in background regions: ", not_used, ".",sep=""))
			bg_reg = bg_reg[bg_reg[,1] %in% test_reg[,1],]
		}			
		if(!(all(test_reg[,1] %in% bg_reg[,1]))){
			wo_bg = paste(unique(test_reg[!(test_reg[,1] %in% bg_reg[,1]),1]), collapse=", ")
			stop(paste("No background region for chromosome: ",  wo_bg, ".",sep=""))
		}
	} else {  # normal blocks option
		# check that biggest bg_region is bigger than biggest test_region
		if (max(bg_reg[,3] - bg_reg[,2]) < max(test_reg[,3] - test_reg[,2])){
			stop("At least one test region is bigger than any background region.")
		}
		# sort candidate regions by length (better chances that random placement works with small bg-regions)
		test_reg = test_reg[order(test_reg[,3] - test_reg[,2], decreasing=T),]
	}
	
	# get genes overlapping mappable-regions
	bg_genes = c()
	for (i in 1:nrow(bg_reg)){
		bg_genes = c(bg_genes, gene_pos[gene_pos[,"chr"]==bg_reg[i,1] & ((gene_pos[,"start"] >= bg_reg[i,2] & gene_pos[,"start"] < bg_reg[i,3]) | (gene_pos[,"end"] >= bg_reg[i,2] & gene_pos[,"end"] < bg_reg[i,3]) |  (gene_pos[,"start"] <= bg_reg[i,2] & gene_pos[,"end"] >= bg_reg[i,3])), "hgnc_symbol"])		
	}
	
	# get genes overlapping test-regions
	test_genes = c()
	for (i in 1:nrow(test_reg)){
		test_genes = c(test_genes, gene_pos[gene_pos[,"chr"]==test_reg[i,1] & ((gene_pos[,"start"] >= test_reg[i,2] & gene_pos[,"start"] < test_reg[i,3]) | (gene_pos[,"end"] >= test_reg[i,2] & gene_pos[,"end"] < test_reg[i,3]) | (gene_pos[,"start"] <= test_reg[i,2] & gene_pos[,"end"] >= test_reg[i,3])), "hgnc_symbol"])		
	}

	# convert to classic "genes" func-input-vector 
	gene_names = unique(c(test_genes, bg_genes))
	genes_vec = rep(0, length(gene_names))
	names(genes_vec) = gene_names
	genes_vec[names(genes_vec) %in% test_genes] = 1
	
	return(list(test_reg, bg_reg, genes_vec))	
}
	
	
	
	
	
	
	
	
