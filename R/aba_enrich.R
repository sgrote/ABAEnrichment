# run "FUNC" with gene expression data from Allen Brain Atlas 
# over flexible expression cutoff
# wilcoxon rank test or hypergeometric test
# with adult microarray data, developmental RNA-seq (5 age categories) or derived developmental effect score


###############
# MAIN function
###############

### Arguments
# genes: vector of 0/1 (hypergeometric) or float (wilcoxon) with gene identifiers as names (ensembl, entrez or hgnc), or chromosomal regions chr:start-stop (only for hypergeometric)
# dataset: "adult", "5_stages" or "dev_effect"
# test: "hyper" or "wilcoxon"
# cutoff_quantiles: numeric values in ]0,1[
# n_randsets: number of random-sets
# gene_len: randomset is dependent on length of genes
# circ_chrom: for regions input: random regions are on same chrom and allowed to overlap multiple bg-regions

#######################################################################################
# 1. check arguments and define parameters
# 2. load expression data
# 3. loop over age categories and cutoffs
#	3a) create input
#	3b) run func
#	3c) summarize func output
# 4. create output
#######################################################################################

aba_enrich=function(genes,dataset="adult",test="hyper",cutoff_quantiles=seq(0.1,0.9,0.1),n_randsets=1000, gene_len=FALSE, circ_chrom=FALSE){


	ref_genome = "grch37" # TODO: add grch38 as parameter and gene_coords_grch38 to sysdata.R 
	
	#####	1. Check arguments and define parameters
		
	## Define dataset dependent parameters
	if (dataset=="adult"){
		folder_ext="adult"
		root_node="Allen:4005"
	} else if (dataset=="5_stages"){
		folder_ext="developmental"
		root_node="Allen:10153"
	} else if (dataset=="dev_effect"){
		folder_ext="developmental"
		root_node="Allen:10153"
	} else {
		stop("Not a valid dataset. Please use 'adult', '5_stages' or 'dev_effect'.")
	}
	
	## Check arguments
	# general
	message("Checking arguments...")
	if (length(genes)==0){
		stop("Please enter genes.")
	}	
	if (length(names(genes))==0){
		stop("Please add gene identifiers as names to 'genes' vector.")
	}
	if (class(cutoff_quantiles)!="numeric"){
		stop("Please enter numeric cutoff_quantiles.")
	}	
	if(min(cutoff_quantiles)<=0 | max(cutoff_quantiles)>=1){
		stop("Please enter cutoff_quantiles between 0 and 1 (exclusive).")
	}			
	valid_datasets = c("adult","5_stages","dev_effect")
	if(!(dataset %in% valid_datasets)){
		stop("Not a valid dataset. Please use 'adult', '5_stages' or 'dev_effect'.")	
	}
	if(length(n_randsets)>1 || !is.numeric(n_randsets) || n_randsets<1){
		stop("Please define 'n_randsets' as a positive integer.")
	}
	if(n_randsets != round(n_randsets)){
		n_randsets = round(n_randsets)
		warning(paste("'n_randsets' is expected to be an integer and was rounded to ",n_randsets,".",sep=""))
	}
	# test-specific arguments
	if (test=="hyper"){
		if(!all(genes %in% c(0,1))){
			stop("Not a valid 'genes' argument for hypergeometric test. Please use a vector of 0/1.")	
		}	
		if(sum(genes)==0){
			stop("Only 0 in genes vector. Please enter test genes.")	
		}
	} else	if (test=="wilcoxon"){
		if(!is.numeric(genes)){
			stop("Not a valid 'genes' argument. Please use a numeric vector.")	
		}	
	} else (stop("Not a valid test. Please use 'hyper' or 'wilcoxon'."))
	
	# Create input and output folder
	directory = tempdir()
	
	# load gene_list and get gene identifier
	gene_symbols = get(paste("gene_symbols",folder_ext,sep="_"))	
	blocks = FALSE
	identifier = detect_identifier(names(genes)[1])	
	
	# load gene_coords
	gene_coords = get(paste("gene_coords_", ref_genome, sep=""))
	
	if (identifier=="blocks"){
		identifier = "hgnc_symbol" # gene-name 
		blocks = TRUE
		# check that background region is specified
		if(sum(genes) == length(genes)){
			stop("All values of the 'genes' input are 1. Using chromosomal regions as input requires defining background regions with 0.")
		}
		# check that test is hyper
		if (test != "hyper"){
			stop("chromosomal regions can only be used with test='hyper'.")
		}	
		# warn if gene_len=TRUE, although regions are used
		if (gene_len == TRUE){
			warning("Unused argument: 'gene_len = TRUE'.")
		}	
		
		# convert coords from genes-vector to bed format, SORT, and extract genes overlapping test regions
		regions = get_genes_from_regions(genes, gene_coords, circ_chrom) # gene_coords from sysdata.rda HGNC
		test_regions = regions[[1]]		
		bg_regions = regions[[2]]
		genes = regions[[3]]

		message("Candidate regions:")
		print(test_regions)
		message("Background regions:")
		print(bg_regions)		
	
		# write regions to files
		write.table(test_regions,file=paste(directory, "/test_regions.bed",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
		write.table(bg_regions,file=paste(directory, "/bg_regions.bed",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")		
	} else if (circ_chrom == TRUE){
		# warn if circ_chrom=TRUE, although individual genes are used
		warning("Unused argument: 'circ_chrom = TRUE'.")
	} 	
	
	gene_list = gene_symbols[,identifier]	
	
	# check which genes dont have expression data annotated
	not_in = names(genes)[!(names(genes) %in% gene_list)]
	remaining_genes = genes[names(genes) %in% gene_list]
	
	if (length(not_in)>0){		
		# is any gene in the data? if not, stop.
		if (length(remaining_genes) == 0) {
			stop("No requested genes in data.")
		}
		# for 0/1 data: are test genes and background genes in the data (background only if background specified)?
		if (test=="hyper" & sum(remaining_genes)==0) {
			stop("No requested test genes in data.")
		}
		if (test=="hyper" & sum(genes)!=length(genes) & sum(remaining_genes)==length(remaining_genes)) {
			stop("No requested background genes in data.")
		}
		if(!blocks){ # this message is usually too long when blocks are used. 
			not_in_string = paste(not_in,collapse=", ")
			warning(paste("No expression data for genes: ",not_in_string,".\n  These genes were not included in the analysis.",sep=""))
		}
	}	


	#####	2. Load expression data
	
	# load pre input
	message("Loading dataset...")
	pre_input = load_by_name(paste("dataset",dataset,sep="_"))
	
	# TODO: save cutoff quantiles for all datasets for default cutoff_quantile values.
		
	 # compute cutoffs (tapply lasts much longer - use only when needed)
	 message("Computing cutoffs... ")
	 cutoff_quantiles = sort(cutoff_quantiles)
	 if(dataset=="5_stages"){
	 	cutoff_list = tapply(pre_input$signal, pre_input$age_category, function(x) quantile(x, probs=cutoff_quantiles),simplify=FALSE)
	 } else {	
	 	cutoff_list = list()
	 	cutoff_list[[1]] = quantile(pre_input$signal,probs=cutoff_quantiles)
	 	names(cutoff_list) = pre_input$age_category[1]
	 }

	# select gene identifier
	pre_input = pre_input[,c("age_category",identifier,"structure","signal")]
	colnames(pre_input)[2] = "gene_id"

	# restrain input to required genes (unless test=hypergeometric and background genes are not defined - all other genes are background in this case)
	# else exlude NAs if entrez-ID, (hgnc and ensembl dont have NAs)
	if(!(test=="hyper" & sum(genes) == length(genes))){		
		message("Select requested genes...")
		pre_input = pre_input[pre_input[,2] %in% names(genes),]
	} else if (identifier == "entrezgene"){
		message("Exclude NAs...")	
		pre_input = pre_input[!is.na(pre_input$gene_id),]
	}			
	# aggregate expression per gene (duplicates due to different identifiers)
	message("Checking for gene duplicates to aggregate... ")
	# simple aggregate would take very long for large input vectors
	# number of structure-age-combinations: if all structures are present in all ages they can simply be multiplied -> thats the case
	# n=nrow(unique(pre_input[,c("age_category","structure")]))
	n = length(unique(pre_input$age_category))*length(unique(pre_input$structure))
	ngenes = table(pre_input$gene_id)
	bose = ngenes[ngenes>n]
	if(length(bose)>0){
		message(" Aggregate expression per gene...")
		part1 = pre_input[pre_input$gene_id %in% names(bose),]
		part2 = pre_input[!(pre_input$gene_id %in% names(bose)),]
		part1 = aggregate(part1$signal,by=list("age_category"=part1$age_category,"gene_id"=part1$gene_id,"structure"=part1$structure),mean)
		colnames(part1)[4]="signal"
		pre_input = rbind(part1,part2)
		message(" Done.")
	} else {message(" No gene duplicates - no aggregation needed.")}		
	
	# load ontology and write files to tmp
	term = get(paste("term",folder_ext,sep="_"))
	term2term = get(paste("term2term",folder_ext,sep="_"))
	graph_path = get(paste("graph_path",folder_ext,sep="_"))
	write.table(term,file=paste(directory, "/term.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	write.table(term2term,file=paste(directory, "/term2term.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	write.table(graph_path,file=paste(directory, "/graph_path.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

	# initialize for warning message if candidate genes have no coordinates
	candi_no_coords = c()

	#####	3. loop over age categories and cutoffs
	# initialize output
	out = data.frame()
	for (s in unique(pre_input$age_category)){
		stage_input = pre_input[pre_input$age_category==s,2:4]	
		# get cutoffs		
		cutoff = cutoff_list[[as.character(s)]]
		colnames(stage_input)[3] = "signal"		
		for (i in 1:length(cutoff)){
			
			##	3a) create input

			if (dataset=="5_stages"){
				message(paste("Preparing age category ", s," with ",names(cutoff[i]),"-quantile expression cutoff for FUNC...",sep=""))			
			} else {
				message(paste("Preparing data with ",names(cutoff[i]),"-quantile expression cutoff for FUNC...",sep=""))	
			}
			# restrain input with cutoff
			message(" Apply cutoff...")
			input = stage_input[stage_input$signal >= cutoff[i],]			
		
			# for Hypergeometric Test: 0 for all genes, then convert to 1 for interesting genes (merge not possible if no background genes are defined)
			message(" Check that there are sufficient genes above cutoff...")
			breaky = FALSE
			if (test=="hyper"){	
				# do input data have both 0 and 1? else break (FUNC would throw error and summary-file would not be generated)
				# now that cutoffs get computed for all genes (not just input) both could be missing
				if (nrow(input)==0){ 
					message(paste("  Warning: No statistics for cutoff >= ",cutoff_quantiles[i],". No input gene expression above cutoff.",sep=""))
					breaky = TRUE
				} else {	
					input$signal = 0
					interesting_genes = names(remaining_genes)[remaining_genes==1]
					input[input[,1] %in% interesting_genes,"signal"] = 1
					if (sum(input$signal)==0){
						message(paste("  Warning: No statistics for cutoff >= ",cutoff_quantiles[i],". No test gene expression above cutoff.",sep=""))
						breaky = TRUE
					} else if (sum(input$signal)==nrow(input)){
						message(paste("  Warning: No statistics for cutoff >= ",cutoff_quantiles[i],". No background gene expression above cutoff.",sep=""))
						breaky = TRUE
					}
				}		
			# for Wilcoxon test: replace expression signal with score
			} else	if (test=="wilcoxon"){
				if(length(unique(input[,1]))<2){
					message(paste("  Warning: No statistics for cutoff >= ",cutoff_quantiles[i],". Less than 2 genes show expression above cutoff.",sep=""))
					breaky = TRUE
				}
				if (!breaky){
					score = genes[match(input$gene_id,names(genes))]
					input$signal = score
				}		
			}	
			# if first cutoff fails, no output produced, else return output from previous cutoffs
			if (breaky){
				if(i==1){
					stop ("Expression cutoffs too high.")
				} else {
					break
				}	
			}
			message(" Rearrange input for FUNC...")
			# add "Allen:" string to structures (shouldn't be numeric)
			input$structure=paste("Allen:",input[,2],sep="")	
			
			# prepare input data (infile-data and root like in separate_taxonomies.pl)
			
			# "infile-data": genes and associated scores (wilcox) 
			# "Allen:4005-changed": one column with test genes (hyper)
			if (test=="hyper"){
				# subset to test genes
				infile_data = data.frame(genes=unique(input[input[,3]==1,1]))
			} else if(test=="wilcoxon"){
				infile_data = unique(input[,c(1,3)])
			}
	
			# create "root" dataframe, add genomic positions if block option is used
			# add gene-coordinates, despite for classic FUNC option (user defined genes might get lost because they have no annotated position - there are some from ABA which are not in gene_coords)			
			if (blocks || gene_len){
				# each row contains gene and a string with all annotaded brain regions
				xx=tapply(input[,2],input[,1],function(x){paste(x,collapse=" ")}) # paste annotations
				gene = as.character(names(xx))
				# add coordinates
				gene_position = gene_coords[match(gene, gene_coords[,identifier]),4:6]	
				root = data.frame(genes=gene, gene_position ,goterms=as.character(xx))				
			} else {
				# each row contains gene and a string with all annotaded brain regions
				xx = tapply(input[,2],input[,1],function(x){paste(x,collapse=" ")}) # paste annotations
				root = data.frame(genes=as.character(names(xx)),goterms=as.character(xx))
			}
			# remove genes with unknown coordinates (possible in gene_len-option) and 
			if (gene_len){					
				# warn if this affects test genes (background genes might be undefined and then its weird to get a warning about them)
				no_coords = root[is.na(root[,3]),1]
				candi_no_coords = c(candi_no_coords, as.character(no_coords[no_coords %in% infile_data[,1]]))
				candi_no_coords = unique(candi_no_coords)	
				root = root[!is.na(root[,3]),] 			
			}	
		
			# write input-files to tmp-directory
			# last priority rename root and infile-data files to all_genes_annotation and test_genes
			write.table(infile_data,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste(directory,"/infile-data",sep=""))
			write.table(root,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste(directory,"/",root_node,sep=""))
			
			###	3b) run func
			if (test=="hyper"){
				# randomset
				if (blocks & circ_chrom){
					hyper_randset(paste(directory,"/",root_node,sep=""), n_randsets, directory, root_node, "roll")
				} else if (blocks){
					hyper_randset(paste(directory,"/",root_node,sep=""), n_randsets, directory, root_node,"block")				
				} else if (gene_len){						
					hyper_randset(paste(directory,"/",root_node,sep=""), n_randsets, directory, root_node, "gene_len")	
				} else {						
					hyper_randset(paste(directory,"/",root_node,sep=""), n_randsets, directory, root_node, "classic")
				}	
#				stop("only randomset")
				# category test
				hyper_category_test(paste(directory, "/randset_out",sep=""), paste(directory,"/category_test_out", sep=""), 1, root_node)
			} else if (test=="wilcoxon"){
				wilcox_randset(paste(directory,"/",root_node,sep=""), n_randsets, directory, root_node)
				# category test
				wilcox_category_test(paste(directory, "/randset_out",sep=""), paste(directory,"/category_test_out", sep=""), 1, root_node)
			}
		
			## 3c) summarize func-output
			groupy = read.table(paste(directory,"/category_test_out",sep=""))
			age_category = rep(s,nrow(groupy))			
			cutoff_quantile = rep(cutoff_quantiles[i],nrow(groupy))
			cutoff_value = rep(cutoff[i],nrow(groupy))
			groups = cbind(age_category,cutoff_quantile,cutoff_value,groupy)
				
			# combine with output from previous cutoffs
			out = rbind(out,groups)			
		} # end cutoffs		
	} # end ages
	
	if (length(candi_no_coords) > 0){
			no_coords_string = paste(candi_no_coords,collapse=", ")
			warning(paste("No coordinates available for candidate genes: ",no_coords_string,".\n  These genes were not included in the analysis.",sep=""))
	}	
	
	#####. 4 rearrange output
	message("Creating output...")
	
	if (test == "hyper"){
		colnames(out)[4:8] = c("node_id","raw_p_underrep","raw_p_overrep","FWER_underrep","FWER_overrep")
	} else if (test  == "wilcoxon"){
		colnames(out)[4:8] = c("node_id","raw_p_low_rank","raw_p_high_rank","FWER_low_rank","FWER_high_rank")
	}
	# remove redundant rows (nodes with no data and only one child)
	cluster = get(paste("node_clusters",folder_ext,sep="_"))
	results = rearrange_output(out,cluster,term)
	# rearrange cutoff-list
	cutoffs = do.call(cbind.data.frame, cutoff_list)
	colnames(cutoffs) = paste("age_category",colnames(cutoffs),sep="_")

	# save in package environment: dataframes for test and background-gene expression, test, dataset	
	# (to be returned with get_expression(structure_ids))	
	requested_gene_expression = pre_input
	rownames(requested_gene_expression) = NULL
	colnames(requested_gene_expression)[3] = "structure_id"
	requested_gene_expression[,3] = paste("Allen:",requested_gene_expression[,3],sep="")
	remember = list(rge=requested_gene_expression, test=c(test,dataset), genes=genes)
	# save to package environment
	aba_env = as.environment("package:ABAEnrichment")
	unlock_environment(aba_env)
	if (exists("remember",where=aba_env)){
		unlockBinding("remember",aba_env)
	}	
	aba_env$remember = remember
	lockEnvironment(aba_env, bindings=TRUE)

	final_output = list(results=results, genes=remaining_genes, cutoffs=cutoffs)
	
	message("\nDone.\n")	
	return(final_output)	
}

	
