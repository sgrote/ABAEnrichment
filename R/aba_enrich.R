# run "FUNC" with gene expression data from Allen Brain Atlas 

#######################################################################################
# 1. check arguments and define parameters
# 2. load expression data
# 3. loop over age categories and cutoffs
#	3a) create input
#	3b) run func
#	3c) summarize func output
# 4. create output
#######################################################################################

aba_enrich=function(genes, dataset="adult", test="hyper", cutoff_quantiles=seq(0.1,0.9,0.1), n_randsets=1000, gene_len=FALSE, circ_chrom=FALSE, ref_genome="grch37", silent=FALSE){


	#####	1. Check arguments and define parameters
	
	## still allow vector as input for hyper and wilcox (like in older versions)
    if (!((is.vector(genes) && !is.null(names(genes))) || is.data.frame(genes))){
        stop("Please provide a data frame as 'genes' input (also named vector is still accepted for hypergeometric or wilcoxon rank-sum test).")
    }
    if (is.vector(genes)){
        genes = data.frame(gene=names(genes), score=unname(genes), stringsAsFactors=FALSE)
    } else {
        genes[,1] = as.character(genes[,1])
    }
		
	## Define dataset dependent parameters
	if (dataset=="adult"){
		folder_ext="adult"
		root_id="Allen:4005"
	} else if (dataset=="5_stages"){
		folder_ext="developmental"
		root_id="Allen:10153"
	} else if (dataset=="dev_effect"){
		folder_ext="developmental"
		root_id="Allen:10153"
	} else {
		stop("Not a valid dataset. Please use 'adult', '5_stages' or 'dev_effect'.")
	}
	
	## Check arguments
	if (!silent) message("Checking arguments...")
	# genes
    if (anyNA(genes)){
        stop("NAs found in 'genes'-input. Please provide finite values only.")
    }
    if (test=="hyper"){
        if (!(is.data.frame(genes) && ncol(genes)==2 && is.numeric(genes[,2]))){
            stop("Please provide a data frame with 2 columns [gene, 0/1] as input for hypergeometric test.")
        }
        if (!(all(genes[,2] %in% c(1,0)))){
            stop("Please provide only 1/0-values in 2nd column of 'genes'-input for hypergeometric test.")
        }
        if (sum(genes[,2])==0){
            stop("Only background genes defined (only 0s in 'genes[,2]'). Please provide candidate genes denoted by 1.")
        }
    } else  if (test=="wilcoxon"){
        if (!(is.data.frame(genes) && ncol(genes)==2  && is.numeric(genes[,2]))){
            stop("Please provide a data frame with 2 columns [gene, score] as input for wilcoxon rank-sum test.")
        }
        if (nrow(genes) < 2){
            stop("Only one gene provided as input.")
        }
    } else if (test=="binomial"){
        if (!(is.data.frame(genes) && ncol(genes)==3 && all(sapply(genes,is.numeric) == c(0,1,1)))){
            stop("Please provide a data frame with columns [gene, count1, count2] as input for binomial test.")
        }
        if (any(genes[,2:3] < 0) | any(genes[,2:3] != round(genes[,2:3]))){
            stop("Please provide non-negative integers in columns 2-3.")
        }
    } else if (test=="contingency"){
        if (!(is.data.frame(genes) && ncol(genes)==5 && all(sapply(genes,is.numeric) == c(0,1,1,1,1)))){
            stop("Please provide a data frame with columns [gene, count1A, count2A, count1B, count2B] as input for contingency table test.")
        }
        if (any(genes[,2:5] < 0) | any(genes[,2:5] != round(genes[,2:5]))){
            stop("Please provide non-negative integers in columns 2-5.")
        }
    } else {
        stop("Not a valid test. Please use 'hyper', 'wilcoxon', 'binomial' or 'contingency'.")
    }
    # check for multiple assignment for one gene
    genes = unique(genes) # allow for multiple assignment of same value
    multi_genes = sort(unique(genes[,1][duplicated(genes[,1])]))
    if (length(multi_genes) > 0){
        stop(paste("Genes with multiple assignment in input:", paste(multi_genes,collapse=", ")))
    }
    
	# other arguments
	if (class(cutoff_quantiles)!="numeric"){
		stop("Please enter numeric cutoff_quantiles.")
	}	
	if (min(cutoff_quantiles)<=0 | max(cutoff_quantiles)>=1){
		stop("Please enter cutoff_quantiles between 0 and 1 (exclusive).")
	}			
	valid_datasets = c("adult","5_stages","dev_effect")
	if (!(dataset %in% valid_datasets)){
		stop("Not a valid dataset. Please use 'adult', '5_stages' or 'dev_effect'.")	
	}
    if (length(n_randsets)>1 || !is.numeric(n_randsets) || n_randsets<1){
        stop("Please define 'n_randsets' as a positive integer.")
    }
    if (n_randsets != round(n_randsets)){
        n_randsets = round(n_randsets)
        warning(paste("'n_randsets' is expected to be an integer and was rounded to ",n_randsets,".",sep=""))
    }
    if (!is.logical(gene_len)){
        stop("Please set gene_len to TRUE or FALSE.")
    }
    if (!is.logical(circ_chrom)){
        stop("Please set circ_chrom to TRUE or FALSE.")
    }
    if (gene_len & !(test == "hyper")){
        stop("Argument 'gene_len = TRUE' can only be used with 'test = 'hyper''.")
    }
	if (!(ref_genome %in% c("grch37", "grch38"))){
		stop("Not a valid ref_genome. Please use 'grch37' or 'grch38'.")
	}
	
	# Create tempfile prefix (in contrast to tempdir() alone, this allows parallel processing)
	directory = tempfile()
#	dir.create("tempdir"); directory = paste("./tempdir/tempfile",Sys.info()["nodename"],sep="_")
	
	# load gene_list and get gene identifier
	gene_symbols = get(paste("gene_symbols",folder_ext,sep="_"))	
	blocks = FALSE
	identifier = detect_identifier(genes[1,1])	
	
	# load gene_coords
	gene_coords = get(paste("gene_coords_", ref_genome, sep=""))
	
	if (identifier == "blocks"){
		identifier = "hgnc_symbol" # gene-name 
		blocks = TRUE
		# check that background region is specified
		if (sum(genes[,2]) == nrow(genes)){
			stop("All values of the 'genes' input are 1. Using chromosomal regions as input requires defining background regions with 0.")
		}
		# check that test is hyper
		if (test != "hyper"){
			stop("chromosomal regions can only be used with test='hyper'.")
		}	
		# warn if gene_len=TRUE, although regions are used
		if (gene_len){
			warning("Unused argument: 'gene_len = TRUE'.")
		}	
		
		# convert coords from genes-vector to bed format, SORT, and extract genes overlapping test regions
		regions = get_genes_from_regions(genes, gene_coords, circ_chrom) # gene_coords from sysdata.rda HGNC
		test_regions = regions[[1]]		
		bg_regions = regions[[2]]
		genes = regions[[3]]

		# avoid scientific notation in regions (read in c++)
		test_regions = format(test_regions, scientific=FALSE, trim=TRUE)
		bg_regions = format(bg_regions, scientific=FALSE, trim=TRUE)

		if (!silent){
			message("Candidate regions:")
			print(test_regions)
			message("Background regions:")
			print(bg_regions)
		}
	
		# write regions to files
		write.table(test_regions,file=paste(directory, "_test_regions.bed",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
		write.table(bg_regions,file=paste(directory, "_bg_regions.bed",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")		
	} else {
		if (circ_chrom == TRUE){
			# warn if circ_chrom=TRUE, although individual genes are used
			warning("Unused argument: 'circ_chrom = TRUE'.")
		}
	}
	if (ref_genome == "grch38" & !(blocks | gene_len)){
		# warn if ref-genome is switched to hg20 but not used
		warning("Unused argument: 'ref_genome = grch38'.")
	}

	gene_list = gene_symbols[,identifier]	

	# restrict to genes that have expression data annotated and warn about the rest
	gene_values = genes[genes[,1] %in% gene_list, ] # restrict
	not_in = unique(genes[!genes[,1] %in% gene_values[,1], 1]) # removed
	if (length(not_in) > 0 && !blocks){
		not_in_string = paste(not_in,collapse=", ")
		warning(paste("No expression data for genes: ",not_in_string,".\n  These genes were not included in the analysis.",sep=""))
	}
	# restrict to genes that have coordinates and warn about the rest
	if (gene_len){
		# removed
		no_coords = unique(gene_values[!(gene_values[,1] %in% gene_coords[,identifier]), 1]) 
		if (length(no_coords) > 0){
			gene_values = gene_values[!(gene_values[,1] %in% no_coords), ] # restrict
			no_coords_string = paste(no_coords,collapse=", ")
			warning(paste("No coordinates available for genes: ",no_coords_string,".\n  These genes were not included in the analysis.",sep=""))
			not_in = c(not_in, no_coords)
		}
	}
	# after removing genes without expression data or coordinates: are enough genes left?
	if (length(not_in) > 0){
		# is any gene in the data? if not, stop.
		if (nrow(gene_values)==0) {
			stop("No requested genes in data.")
		}
		# for 0/1 data: are test genes and background genes in the data (background only if background specified)?
		if (test=="hyper" && sum(gene_values[,2])==0) {
			stop("No requested test genes in data.")
		}
		if (test=="hyper" && 0 %in% genes[,2] && all(gene_values[,2]==1)) {
			stop("No requested background genes in data.") # although background was defined
		}
		# at least two for wilcoxon
        if (test=="wilcoxon" && nrow(gene_values) < 2) {
            stop(paste("Less than 2 genes have annotated GO-categories.",sep=""))
        }
	}
	# candidate genes	
	if (test=="hyper"){					
		interesting_genes = gene_values[gene_values[,2]==1, 1] 
	}

	#####	2. Load expression data
	
	# load pre input
	if (!silent) message("Loading dataset...")
	gene_expr = load_by_name(paste("dataset",dataset,sep="_"), silent)
	gene_expr = as.data.table(gene_expr)
	
	# compute cutoffs
	if (!silent) message("Computing cutoffs... ")
	cutoff_quantiles = sort(cutoff_quantiles)
	cutoffs = gene_expr[,list(q=quantile(signal, probs=cutoff_quantiles)), by=age_category]
	# write to a list (to be consistent with tapply-approach that was previously used)
	cutoff_list = tapply(cutoffs$q, cutoffs$age_category, function(x){
		names(x) = paste(100*cutoff_quantiles, "%", sep=""); return(x)}, simplify=FALSE)

	# select gene identifier
	gene_expr = gene_expr[,c("age_category",identifier,"structure","signal"), with=FALSE]
	colnames(gene_expr)[2] = "gene_id"

	# restrain input to required genes (unless test=hypergeometric and background genes are not defined - all other genes are background in this case)
	# else exlude NAs if entrez-ID, (hgnc and ensembl dont have NAs)
	if (!(test=="hyper" & sum(genes[,2]) == nrow(genes))){		
		if (!silent) message("Select requested genes...")
		gene_expr = gene_expr[gene_id %in% gene_values[,1],]
	} else if (identifier == "entrezgene"){
		if (!silent) message("Exclude NAs...")
		gene_expr = gene_expr[!is.na(gene_expr$gene_id),]
	}

	# aggregate expression per gene (duplicates due to different identifiers)
	if (!silent) message("Aggregate expression per gene...")
	gene_expr_ag = gene_expr[,mean(signal), by=list(age_category, gene_id, structure )]
	colnames(gene_expr_ag)[4] = "signal"

	# Add "Allen:"-string to brain region IDs
	gene_expr_ag[,structure:=paste("Allen:", gene_expr_ag[,structure], sep="")]
	
	# load ontology and write files to tmp
	term = get(paste("term",folder_ext,sep="_"))
	term2term = get(paste("term2term",folder_ext,sep="_"))
	graph_path = get(paste("graph_path",folder_ext,sep="_"))
	write.table(term,file=paste(directory, "_term.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	write.table(term2term,file=paste(directory, "_term2term.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	write.table(graph_path,file=paste(directory, "_graph_path.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

	#####	3. loop over age categories and cutoffs
	# initialize output
	out = data.frame()
	for (s in unique(gene_expr_ag$age_category)){
		stage_expr = gene_expr_ag[gene_expr_ag$age_category==s,2:4] # gene, structure, expression
		# get cutoffs		
		cutoff = cutoff_list[[as.character(s)]]
		for (i in 1:length(cutoff)){
			
			##	3a) create input
			if (!silent){ 
				if (dataset=="5_stages"){
					message(paste("Preparing age category ", s," with ",names(cutoff[i]),"-quantile expression cutoff for FUNC...",sep=""))			
				} else {
					message(paste("Preparing data with ",names(cutoff[i]),"-quantile expression cutoff for FUNC...",sep=""))	
				}
			}
			# restrain input with cutoff
			if (!silent) message(" Apply cutoff...")
			# expressed_genes = expression for input genes in structures where it's above cutoff
			# gene, structure, expression
			expressed_genes = as.data.frame(stage_expr[stage_expr$signal >= cutoff[i],])
			
	        ### prepare input data ('infile-data' and 'root' in c++ scripts)
			
			# 'infile-data'
			# gene | value1 | value2 ...
			# subset genes-input to genes expressed at this age/cutoff
			infile_data = gene_values[gene_values[,1] %in% expressed_genes[,1],] # gene|value1|value2...
			# check if cutoff too high (wilcox needs 2, other tests 1 gene (but not meaningful...))
			if (!silent){
				message(" Check that there are sufficient genes above cutoff...")
			}
			breaky = FALSE
			if (test != "hyper" && nrow(infile_data) < 2){  # hyper can have one candi-gene and default bg
				breaky = TRUE
				if (!silent){
					warning(paste("No statistics for cutoff >= ",cutoff_quantiles[i],". Less than two genes have expression above cutoff.",sep=""))
				}
			}
			# for hyper infile-data is just a column with candidate genes
			# for hyper also check if candidate and background are present
			if (test=="hyper"){
				infile_data = infile_data[infile_data[,2]==1, 1]
				if (length(infile_data)==0){
					breaky = TRUE
					warning(paste("No statistics for cutoff >= ",cutoff_quantiles[i],". No candidate gene expression above cutoff.",sep=""))
				}
				else if (length(infile_data)==nrow(expressed_genes)){
					breaky = TRUE
					warning(paste("No statistics for cutoff >= ",cutoff_quantiles[i],". No background gene expression above cutoff.",sep=""))
				}
			}
			if (breaky){
				# finish at this cutoff or throw error if this is first (lowest) cutoff
				if (i==1){
					stop("Expression cutoffs too high.")
				} else {
					break
				}	
			}
	
			# 'root'
			# gene | (chrom | start | end) | Allen:1 Allen:2 Allen:3
			anno = tapply(expressed_genes[,2], expressed_genes[,1], function(x) {paste(x,collapse=" ")})
			gene = as.character(names(anno))
			if (blocks || gene_len){
				# add coordinates
				gene_position = gene_coords[match(gene, gene_coords[,identifier]),1:3]
				root = data.frame(genes=gene, gene_position, nodes=as.character(anno))
				# remove genes with unknown coordinates (possible for default bg in gene_len-option) 
				if (gene_len){					
					root = root[!is.na(root[,3]),] 			
				}
				# just in case - avoid scientific notation of gene coordinates	
				root[,3:4] = format(root[,3:4], scientific=FALSE, trim=TRUE)				
			} else {
				root = data.frame(genes=gene, nodes=as.character(anno))
			}
			# write input-files to tmp-directory
			# last priority rename root and infile-data files to all_genes_annotation and test_genes
			write.table(infile_data,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste(directory,"_infile-data",sep=""))
			write.table(root,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste(directory,"_",root_id,sep=""))

			###	3b) run func
			if (!silent) message("Run Func...\n")
			if (test == "hyper"){
				# randomset
				if (blocks & circ_chrom){
					hyper_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, root_id, "roll" , silent)
				} else if (blocks){
					hyper_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, root_id, "block", silent)
				} else if (gene_len){
					hyper_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, root_id, "gene_len", silent)
				} else {
					hyper_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, root_id, "classic", silent)
				}
				# category test
				hyper_category_test(directory, 1, root_id, silent)
			} else if (test == "wilcoxon"){
				wilcox_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, root_id, silent)
				wilcox_category_test(directory, 1, root_id, silent)
			} else if (test == "binomial"){
				binom_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, root_id, silent)
				binom_category_test(directory, 1, root_id, silent)          
			} else if (test == "contingency"){
				conti_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, root_id, silent)
				conti_category_test(directory, 1, root_id, silent)
			}

			## 3c) summarize func-output
			groupy = read.table(paste(directory,"_category_test_out",sep=""))
			# NEW: switch columns for binomial test high B - then high A to be more consistent
			if (test == "binomial"){
				groupy[,2:5] = groupy[,c(3,2,5,4)]
			}
			# NEW remove expected and actual no. of candidate, ranksum (included in c++ files to be the same for go_enrich)
			groupy = groupy[,1:5]
			age_category = rep(s,nrow(groupy))			
			cutoff_quantile = rep(cutoff_quantiles[i],nrow(groupy))
			cutoff_value = rep(cutoff[i],nrow(groupy))
			groups = cbind(age_category,cutoff_quantile,cutoff_value,groupy)
				
			# combine with output from previous cutoffs
			out = rbind(out,groups)			
		} # end cutoffs		
	} # end ages
	 
	#####. 4 rearrange output
	if (!silent) message("Creating output...")
	
	print(head(out))
	
#	if (test == "hyper"){
#		colnames(out)[4:8] = c("node_id","raw_p_underrep","raw_p_overrep","FWER_underrep","FWER_overrep")
#	} else if (test  == "wilcoxon"){
#		colnames(out)[4:8] = c("node_id","raw_p_low_rank","raw_p_high_rank","FWER_low_rank","FWER_high_rank")
#	} else if (test == "binomial"){
#        colnames(out)[4:8]=c("node_id","raw_p_high_B","raw_p_high_A","FWER_high_B","FWER_high_A")
#    } else if (test == "contingency"){
#        colnames(out)[4:8]=c("node_id","raw_p_high_CD","raw_p_high_AB","FWER_high_CD","FWER_high_AB")
#    }
#    colnames(out)[4] = "node_id"
	
	# remove redundant rows (nodes with no data and only one child)
	cluster = get(paste("node_clusters",folder_ext,sep="_"))
	results = rearrange_output(out,cluster,term)
	# rearrange cutoff-list
	cutoffs = do.call(cbind.data.frame, cutoff_list)
	colnames(cutoffs) = paste("age_category",colnames(cutoffs),sep="_")
	# sort genes alphabetically (useful for regions input)
	gene_values = gene_values[mixedorder(gene_values[,1]),]
	
	# save in package environment: dataframes for test and background-gene expression, test, dataset	
	# (to be returned with get_expression(structure_ids))	
	requested_gene_expression = as.data.frame(gene_expr_ag)
	rownames(requested_gene_expression) = NULL
	colnames(requested_gene_expression)[3] = "structure_id"
#	requested_gene_expression[,3] = paste("Allen:",requested_gene_expression[,3],sep="")
	remember = list(rge=requested_gene_expression, test=c(test,dataset), genes=gene_values)
	# save to package environment
	aba_env = as.environment("package:ABAEnrichment")
	unlock_environment(aba_env)
	if (exists("remember",where=aba_env)){
		unlockBinding("remember",aba_env)
	}	
	aba_env$remember = remember
	lockEnvironment(aba_env, bindings=TRUE)

	final_output = list(results=results, genes=gene_values, cutoffs=cutoffs)
	
	if (!silent) message("\nDone.\n")	
	return(final_output)
}

	
