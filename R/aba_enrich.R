# run "FUNC" with gene expression data from Allen Brain Atlas 
# over flexible expression cutoff
# wilcoxon rank test or hypergeometric test
# with adult microarray data, developmental RNA-seq (5 age categories) or derived developmental effect score


###############
# MAIN function
###############

### Arguments
# genes: vector of 0/1 (hypergeometric) or float (wilcoxon) with gene identifiers as names (ensembl, entrez or hgnc)
# dataset: "adult", "5_stages" or "dev_effect"
# test: "hyper" or "wilcoxon"
# cutoff_quantiles: numeric values in ]0,1[
# n_randsets: number of random-sets

#######################################################################################
# 1. check arguments and define parameters
# 2. load expression data
# 3. loop over age categories and cutoffs
#	3a) create input
#	3b) run func
#	3c) summarize func output
# 4. create output
#######################################################################################

aba_enrich=function(genes,dataset="adult",test="hyper",cutoff_quantiles=seq(0.1,0.9,0.1),n_randsets=1000)
{

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
	valid_datasets=c("adult","5_stages","dev_effect")
	if(!(dataset %in% valid_datasets)){
		stop("Not a valid dataset. Please use 'adult', '5_stages' or 'dev_effect'.")	
	}
	if(length(n_randsets)>1 || !is.numeric(n_randsets) || n_randsets<1){
		stop("Please define 'n_randsets' as a positive integer.")
	}
	if(n_randsets != round(n_randsets)){
		n_randsets=round(n_randsets)
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
	
	## are all genes in data? if not, warn, list missing genes and ask to continue
	# load gene_list and get gene identifier
	gene_symbols=get(paste("gene_symbols",folder_ext,sep="_"))	
	identifier=detect_identifier(names(genes)[1])	
	gene_list=gene_symbols[,identifier]	
	not_in=names(genes)[!(names(genes) %in% gene_list)]
	remaining_genes=genes[names(genes) %in% gene_list]
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
		not_in_string=paste(not_in,collapse=", ")
#		message(paste("Warning: no expression data for genes: ",not_in_string,sep=""))
#		repeat{		
#			answer=readline("Continue anyway? Y/n: ")
#			if(answer %in% c("n","N")) return(message("Analysis stopped."))
#			if(answer %in% c("y","Y","")) break
#		}
	}	
	# Create input and output folder and write file with missing genes
	directory=tempdir()

	#####	2. Load expression data
	
	# load pre input
	message("Loading dataset...")
	pre_input=load_by_name(paste("dataset",dataset,sep="_"))
	
	# (filter dev_effect_score by max fraction of zero expression)
	#if(dataset=="dev_effect") pre_input=pre_input[(pre_input$zeros/pre_input$total)<0.5,]	
	
	# select gene identifier
	pre_input = pre_input[,c("age_category",identifier,"structure","signal")]
	colnames(pre_input)[2]="gene_id"

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
	# number of structure-age-combinations: if all structures are present in all ages	they can simply be multiplied -> thats the case
	# n=nrow(unique(pre_input[,c("age_category","structure")]))
	n=length(unique(pre_input$age_category))*length(unique(pre_input$structure))
	ngenes=table(pre_input$gene_id)
	bose=ngenes[ngenes>n]
		if(length(bose)>0){
		message(" Aggregate expression per gene...")
		part1=pre_input[pre_input$gene_id %in% names(bose),]
		part2=pre_input[!(pre_input$gene_id %in% names(bose)),]
		part1=aggregate(part1$signal,by=list("age_category"=part1$age_category,"gene_id"=part1$gene_id,"structure"=part1$structure),mean)
		colnames(part1)[4]="signal"
		pre_input=rbind(part1,part2)
		message(" Done.")
	} else {message(" No gene duplicates - no aggregation needed.")}		
	
	# compute cutoffs (tapply lasts much longer - use only when needed)
	message("Computing cutoffs... ")
	if(dataset=="5_stages"){
		cutoff_list = tapply(pre_input$signal, pre_input$age_category, function(x) quantile(x, probs=cutoff_quantiles),simplify=FALSE)
	} else {	
		cutoff_list = list()
		cutoff_list[[1]] = quantile(pre_input$signal,probs=cutoff_quantiles)
		names(cutoff_list) = pre_input$age_category[1]
	}
	
	# load ontology and write files to tmp
	term=get(paste("term",folder_ext,sep="_"))
	term2term=get(paste("term2term",folder_ext,sep="_"))
	graph_path=get(paste("graph_path",folder_ext,sep="_"))
	write.table(term,file=paste(directory, "/term.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	write.table(term2term,file=paste(directory, "/term2term.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	write.table(graph_path,file=paste(directory, "/graph_path.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

	#####	3. loop over age categories and cutoffs
	# initialize output
	out = data.frame()
	for (s in unique(pre_input$age_category)){
		if (dataset=="5_stages"){
			stage_input = pre_input[pre_input$age_category==s,2:4]	
		} else 	stage_input = pre_input[,2:4]
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
			if (test=="hyper"){		
				input$signal = 0
				interesting_genes = names(remaining_genes)[remaining_genes==1]
				input[input[,1] %in% interesting_genes,"signal"] = 1
				# do input Input data have both 0 and 1? else break (FUNC would throw error and summary-file would not be generated)
				if(sum(input$signal)==0){
					message(paste("  Warning: No statistics for cutoff >= ",cutoff_quantiles[i],". No test gene expression above cutoff.",sep=""))
				break
				}	
				if(sum(input$signal)==nrow(input)){
					message(paste("  Warning: No statistics for cutoff >= ",cutoff_quantiles[i],". No background gene expression above cutoff.",sep=""))
				break
				}	
			# for Wilcoxon test: replace expression signal with score
			} else	if (test=="wilcoxon"){
				if(length(unique(input[,1]))<2){
					message(paste("  Warning: No statistics for cutoff >= ",cutoff_quantiles[i],". Less than 2 genes show expression above cutoff.",sep=""))
					break
				}
				score = genes[match(input$gene_id,names(genes))]
				input$signal = score
			}	
			message(" Rearrange input for FUNC...")
			# add "Allen:" string to structures (shouldn't be numeric)
			input$structure=paste("Allen:",input[,2],sep="")	
			# generate input and output name (cutoff and age_category)
			input_name=paste("cutoff",names(cutoff)[i],s,sep="_")	
			# write FUNC-input to Input-directory (how it would be for FUNC standalone)		
			write.table(input,paste(directory,"/",input_name,".tsv",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")	
			
			# prepare input data (infile-data and root like in separate_taxonomies.pl)
			# paste annotations
			xx=tapply(input[,2],input[,1],function(x){paste(x,collapse=" ")})
			# create "root" dataframe
			root=data.frame(genes=as.character(names(xx)),goterms=as.character(xx))
			# "infile-data" (wilcox) / "Allen:4005-changed" (hyper)
			if (test=="hyper"){
				# subset to test genes
				input=input[input[,3]==1,]	
				# change factor to character			
				input[,1]=as.character(input[,1])
				input[,2]=as.character(input[,2])
				# paste test annotations
				yy=tapply(input[,2],input[,1],function(x){paste(paste(x,collapse=" "),"")})
				# create infile-data dataframe
				infile_data=data.frame(genes=as.character(names(yy)),goterms=as.character(yy))
			} else if(test=="wilcoxon"){
				infile_data = unique(input[,c(1,3)])
			}
	
			# write input-files to tmp-directory
			write.table(infile_data,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste(directory,"/infile-data",sep=""))
			write.table(root,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste(directory,"/",root_node,sep=""))
								
			##	3b) run func
			if(test=="hyper"){
				run_func(hyper_randset, hyper_category_test, directory, root_node, n_randsets)			
			} else if (test=="wilcoxon"){				
				run_func(wilcox_randset, wilcox_category_test, directory, root_node, n_randsets)
			}	
		
			## 3c) summarize func-output
			groupy=read.table(paste(directory,"/category_test_out",sep=""))
			age_category=rep(s,nrow(groupy))			
			cutoff_quantile=rep(cutoff_quantiles[i],nrow(groupy))
			cutoff_value=rep(cutoff[i],nrow(groupy))
			groups=cbind(age_category,cutoff_quantile,cutoff_value,groupy)
				
			# combine with output from previous cutoffs
			out = rbind(out,groups)			
		}		
	}
	
	#####. 4 rearrange output
	message("Creating output...")
	
	if (test == "hyper"){
		colnames(out)[4:8]=c("node_id","raw_p_underrep","raw_p_overrep","FWER_underrep","FWER_overrep")
	} else if (test == "wilcoxon"){
		colnames(out)[4:8]=c("node_id","raw_p_low_rank","raw_p_high_rank","FWER_low_rank","FWER_high_rank")
	}
	# remove redundant rows (nodes with no data and only one child)
	cluster=get(paste("node_clusters",folder_ext,sep="_"))
	results=rearrange_output(out,cluster,term)
	# rearrange cutoff-list
	cutoffs=do.call(cbind.data.frame, cutoff_list)
	colnames(cutoffs)=paste("age_category",colnames(cutoffs),sep="_")

	# save in package environment: dataframes for test and background-gene expression, test, dataset	
	# (to be returned with get_expression(structure_ids))	
	requested_gene_expression = pre_input
	rownames(requested_gene_expression)=NULL
	colnames(requested_gene_expression)[3]="structure_id"
	requested_gene_expression[,3]=paste("Allen:",requested_gene_expression[,3],sep="")
	remember=list(rge=requested_gene_expression,test=c(test,dataset),genes=genes)
	# save to package environment
	aba_env=as.environment("package:ABAEnrichment")
	unlock_environment(aba_env)
	if (exists("remember",where=aba_env)){
		unlockBinding("remember",aba_env)
	}	
	aba_env$remember=remember
	lockEnvironment(aba_env, bindings = TRUE)

	final_output=list(results=results,genes=remaining_genes,cutoffs=cutoffs)
	message("\nDone.\n")
	if (length(not_in)>0){
		warning(paste("No expression data for genes: ",not_in_string,".\n  These genes were not included in the analysis.",sep=""))
	}		
	return(final_output)	
}

	
