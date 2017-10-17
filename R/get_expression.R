# reads list with expression info from last aba_enrich run, mostly to avoid recomputing of aggregate() which can be time consuming, i.e. for a large set of background genes
# alternatively genes and dataset can be specified and data are loaded (if necessary) and processed again
# returns a list with a data_frame for each age category that is ready to heatmap()-plot

# remember$rge: expression data in long format
# remember$test: the test and dataset used 
# remember$genes

# change from long to wide format
to_wide=function(suby){
	suby=suby[order(suby$gene_id,as.character(suby$structure_id)),]	
	# matrix is filled column-wise, expression data have genes as blocks, i.e. genes=columns
	tab=matrix(nrow=length(unique(suby$structure_id)), ncol=length(unique(suby$gene_id)),data=suby$signal)
	rownames(tab)=unique(suby$structure_id)
	colnames(tab)=unique(suby$gene_id)
	return(tab)
}	

get_expression=function(structure_ids, gene_ids=NA, dataset=NA, background=FALSE)
{	
	if(xor(is.na(gene_ids)[1], is.na(dataset))) stop("Please specifiy gene_ids and dataset. Alternatively specify none of them which will use data from last aba_enrich() execution.")
	
	# A) use precomputed data from last FUNC run
	if(is.na(dataset)){
		new_data=FALSE
		aba_env=as.environment("package:ABAEnrichment")
		if(!exists("remember",where=aba_env)) {
			stop("This function requires previous execution of aba_enrich(). Alternatively specify gene_ids and dataset")
		}
		expr=aba_env$remember$rge
		genes=aba_env$remember$genes
		dataset=aba_env$remember$test[2]
		# delete background genes if not wanted
		if(aba_env$remember$test[1]=="hyper" & !background){
			expr=expr[expr$gene_id %in% genes[genes[,2]==1, 1] ,]
		}
	# B) take user defined genes and dataset and load 'raw' pre_input_data
	} else {
		new_data=TRUE
		valid_datasets=c("adult","5_stages","dev_effect")
		if(!(dataset %in% valid_datasets)){
			stop("Not a valid dataset. Please use 'adult', '5_stages' or 'dev_effect'.")	
		}
		expr=load_by_name(paste("dataset",dataset,sep="_"))
		expr=expr[,c("age_category",detect_identifier(gene_ids[1]),"structure","signal")]
		# warn if any genes do not have data 		
		not_in=gene_ids[!(gene_ids %in% expr[,2])]
		if(length(not_in)>0){
			not_in_string=paste(not_in,collapse=", ")
			warning(paste("No expression data for genes: ",not_in_string,".",sep=""))
		}
		# restrict to required genes and aggregate
		expr=expr[expr[,2] %in% gene_ids,]
		expr=aggregate(expr$signal,by=list("age_category"=expr[,1],"gene_id"=expr[,2],"structure_id"=expr[,3]),mean)
		colnames(expr)[4]="signal"
		expr$structure_id=paste("Allen:",expr$structure_id,sep="")
	}	
		
	# get all structures that provide data to the structures
	daty_children=c()
	structure_ids = as.character(structure_ids)
	for (i in 1:length(structure_ids)){
		daty_children=c(daty_children, get_sampled_substructures(structure_ids[i]))}	
	# subset to queried structures (the structures that provide data)	
	expr=expr[expr$structure_id %in% daty_children,]
	if(nrow(expr)==0) stop("No data for those genes-structure-dataset-combinations.")
	# rearrange data to wide format and one df per age_category (for easy heatmap()-calling)	
	expr_list=by(expr, expr$age_category, to_wide, simplify=FALSE)	
	
	# SORT GENES
	# A) data from aba_enrich()
	if(!(new_data)){
		# cluster background and test genes / sort wilcox-scores (no sense for hyper test-genes-only)
		if(!(aba_env$remember$test[1]=="hyper" & !background)){		
			# exclude duplicates for sorting
			genes=unique(genes)		
			# cant order with genes vector if only contains test-genes and bg is requested
			if(aba_env$remember$test[1]=="hyper" & sum(genes[,2])==nrow(genes) & background){
				expr_list=lapply(expr_list,function(x) { return(cbind(x[,colnames(x) %in% genes[,1]], x[,!(colnames(x) %in% genes[,1])]))})
			} else {
				# order with genes vector
				genes=genes[genes[,1] %in% colnames(expr_list[[1]]),]
				genes=genes[order(genes[,2],genes[,1]),]
				expr_list=lapply(expr_list,function(x) return(x[,match(genes[,1],colnames(x)),drop=FALSE]))
			}	
		}
	# B) user-defined genes and dataset	(sort according to order in which requested)
	} else {
		gene_ids=unique(gene_ids)
		gene_ids=gene_ids[gene_ids %in% colnames(expr_list[[1]])]
		expr_list=lapply(expr_list,function(x) return(x[,match(gene_ids,colnames(x)),drop=FALSE]))
	}
		
	# SORT REGIONS according to how they were requested
	daty_children=unique(daty_children)
	expr_list=lapply(expr_list,function(x) return(x[match(daty_children,rownames(x)),,drop=FALSE]))
	# if adult or dev_effect, return only data.frame; else give names to list items
	if(length(expr_list)==1){
		expr_list=expr_list[[1]]
	} else {
		names(expr_list)=paste("age_category",names(expr_list),sep="_")
	}
	# inform about dataset and units
	if(dataset=="adult") {
		units="z-score-normalized microarray expression data"
	} else if (dataset=="5_stages"){
		units="RPKM from RNA-seq"
	} else if (dataset=="dev_effect"){
		units="developmental effect score"
	}
	message(paste("Returning data from: ", dataset," (",units,").", sep=""))
	return(expr_list)
}	



