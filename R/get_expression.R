# get matrix of expression data
# given structures, genes and dataset
# returns a data_frame for each age category that is ready to heatmap()-plot


# change from long to wide format
to_wide=function(suby){
    suby=suby[order(suby$gene_id,as.character(suby$structure_id)),] 
    # matrix is filled column-wise, expression data have genes as blocks, i.e. genes=columns
    tab=matrix(nrow=length(unique(suby$structure_id)), ncol=length(unique(suby$gene_id)),data=suby$signal)
    rownames(tab)=unique(suby$structure_id)
    colnames(tab)=unique(suby$gene_id)
    return(tab)
}   

get_expression=function(structure_ids, gene_ids, dataset="adult")
{   
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
    # sort according to order in which requested
    gene_ids=unique(gene_ids)
    gene_ids=gene_ids[gene_ids %in% colnames(expr_list[[1]])]
    expr_list=lapply(expr_list,function(x) return(x[,match(gene_ids,colnames(x)),drop=FALSE]))
        
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



