
# return a dataframe with genes annotated to (enriched) brain regions 

# input: result list produced by aba_enrich, optional: FWER, background
# alternative input: structure_ids, cutoff-quantiles, dataset, (optional: genes)

# output: age, strcuture_id, cutoff_quantile, gene, (FWER, candiate/score)
# additional output: structure name? (sampled substructure? - vielleicht nur verwirrend)


get_annotated_genes = function(res, fwer_threshold=0.05, background=FALSE, structure_ids=NULL, dataset=NULL, cutoff_quantiles=NULL, genes=NULL){
    
    ## given res as input:
    if(!(missing(res))){
        ## checks
        # check that it is really a res-object
        if(!is.list(res) || is.null(names(res)) || !(all(names(res) == c("results","genes","cutoffs")))){
            stop("Please use an object returned from aba_enrich as input (list with 3 elements).\n Alternatively structure_ids, cutoff_quantiles and dataset need to be defined.")
        }
        # outcommented to also allow a threshold > 1 returning all annotations (especi. including root node)
#       if(fwer_threshold > 1 | fwer_threshold <= 0){
#           stop("Please set a valid fwer_threshold (0 < fwer_threshold <= 1).")
#       }
        if(!is.null(structure_ids)){
            warning("Unused argument: structure_ids")
        }       
        if(!is.null(cutoff_quantiles)){
            warning("Unused argument: cutoff_quantiles")
        }       
        if(!is.null(dataset)){
            warning("Unused argument: dataset")
        }
        if(!is.null(genes)){
            warning("Unused argument: genes")
        }
        ## get significant age-structure-cutoff combinations given the fwer_threshold
        message("extract enriched brain regions...")
        combis = get_signi_combis(res, fwer_threshold)
        # extract dataset (given the age_categories: 5_stages: 1-5, adult:5, dev_effect:0)
        if (length(unique(res[[1]][,1])) == 1){
            if (res[[1]][1,1] == 5) dataset = "adult"
            if (res[[1]][1,1] == 0) dataset = "dev_effect"
        } else if (length(unique(res[[1]][,1])) == 5) dataset = "5_stages"
    } else {
        ## checks
        # structures
        if(is.null(structure_ids)){
            stop("Please define structure_ids.")
        }
#       structure_ids = as.character(structure_ids)
        allen = TRUE
        if(!grepl("Allen:", structure_ids[1])){
            allen = FALSE
            structure_ids = paste("Allen:", structure_ids, sep="")
        }
        # cutoff_quantiles
        if(is.null(cutoff_quantiles)){
            stop("Please define cutoff_quantiles.")
        }
        if (class(cutoff_quantiles)!="numeric"){
            stop("Please enter numeric cutoff_quantiles.")
        }   
        if (min(cutoff_quantiles) < 0 | max(cutoff_quantiles) > 1){
            stop("Please enter cutoff_quantiles between 0 and 1.")
        }
        cutoff_quantiles = sort(unique(cutoff_quantiles))               
        # dataset
        if(is.null(dataset)){
            stop("Please define a dataset.")
        }   
        valid_datasets=c("adult","5_stages","dev_effect")
        if(!(dataset %in% valid_datasets)){
            stop("Not a valid dataset. Please use 'adult', '5_stages' or 'dev_effect'.")    
        }
        # background
        if(background){
            warning("Unused argument: background")
        }
    }
    
    # load expression dataset
    expr = load_by_name(paste("dataset",dataset,sep="_"))
    
    ## else given alternative input, make 'combis' dataframe like for standard input (but without FWER) 
    if(missing(res)){
        age_category = rep(unique(expr$age_category), length(structure_ids))
        structures = data.frame(structure_id = structure_ids, structure = get_name(structure_ids), row.names=NULL)
        structures = structures[rep(seq_len(nrow(structures)),each=length(unique(age_category))),]
        pre_combis = cbind(age_category, structures)
        # add cutoffs
        combis = pre_combis[rep(seq_len(nrow(pre_combis)),length(unique(cutoff_quantiles))),]
        combis$cutoff = rep(cutoff_quantiles, each=nrow(pre_combis))
        # get cutoff-values     
        message("compute cutoffs...")
        cutoff_list = tapply(expr$signal, expr$age_category, function(x) quantile(x, probs=cutoff_quantiles),simplify=FALSE)
        cutoff_df = do.call(cbind.data.frame, cutoff_list)
        combis$cutoff_value = apply(combis, 1, function(x) cutoff_df[paste(as.numeric(x[4])*100,"%",sep=""), as.character(x[1])]) # x[4]=quantile, x[1]=age
        combis$structure_id = as.character(combis$structure_id)
    }
    
    ## add annotated genes
    # for every structure_id, get IDs of sampled substructures
    message("get sampled substructures...")
    anno_children = sapply(unique(combis$structure_id), get_sampled_substructures, simplify=FALSE)
    # restrict expression data to those substructures
    expr$structure = paste("Allen:",expr$structure,sep="")
    expr = expr[expr$structure %in% unique(unlist(anno_children)), ]
    # restrict expression data to gene-identifier (and to input genes if present)
    if(missing(res)){
        if(is.null(genes)){
            expr = expr[,c(6,1,4,5)]
        } else {
            expr = expr[,c("age_category",detect_identifier(genes[1]),"structure","signal")]
            expr = expr[expr[,2] %in% genes,]
        }
    } else {
        # restrict to (candidate) genes
        genes = res[[2]]
        if(!background && all(genes[,2] %in% c(0,1))){
            genes = genes[genes[,2]==1,]
        }
        expr = expr[,c("age_category",detect_identifier(genes[1,1]),"structure","signal")]
        # don't restrict when only candidate defined for aba_enrich, but background=TRUE
        if(!(background && all(genes[,2] == 1))){
            expr = expr[expr[,2] %in% genes[,1],]
        }
    }
    if(nrow(expr)==0){
        warning("No expression data found for the input.")
        return(NULL)
    }
    colnames(expr)[2:3] = c("gene_id", "structure_id")
    
    # aggregate genes per structure and age (duplicates possible due to 3 different gene-IDs)
    message("aggregate gene expression...")
    expr = as.data.table(expr)
    expr_ag = expr[,mean(signal), by=list(age_category, gene_id, structure_id)]
    colnames(expr_ag)[4] = "signal"

    # subset expression to to signal>min(cutoff)
    expr_ag2 = expr_ag[signal > min(combis$cutoff_value),]

    # for every significant age-structure-cutoff combi:
    message("get annotated genes...")
    out = vector(mode = "list", length = nrow(combis)) # initialize list
    li = 1 # list index
    for (age in unique(combis$age_category)){
        expr_ag3 = expr_ag2[age_category == age,]
        combis_age = combis[combis$age_category == age,]
        for (i in 1:nrow(combis_age)){
            child_ids = anno_children[[combis_age[i,"structure_id"]]]
            # select genes above cutoff
            express = expr_ag3[structure_id %in% child_ids & signal > combis_age[i,"cutoff_value"],]
            # make unique
            expre_genes = unique(express$gene_id)
            # duplicate the line for every annotated gene
            to_out = combis_age[rep(i, length(expre_genes)), -c(3,5)] # remove structure-name and cutoff-value
            to_out$anno_gene = expre_genes
            out[[li]] = to_out
            li = li + 1
        }
    }
    out = do.call(rbind,out)

    if(nrow(out) == 0){
        warning("No GO-annotations found for the input genes.")
        return(out)
    }

    out = out[mixedorder(out$anno_gene),] # NEW: mixedorder genes instead of normal order
    ## given res as input:
    if(!(missing(res))){
        # add gene-associated variables
        out = data.frame(out, genes[match(as.character(out$anno_gene), genes[,1]), 2:ncol(genes)])
        if(ncol(out) == 6){ # hyper or wilcox
            # fix colname (lost through subset to one column)
            colnames(out)[6] = "score"
            # replace NA with 0 for background genes (needed for hyper and default background)
            out[is.na(out[,6]), 6] = 0
        }
        out = out[order(out$FWER, out$age_category, out$cutoff, out$structure_id, out[,6]),]
    } else {
        out = out[order(out$cutoff, out$age_category, out$structure_id),]
        if(!allen){
            out$structure_id = sapply(out$structure_id, function(name) {strsplit(name,":")[[1]][2]})
        }
    }
    # ugly row.names due to duplications
    row.names(out) = 1:nrow(out)
    
    return(out)
}


# significant age-structure-cutoff combinations (with cutoff-values)
get_signi_combis = function(res, fwer_threshold){
    fwers = res[[1]]
    cutoffs = res[[3]]
    # remove lines were min_FWER >= threshold
    fwers = fwers[fwers$min_FWER < fwer_threshold, ]
    if(nrow(fwers)==0){
        warning(paste("No significantly enriched brain regions at FWER-threshold ", fwer_threshold,sep=""))
        return(NULL)
    }
    # split FWERs
    # NEW erkenntnis: some structures have no expression at high cutoffs and are therefor missing in Func-output. LOESUNG: fill with FWER=1 -> will be below no cutoff 
    single_fwers = strsplit(fwers$FWERs, ";")
    # fill up rows with 1 that dont have the maximum number of FWERs (due to no expression at high cutoff)
    n_fwers = max(sapply(single_fwers, length))
    single_fwers = lapply(single_fwers, function(x){ 
        if(length(x) < n_fwers) x = c(x, rep("1", n_fwers-length(x)))
        return(x)
    })
    single_fwers = lapply(single_fwers, as.numeric)
    # restrict cutoffs to length of FWERs (remove unused (too high) cutoffs)
    used_cutoffs = row.names(cutoffs)[1:n_fwers]  # 10%, 20%, ...
    # for each line (age and structure), make a line for every significant cutoff
    out = fwers[,1:3]
    n_sig = unlist(lapply(single_fwers, function(x) sum(x < fwer_threshold)))
    # replicate each line the number of signi. cutoffs
    out = out[rep(seq_len(nrow(out)), n_sig),]
    # add cutoff-quantiles
    out$cutoff = unlist(lapply(single_fwers, function(x) used_cutoffs[x < fwer_threshold]))
    # add cutoff-values
    colnames(cutoffs) = substr(colnames(cutoffs), nchar(colnames(cutoffs)), nchar(colnames(cutoffs)))
    # (x[4]=quantile, x[1]=age)
    out$cutoff_value = apply(out, 1, function(x) cutoffs[x[4],as.character(x[1])]) 
    # transform cutoffs to [0,1]
    out$cutoff = as.numeric(unlist(strsplit(out$cutoff, "%"))) / 100    
    # add FWERs
    out$FWER = unlist(lapply(single_fwers, function(x) x[x < fwer_threshold]))
    
    return(out)
}
