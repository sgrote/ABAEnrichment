
# leave only FWER
# delete lines with structures that have same info as inferior nodes
# add string with brain region names for regions of the same info cluster
# compute min and mean FWER per info cluster
# put FWER for different cutoffs per structure and age into one field

rearrange_output=function(summary,cluster,term){

    # add region name
    summary$node_name=term[match(summary[,4],term$acc),"name"]
    # add category 
    summary$cate=cluster[match(summary[,4], cluster$samples),"lables"]
    # confirm that inside clusters FWERs are the same
    # xx=aggregate(summary[,8],by=list(summary$cate,summary$age_category,summary$cutoff_quantile), function(x) length(unique(x)))
    # if(!all(xx$x==1)) stop("Brain regions grouped have different FWERs.")

    # exclude redundant structures and stay only with FWER_overrep
    pure=unique(summary[,c(1:3,6,8,10)])
    # how often is FWER < 0.05
    pure$age_category=as.factor(pure$age_category)
    pure$cate=as.factor(pure$cate)
    freqs=as.data.frame(table(pure[pure[,5] < 0.05,c("age_category","cate")]))
    colnames(freqs)[3]="n_significant"
    # lowest FWER
    lowest=aggregate(pure[,5], by=list(age_category=pure$age_category,cate=pure$cate),min)
    colnames(lowest)[3]="min_FWER"
    # mean FWER
    meany=aggregate(pure[,5], by=list(age_category=pure$age_category,cate=pure$cate),mean)
    colnames(meany)[3]="mean_FWER"  
    
    ## collapse FWER over cutoffs
    fwers=aggregate(1:nrow(pure), by=list(cate=pure$cate,age_category=pure$age_category), function(x) paste(pure[x,5][order(pure[x,"cutoff_quantile"])],collapse=";"))
    colnames(fwers)[3]="FWERs"
    
#   ## NEW: collapse p_vals over cutoffs (erstmal nicht)
#   pvals=aggregate(1:nrow(pure), by=list(cate=pure$cate,age_category=pure$age_category), function(x) paste(pure[x,4][order(pure[x,"cutoff_quantile"])],collapse=";"))
#   colnames(pvals)[3]="raw_p_vals"
    # TODO: p-vals can be quite long to be pasted in one column, maybe use signif(x, digits = 3) rounding to significant decimal places

    ## add flagship name
    cluster$region=summary[match(cluster$flagship,summary[,4]),"node_name"]
    
    # collapse structure_names for categories
    #cluster$name=summary[match(cluster$samples,summary[,4]),"node_name"]
        
    names_string=aggregate(1:nrow(cluster),by=list(cate=cluster$lables),function(x) paste(cluster[x,"samples"][order(cluster[x,"samples"],decreasing=TRUE)],collapse=";"))
    
    colnames(names_string)[2]="equivalent_structures"
        
    ## merge together
    preout=data.frame(cluster[match(freqs$cate,cluster$lables),3:4],freqs)
    preout=merge(preout,meany)
    preout=merge(preout,lowest)
    preout=merge(preout,names_string)
    preout=merge(preout,fwers)
#   preout=merge(preout,pvals)
    colnames(preout)[3:4]=c("structure_id","structure")
#   # preorder with mixedorder on structure name first (mixedorder does not take multiple columns)
#   preout=preout[mixedorder(preout$structure),]
    # order
    out=preout[order(-1*(preout$n_significant),preout$min_FWER,preout$mean_FWER,preout$age_category,preout$structure_id),-1]
    # order rownames 
    rownames(out)=1:nrow(out)
                    
    return(out)
}
