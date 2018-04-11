## plot expression

# use get_expression to retrieve data
# take data either from last aba_func()-call or specify genes and dataset explicitly


plot_expression=function(structure_ids, gene_ids=NA, dataset=NA, background=FALSE, dendro=TRUE, age_category=1){

    # take data from last aba_enrich()-call or retrieve them from pre_input?
    if(is.na(dataset)){
        aba_env=as.environment("package:ABAEnrichment")
        if(!exists("remember",where=aba_env)) {
            stop("This function requires previous execution of aba_enrich(). Alternatively define dataset and gene_ids.")
        }
        dataset_type=aba_env$remember$test[2]
        test=aba_env$remember$test[1]
        genes=aba_env$remember$genes
    } else {
        dataset_type=dataset
        test=NULL
    }   
    
    # check that at least 2 genes are requested (one structure is ok, it may have several sampled substructures)
    if(!(is.na(gene_ids)) && length(gene_ids)<2){
        stop("Please enter at least two genes.")
    }
    
    # check age_category (remaining checks performed in get_expression())
    age_category = as.numeric(age_category)
    if (dataset_type=="5_stages" & (length(age_category)!=1 || !(age_category %in% 1:5))){
        stop("Please specify an age_category between 1 and 5.")
    }
    
    # get expression data
    expr=get_expression(structure_ids, gene_ids, dataset, background)
    if (dataset_type=="5_stages") {
        expr=expr[[age_category]]
        main_append=paste("age_category: ",age_category,sep="")
    } else { main_append="" }
    
    # ensure there are at least 2 columns and rows (heatmap won't accept less)
    if(nrow(expr)<2){
        stop("At least two brain structures are needed for plotting.")
    }
    if(ncol(expr)<2){
        stop("At least two genes are needed for plotting.")
    }
    
    # add acronym
    full_names=c()
    for (id in rownames(expr)){ 
        full_names=c(full_names ,get_name(id))
    }
    acronyms=sapply(strsplit(full_names,"_"),"[[",1)
    rownames(expr)=paste(rownames(expr), " (",acronyms,")",sep="")  
    
    # create plot   
    cexRow = min(1.5,0.2 + 1/log10(nrow(expr)))
    cexCol = min(1.5,0.2 + 1/log10(ncol(expr)))
#   colramp=colorRamps::matlab.like(100)
    colramp=rev(heat.colors(100))
    
    if (dendro==TRUE){
        gplots::heatmap.2(expr,scale="none",col=colramp,margins=c(13.5,15), main=paste(test,dataset_type,main_append,sep=" "),density.info="none",trace="none",keysize=1.2,cexRow=cexRow,cexCol=cexCol)
    } else {            
        ## define colors of sidebar
        if (is.null(test)){
            sidebar=rep("white",ncol(expr)) # none for no enrichment-test input
        } else if (test=="hyper"){
            coly=c("black","red")
            sidebar=coly[genes[match(colnames(expr),genes[,1]), 2] + 1]
        } else { # binom, conti, wilcox
            # for wilcox plot scores, for the others combined value
            if (test == "binomial"){ # A/(A+B) for binomial
                genes[,2] = genes[,2] / (genes[,2] + genes[,3])
            } else if (test == "contingency"){ # (A/B)/(C/D) for binomial, (add 1 to prevent division by 0)
                genes[,2] = (genes[,2]+1 / genes[,3]+1) / (genes[,4]+1 / genes[,5]+1)
            }
            coly=rev(rainbow(50,start=0,end=0.5))
            genes[,2]=genes[,2]-min(genes[,2])
            genes[,2]=round(genes[,2]/ max(genes[,2]) * 49) 
            sidebar=coly[genes[match(colnames(expr),genes[,1]), 2] + 1] 
        }
        gplots::heatmap.2(expr,scale="none",col=colramp,margins=c(13.5,15),Colv=NA,Rowv=NA,dendrogram="none",ColSideColors=sidebar, main=paste(test,dataset_type,main_append,sep=" "),density.info="none",trace="none",keysize=1.2,cexRow=cexRow,cexCol=cexCol)
    }
}



