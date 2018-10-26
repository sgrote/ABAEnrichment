
## plot expression given matrix as returned by get_expression

plot_expression=function(expr_mat, dendro=TRUE, gene_vars=NULL, main=""){
    
    # check arguments    
    if (!(is.matrix(expr_mat))){
        stop("expr_mat needs to be a matrix.")
    }
    # ensure there are at least 2 columns and rows (heatmap won't accept less)
    if(nrow(expr_mat)<2){
        stop("At least two brain structures are needed for plotting.")
    }
    if(ncol(expr_mat)<2){
        stop("At least two genes are needed for plotting.")
    }
    
    # add acronym
    full_names = get_name(rownames(expr_mat))
    acronyms=sapply(strsplit(full_names,"_"),"[[",1)
    rownames(expr_mat)=paste(rownames(expr_mat), " (",acronyms,")",sep="")  
    
    # create plot   
    cexRow = min(1.5,0.2 + 1/log10(nrow(expr_mat)))
    cexCol = min(1.5,0.2 + 1/log10(ncol(expr_mat)))
#   colramp=colorRamps::matlab.like(100)
    colramp=rev(heat.colors(100))
    
    if (dendro){
        Colv = TRUE
        Rowv = TRUE
        dendrogram = "both"
    } else {
        Colv = NA
        Rowv = NA
        dendrogram = "none"
    }       
    if (is.null(gene_vars)){
        gplots::heatmap.2(expr_mat, scale="none", col=colramp, margins=c(13.5,15), Colv=Colv, Rowv=Rowv, dendrogram=dendrogram, main=main, density.info="none", trace="none", keysize=1.2,cexRow=cexRow, cexCol=cexCol)
    } else {
        # extract needed gene_vars
        gene_vars = gene_vars[match(colnames(expr_mat),gene_vars[,1]),]
        if (any(is.na(gene_vars))){
            stop("Not all genes are present in gene_vars[,1], or gene_vars contains NA.")
        }
        if (ncol(gene_vars) == 2 && all(gene_vars[,2] %in% c(0,1))){
            coly=c("black","red")
            leg_main = "hyper"
        } else { # binom, conti, wilcox
            coly=rev(rainbow(50,start=0,end=0.5))
            # for wilcox plot scores, for the others combined value
            if (ncol(gene_vars) == 3){ # A/(A+B) for binomial
                leg_main = "A/(A+B)\n"  # \n to move title up and use small y.intersp in legend
                gene_vars[,2] = gene_vars[,2] / (gene_vars[,2] + gene_vars[,3])
            } else if (ncol(gene_vars) == 5){ # (A/B)/(C/D) for conti, (add 1 to prevent division by 0)
                leg_main = "(A+1/B+1)/(C+1/D+1)\n"
                gene_vars[,2] = ((gene_vars[,2]+1) / (gene_vars[,3]+1)) / ((gene_vars[,4]+1) / (gene_vars[,5]+1))
            } else {
                leg_main = "gene-score\n"
            }
            gene_vars[,2]=gene_vars[,2]-min(gene_vars[,2])
            gene_vars[,2]=round(gene_vars[,2]/ max(gene_vars[,2]) * 49)
        }
        sidebar=coly[gene_vars[,2] + 1]
        # reorder based on gene_vars if dendrogram is omitted
        if (!dendro){
            expr_mat = expr_mat[,order(gene_vars[,2])]
            sidebar = sidebar[order(gene_vars[,2])]
        }
        gplots::heatmap.2(expr_mat, scale="none", col=colramp, margins=c(13.5,15), Colv=Colv, Rowv=Rowv, dendrogram=dendrogram, ColSideColors=sidebar, main=main, density.info="none", trace="none", keysize=1.2,cexRow=cexRow, cexCol=cexCol)
        
        # legend for gene-vars
        lc = 0.8
        if (leg_main=="hyper"){
            legend("topright", legend=c("candidate","background"), fill=c("red","black"), bty="n", cex=lc)
        } else {
            ncol = 10
            icol = round(seq(50,1,length.out=ncol))
            legend("topright", legend=c("high",rep("",ncol-2),"low"), fill=c(coly[icol]), bty="n",          border=c(coly[icol]), cex=lc, title=leg_main, y.intersp=0.4, xpd=TRUE)
        }           
    }    
}



