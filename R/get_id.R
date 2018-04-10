
## return the names of structure_ids of and names of brain regions that match ONE name

get_id=function(structure_name){
    if(length(structure_name)!=1){
        stop("Please use a single brain region name.")
    }
    terms = rbind(term_developmental, term_adult)
    # restrict to brain regions that have expression data
    terms = terms[terms$acc %in% paste("Allen:", daty_strucs$id,sep=""), ]
    terms = terms[terms$is_relation==0,]    
    terms$name = as.character(terms$name)
    terms$term_type = as.character(terms$term_type)
    
    # grep for name in combined adult and developmental term-tables
    out = terms[grepl(structure_name, terms$name, ignore.case=TRUE),2:4]
    if(nrow(out)==0){
        stop("No matches found.")
    }
    rownames(out) = 1:nrow(out)
    colnames(out) = c("structure", "ontology", "structure_id")
    out[out[,"ontology"]=="neural plate","ontology"] = "developmental"  
    out[out[,"ontology"]=="Brain","ontology"] = "adult"  
    
    return(out)
}

