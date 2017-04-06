
## return the names of structures given the structure_id

get_one_name=function(structure_id, no_allen){
	if(!(structure_id %in% term_developmental$acc | structure_id %in% term_adult$acc)){
		if(no_allen){
			structure_id = strsplit(structure_id,":")[[1]][2]
		}
		stop(paste("Invalid structure_id: ",structure_id,".",sep=""))
	}	
	xx=rbind(term_developmental,term_adult)
	xx=xx[xx$is_relation==0,]	
	xx$name=as.character(xx$name)
	name=xx[xx$acc==structure_id, "name"]
	return(name)
}

get_name=function(structure_ids){
	structure_ids = as.character(structure_ids)
	# add and then remove again Allen: string if not present
	no_allen = FALSE
	if(!(grepl("Allen:", structure_ids[1]))){
		structure_ids = paste("Allen:",structure_ids,sep="")
		no_allen=TRUE
	}	
	out = sapply(structure_ids,get_one_name,no_allen)
	if(no_allen){
		names(out) = sapply(names(out), function(name) {strsplit(name,":")[[1]][2]})
	}
	return(out)
}
