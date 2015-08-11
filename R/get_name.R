
## return the names of structures given the structure_id

get_one_name=function(structure_id){
	if(!(structure_id %in% term_developmental$acc | structure_id %in% term_adult$acc)){
		stop(paste("Invalid structure_id: ",structure_id,".",sep=""))
	}	
	xx=rbind(term_developmental,term_adult)
	xx=xx[xx$is_relation==0,]	
	xx$name=as.character(xx$name)
	name=xx[xx$acc==structure_id, "name"]
	return(name)
}

get_name=function(structure_ids){
	return(sapply(structure_ids,get_one_name))
}
