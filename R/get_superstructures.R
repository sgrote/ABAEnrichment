
# get all structures in the hierarchy from a given structure to the root (brain)

get_superstructures=function(structure_id){
	if(length(structure_id)!=1){
		stop("Please use a single structure id.")
	}
	structure_id = as.character(structure_id)
	# remove Allen: string if present
	if(grepl("Allen:", structure_id)){
		struc=strsplit(structure_id,":")[[1]][2]
	} else {
		struc=structure_id	
	}	
	# check if valid structure_id
	if(!(struc %in% daty_strucs[,1])){
		stop(paste("No data for structure_id ",structure_id,".",sep=""))
	}
	# get path
	path=as.character(daty_strucs[daty_strucs[,1]==struc,2])
	path=unlist(strsplit(path,"/"))
	# remove empty (first) entry
	path=path[path!=""]
	# add Allen: string again if it was in input	
	if(grepl("Allen:", structure_id)){
		path=paste("Allen:",path,sep="")
	}
#	# turn around, to go frome from current node to root
#	path=rev(path)
	return(path)
}


