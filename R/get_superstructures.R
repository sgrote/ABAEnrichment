
# get all structures in the hierarchy from a given structure to the root (brain)

get_superstructures=function(structure_id){
	if(length(structure_id)!=1){
		stop("Please use a single structure id.")
	}
	# remove Allen: string
	struc=strsplit(structure_id,":")[[1]][2]
	# check if valid structure_id
	if(!(struc) %in% daty_strucs[,1]){
		stop(paste("No data for structure_id ",structure_id,".",sep=""))
	}
	# get path
	path=as.character(daty_strucs[daty_strucs[,1]==struc,2])
	path=unlist(strsplit(path,"/"))
	# remove empty (first) entry
	path=path[path!=""]
	# add Allen: string again
	path=paste("Allen:",path,sep="")
#	# turn around, to go frome from current node to root
#	path=rev(path)
	return(path)
}


