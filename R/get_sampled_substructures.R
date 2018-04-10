
# get all structure that provide data for a given structure

get_sampled_substructures=function(structure_id){
    if(length(structure_id)!=1){
        stop("Please use a single structure id.")
    }
    structure_id = as.character(structure_id)
    # subset ontology to directly sampled leaves and their ancestors
    directly_sampled=daty_strucs[daty_strucs$directly_sampled==1,]
    # remove Allen: string if present
    if(grepl("Allen:", structure_id)){
        struc=strsplit(structure_id,":")[[1]][2]
    } else {
        struc=structure_id  
    }   
    # grep structure_id, avoid getting "4008" when required "400"
    daty_children=directly_sampled[grep(paste("/",struc,"/",sep=""),directly_sampled[,2]),1]
    if(length(daty_children)==0) stop(paste("No data for structure_id ",structure_id,".",sep=""))
    # re-add Allen:-prefix if it was in input
    if(grepl("Allen:", structure_id)){
        daty_children=paste("Allen:",daty_children,sep="")  
    }   
    return(as.character(daty_children))
}
