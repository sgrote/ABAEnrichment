
# load dataset from data package into package environment 
# to prevent problems with existing objects in the global environment
# load only if not already loaded 
# return it so that it can be renamend

load_by_name = function(data_name, silent=FALSE){
	aba_env = as.environment("package:ABAEnrichment")
	already_loaded = exists(data_name, where=aba_env)
	if(!silent) message(paste(data_name,"already exists in package environment:", already_loaded))
	if(!already_loaded){
		if(!silent) message(paste(" Load ",data_name,"...", sep=""))
		unlock_environment(aba_env)
		data(list=data_name, package="ABAData" ,envir=aba_env)	
		lockEnvironment(aba_env, bindings = TRUE)
		if(!silent) message(" Done.")
	}	
	return(get(data_name, envir=aba_env))
}	
				
	
