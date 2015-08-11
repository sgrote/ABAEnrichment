## run FUNC with either "hyper" or "wilcoxon"
#(input are files except for n_randsets)
#(output is file)

run_func=function(FUN1, FUN2, directory, root_node, n_randsets)
{
	# randset
	FUN1(
		paste(directory,"/",root_node,sep=""),
		paste(directory,"/infile-data",sep=""),
		n_randsets,
		paste(directory, "/randset_out",sep=""),
		paste(directory, "/term.txt",sep=""), 
		paste(directory, "/graph_path.txt",sep=""),
		paste(directory, "/term2term.txt",sep=""),
		root_node
	)	
	# category test			
	FUN2(
		paste(directory, "/randset_out",sep=""),
		paste(directory,"/category_test_out",sep=""),
		1,
		root_node
	)	
}	
	

