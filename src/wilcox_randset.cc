
#include <time.h>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <memory>

#include "go_graph.h"
#include "idmap.h"
#include "transitions.h" 
#include "genes.h"

#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

#define MAX_LINE_LENGTH 20000

//[[Rcpp::export]]
void wilcox_randset(std::string nodes_per_gene ,int number_of_randomsets, std::string directory, std::string root) 
{
		
	/*****************
         * read graph-structure and create graph
	 *******************/
	// steffi: term.txt lesen
	string term = directory + "/term.txt";
	std::ifstream terms( term.c_str() ) ;
	if ( ! terms ) {
		Rcerr << "Cannot open " << term << "." << endl ;
		//exit( 1 ) ;
	}
	idmap id_to_go( terms ) ;
	terms.close(  ) ;
	Rcerr << "Read " << id_to_go.size() << " terms." << endl ;
	
	// steffi: graph_path.txt lesen	
	string graph_path = directory + "/graph_path.txt";
	std::ifstream transition_graph( graph_path.c_str() ) ;
	if ( ! transition_graph ) {
		Rcerr << "Cannot open " << graph_path << "." << endl ;
		//exit( 1 ) ;
	}
	//steffi:
	string parent_go( root ) ;
	string parent_id = id_to_go.get_id_for_go( parent_go ) ;
	transitions trans( parent_id, transition_graph ) ;
	transition_graph.close(  ) ;
	Rcerr << "Found " << trans.size() << " nodes." << endl ;
	
	// steffi: term2term lesen
	string termtoterm = directory + "/term2term.txt";
	std::ifstream term2term( termtoterm.c_str() ) ;
	if ( ! term2term ) {
		Rcerr << "Cannot open " << termtoterm << "." << endl ;
		//exit( 1 ) ;
	}
	go_graph graph( trans, term2term, id_to_go ) ;
	term2term.close(  ) ;
	Rcerr << "Graph created." << endl ;

	/*****************
         * read gene information and annotate to graph
	 *******************/
    //steffi: 
	std::ifstream annf( nodes_per_gene.c_str() ) ;
	if ( ! annf ) {
		Rcerr << "Cannot open " << nodes_per_gene << "." << endl ;
		//exit( 1 ) ;
	}
	//steffi:
	string gene_scores = directory + "/infile-data";
	std::ifstream dataf( gene_scores.c_str() ) ;
	if ( ! dataf ) {
		Rcerr << "Cannot open " << gene_scores << "." << endl ;
		//exit( 1 ) ;
	}

	genes gns( graph, annf, dataf ) ;
	Rcerr << "Data and annotation file parsed." << endl ;

	Rcout << "Number of randomsets: " << number_of_randomsets << "." <<endl;
	
	// steffi:
	Rcout << "Computing randomsets..." << number_of_randomsets << "." <<endl;

	string outfile = directory + "/randset_out";
	ofstream out;
	out.open ( outfile.c_str() );


	// force full length of all written floats...
	out.precision( 100 ) ; 

	/*****************
         * save original data, create and save randdata to randomsetfile, 
	 *******************/
	out << gns.sumnties() << endl ;
	graph.print_header( out ) ;
	graph.print_sumranks( out ) ;

	for ( int i=1 ; i <= number_of_randomsets ; ++i ) {
		graph.clear_genes(  ) ;
		gns.create_random_set(  ) ;
		graph.print_sumranks( out ) ;
	}

}
