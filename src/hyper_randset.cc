/*
 * FUNC - Functional Analysis of Gene Expression Data
 * Copyright (C) 2002  Bjoern Muetzel, Kay Pruefer
 * 
 * This program is modifiable/redistributable under the terms of the
 * GNU General Public License.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

//steffi:
//using namespace std ;

#include <time.h>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <memory>

// steffi:
#include "go.h"
#include "go_graph_hyper.h"
#include "idmap.h"
#include "transitions.h"

#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

#define MAX_LINE_LENGTH 20000

// steffi:
// int main( int argc, char *argv[] ) 

//[[Rcpp::export]]
void hyper_randset(std::string detected, std::string changed ,int number_of_randomsets, std::string outfile, std::string term, std::string graph_path, std::string termtoterm, std::string root) 

{
	/*steffi:
	if (  argc != 9 ) {
		cerr << "Usage: " << argv[0] << " detected "
				"changed number_of_sets outfile term.txt" << endl
				<< "       graph_path.txt term2term.txt GO_ID" << endl ;
		//exit( 1 ) ;
	}*/
	
	// Build GO-Graph using different files from go_date_termdb-tables.tar.gz
	// steffi: term.txt lesen
	std::ifstream terms( term.c_str() ) ;
	if ( ! terms ) {
		Rcerr << "Cannot open " << term << "." << endl ;
		//exit( 1 ) ;
	}
	idmap id_to_go( terms ) ;
	terms.close(  ) ;
	Rcerr << "Read " << id_to_go.size() << " terms." << endl ;

	// steffi: graph_path.txt lesen
	std::ifstream transition_graph( graph_path.c_str() ) ;
	if ( ! transition_graph ) {
		Rcerr << "Cannot open " << graph_path << "." << endl ;
		//exit( 1 ) ;
	}
	// steffi:
	string parent_go( root ) ;
	string parent_id = id_to_go.get_id_for_go( parent_go ) ;
	transitions trans( parent_id, transition_graph ) ;
	transition_graph.close(  ) ;
	Rcerr << "Found " << trans.size() << " nodes." << endl ;

	// steffi: term2term lesen
	std::ifstream term2term( termtoterm.c_str() ) ;
	if ( ! term2term ) {
		Rcerr << "Cannot open " << termtoterm << "." << endl ;
		//exit( 1 ) ;
	}
	go_graph_hyper graph( trans, term2term, id_to_go ) ;
	term2term.close(  ) ;
	Rcerr << "Graph created." << endl ;
	
	// gos-object will be used to get one int* per GO and to print 
	// the results.
	go gos ;
	
	/* steffi: output immer file - delete am ende immer möglich
	ostream *out ;
	if ( string( argv[4] ) == "-" ) {
		out = &cout ;
	} else {
		out = new ofstream( argv[4] ) ;
	}
	if ( !*out ) {
		cerr << "Cannot open output file " << argv[4] << endl ;
		//exit( 1 ) ;
	}*/
	std::ostream *out ;
	out = new std::ofstream( outfile.c_str() ) ;

	// steffi:	
	std::ifstream in( detected.c_str() ) ;

	// steffi:
	Rcerr << "Reading detectedfile... " 
			<< endl ;
	
	// gens == genes. This vector is a simple representation of the go tree.
	// every gene is 1 vector of int*, where the int represents one go-node.
	vector<vector<int*> > gens;


	map<string,int> genename_to_index ;
	int index = 0 ;

//	char line[MAX_LINE_LENGTH] ;
	string line ;
	while ( in ) {
		getline( in, line ) ;
//		in.getline(line, MAX_LINE_LENGTH, '\n' ) ;
		std::istringstream is( line.c_str() ) ; 
		string gen_name ;
		is >> gen_name ; 
#ifdef DEBUG
		// Steffi:
		Rcout << "gen_name: " << gen_name << endl ;
#endif

		vector<int*> gen_vec;
		string go_name ;
		set<string> parents ;
		while ( is >> go_name ) {
			// Get the names of all nodes that are parents of go_name
			graph.get_parents( go_name, &parents ) ;
		}
		for ( set<string>::const_iterator it = parents.begin() ; 
				it != parents.end() ; ++it ) 
		{
			// gos.add returns a unique int* for every string.
			gen_vec.push_back( gos.add( *it ) ) ; 
#ifdef DEBUG
			// Steffi:
			Rcout << "go: " << *it << endl ;
#endif
		}
		// if the gene is annotated, add it to the genes-vector
		if ( gen_vec.size() > 0 ) {
			gens.push_back( gen_vec ) ;
			genename_to_index[gen_name] = index ;
			index++ ;	
		}		

	}
	// steffi:
	Rcerr << "Found " << gens.size() << " usable entrys in " << detected
		<< " with " << gos.size() << " GOs" << endl ;
	
	*out << "Genes:\t" << gens.size() << endl ;
	*out << "GOs:\t" << gos.size() << endl ;
	
	gos.print_names( *out ) ;

	// add changedfile...
	// steffi:
	std::ifstream changed_in( changed.c_str() ) ;
	if ( ! changed_in ) {
		Rcerr << "Cannot open " << changed << endl ;
		//exit( 1 ) ;
	}

	string line_s ;

	int size_of_random_sets = 0 ;

	while( changed_in ) {
		getline( changed_in, line_s ) ;
		std::istringstream is( line_s.c_str() ) ;
		string gen_name ;
		is >> gen_name ;
		if ( genename_to_index.find( gen_name ) != genename_to_index.end() ) {
			for ( vector<int*>::iterator it = gens[genename_to_index[gen_name]].begin() ;
					it != gens[genename_to_index[gen_name]].end() ; ++it ) {
				(*(*it))++ ;
			}
			size_of_random_sets++ ;
		}
	}
	gos.print_sum( *out ) ;
	gos.clear() ;

	/*steffi: number_of_randomsets schon als int-argument 
	int number_of_randomsets ;	
	{
		istringstream is( argv[3] ) ;
		is >> number_of_randomsets ;
	}*/

/*	int size_of_random_sets ;
	{
		istringstream is( argv[2] ) ;
		is >> size_of_random_sets;
	}
*/
	//steffi:
	Rcerr << "Creating " << number_of_randomsets << " randomsets with "
			"size " << size_of_random_sets << endl ;

/*	*out << "Randomsets:\t" << number_of_randomsets << endl ;
	*out << "Genes per randomset:\t" << size_of_random_sets << endl ;
*/




	// forall randomsets
	for ( int i = 1 ; i <= number_of_randomsets ; ++i ) {
		// create randomset
		set<int> random_numbers ;
		double max = gens.size() ;
		for ( int randi = 1 ; randi <= size_of_random_sets ; ++randi ) {
		
			int rand_num = static_cast<int>(max*R::runif(0,1)) ;
			while ( random_numbers.find( rand_num ) !=
												random_numbers.end() ) {
				rand_num = static_cast<int>(max*R::runif(0,1)) ;
			}
			random_numbers.insert( random_numbers.begin(), rand_num ) ;
		}
		// reset all go-nodes
		gos.clear(  ) ;

		// add 1 to every GO that a randomly choosen gene is part of
		for ( set<int>::const_iterator it = random_numbers.begin() ;
				it != random_numbers.end() ; ++it ) {
			for ( vector<int*>::iterator it2 = gens[*it].begin() ;
					it2 != gens[*it].end() ; ++it2 ) {
				(*(*it2))++ ;
			}
		}
		random_numbers.clear() ;
		
		// print a line with the values for every go
		gos.print_sum( *out ) ;

	}
	Rcerr << "\rFinished" << endl ;
	// steffi:
	delete out;
}
