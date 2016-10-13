
// draw unique integers [1,total length of genes] and select genes dependent on gene length

#include <set>
//#include <stdlib.h>     // srand, rand 
#include <vector>
#include <map>
#include <iostream>
#include "structures.h"

#include <Rcpp.h>

std::set<int> rannum_genelen(int n_candidate, const std::map<std::string,int> &genename_to_index, std::vector<gen_pos_str> genes_pos){
	
	// go through all genes and add length and cumulative length
	//Rcpp::Rcout << std::endl << "ran_genelen option:" << std::endl << std::endl;
	long total_length = 0;
	for (int i=0; i < genes_pos.size(); i++){
		total_length += genes_pos[i].end - genes_pos[i].start;
		genes_pos[i].cumu_len = total_length;
		//if (i < 50) {
			//Rcpp::Rcout << genes_pos[i].name << ", len: " << genes_pos[i].end - genes_pos[i].start << ", cumulen: " << genes_pos[i].cumu_len << std::endl;
		//}	
	}		
	// get random genes and their index, without duplicates
	std::set<int> random_numbers;
	// draw as many unique numbers of candidate genes
	while (random_numbers.size() < n_candidate) {	
		// get random number [1,total length of genes]
		//long ran = rand() % total_length + 1;	
		long ran = R::runif(0,1) * (total_length) + 1; 
		// select gene
		int k = 0;
		while (ran > genes_pos[k].cumu_len){
			k++;
		}	
		random_numbers.insert(genename_to_index.find(genes_pos[k].name)->second);
	}	
	//Rcpp::Rcout << "Candidate genes: " << n_candidate << ", Random genes: " << random_numbers.size() << std::endl;
	return(random_numbers);
}
