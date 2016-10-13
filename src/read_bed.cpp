
// Read bedfile and store in vector of regions
	
#include "structures.h"
#include <fstream>
#include <sstream>
#include <vector>


std::vector< bed_str > read_bed(std::string bed_file){
	std::vector< bed_str > desert_bed;
	std::ifstream desert (bed_file.c_str());
	std::string line;
	while(std::getline( desert, line )){
		bed_str bed;
		std::istringstream is( line.c_str()) ;
		is >> bed.chrom >> bed.start >> bed.end;	
		bed.len = bed.end - bed.start; // lengths of deserts			
		desert_bed.push_back(bed);  
	}
	desert.close();
	return(desert_bed);
}
