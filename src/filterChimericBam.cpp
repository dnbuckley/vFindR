#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sstream>
#include <vector>

bool checkChrs(std::vector<std::string> chrs) {
	bool human, virus = 0;
	for (size_t i = 0; i < chrs.size(); i++) {
		if (chrs[i].substr(0, 3) == "chr") human = 1;
		if (chrs[i].substr(0, 3) == "NC_") virus = 1;
	}
	return (human && virus);
}

// filters for reads mapped to human and viral chromosomes
int main() {
	std::string read, chr, flag, curReadName, prevReadName;
	std::vector<std::string> reads, chrs;
	while(getline(std::cin, read)) {
		// check for header
		if (read.front() == '@') {
			std::cout << read << std::endl;
			prevReadName = read;
			continue;
		}
		std::stringstream ss(read);
		ss >> curReadName;
		ss >> flag;
		ss >> chr;
		if (curReadName != prevReadName) {
			if (checkChrs(chrs)) {
				for (size_t i = 0; i < reads.size(); i++) {
					std::cout << reads[i] << "\n";
				}
			}
			reads.clear();
			chrs.clear();
		}
		reads.push_back(read);
		chrs.push_back(chr);
		prevReadName = curReadName;
	}
	return 0;
}