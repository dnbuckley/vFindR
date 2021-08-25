#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sstream>
#include <vector>
#include <thread>
#include <algorithm>

using std::string;
using std::stringstream;
using std::cerr;
using std::cout;
using std::endl;

// only outputs reads in text file - argv[1]
// sam file and readnames MUST be name sorted
int main (int argc, char* argv[]) {
	if (argc != 2) {
		cerr << "provide read names." << endl;
		return 1;
	}
	std::ifstream readNameFile(argv[1]);
	string samReadName, refReadName, samLine;
	// get first readname
	getline(readNameFile, refReadName);
	cerr << "refRead = " << refReadName << endl;
	// go through sam file
	while (getline(std::cin, samLine)) {
		// check for header
		if (samLine.front() == '@') {
			cout << samLine << endl;
			continue;
		}
		std::stringstream ss(samLine);
		ss >> samReadName;
		// cerr << samReadName << "\t" << refReadName << endl;
		// if matched output name
		if (samReadName == refReadName) {
			cout << samLine << endl;
		} else {
			// if no match check if sam > ref
			if (samReadName > refReadName) {
				getline(readNameFile, refReadName);
				// check if next read is a match
				if (samReadName == refReadName) {
					cout << samLine << endl;
				}
			}
		}
		// get the next sam line
	}
	if (samReadName < refReadName) {
		cerr << "this should not be possible" << endl;
		return 1;
	}
	return 0;
}
