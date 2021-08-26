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

static int strnum_cmp(const char *_a, const char *_b)
{
    const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
    const unsigned char *pa = a, *pb = b;
    while (*pa && *pb) {
        if (isdigit(*pa) && isdigit(*pb)) {
            while (*pa == '0') ++pa;
            while (*pb == '0') ++pb;
            while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
            if (isdigit(*pa) && isdigit(*pb)) {
                int i = 0;
                while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
                return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
            } else if (isdigit(*pa)) return 1;
            else if (isdigit(*pb)) return -1;
            else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
        } else {
            if (*pa != *pb) return (int)*pa - (int)*pb;
            ++pa; ++pb;
        }
    }
    return *pa? 1 : *pb? -1 : 0;
}

bool sortStrNum(const string a, const string b) {
    int t = strnum_cmp(a.c_str(), b.c_str());
    return t < 0;
}

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
			if (!sortStrNum(samReadName, refReadName)) {
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
