#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctype.h>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

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

// uses function from samtools to sort read names, idfk why they don't do a simple 
// string comparison
int main () {
    vector <string> readNames;
    string line;
    while (getline(std::cin, line)) {
        readNames.push_back(line);
    }
    sort(readNames.begin(), readNames.end(), sortStrNum);
    for(size_t i = 0; i < readNames.size(); i ++) {
        // cerr << sortStrNum(readNames[0], readNames[i]) << endl;
        cout << readNames[i] << endl;
    }
    return 0;
}