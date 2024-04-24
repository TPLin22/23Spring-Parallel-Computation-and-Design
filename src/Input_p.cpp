#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <algorithm>
#include "Input_p.h"
using namespace std;

Input_p::Input_p(const char* FilePath){
    ifstream InputFile;
    stringstream InputStream;
    try{
        InputFile.open(FilePath);
        InputStream << InputFile.rdbuf();
        //LineStream << InputFile.rdbuf();
        InputFile.close();
    }
    catch (ifstream::failure& e) { 
        cout << "ERROR::FILE_NOT_SUCCESSFULLY_READ:" << e.what() << endl;
    }
    string content = InputStream.str();
    stringstream LineStream(content);
    n = 0;
    string ASentence;
    while (getline(LineStream,ASentence)) n ++;
    pos = new double*[n+1];
    for (int i = 0 ; i <= n ; i ++)
        pos[i] = new double[4];
    for (int i = 0 ; i < n ; i ++){
        getline(InputStream,ASentence);
        istringstream iss(ASentence);
        char leftk,comma1,comma2,rightk;
        iss>>leftk>>pos[i][1]>>comma1>>pos[i][2]>>comma2>>pos[i][3]>>rightk;
    }
}

Input_p::~Input_p(){
    for (int i = 0 ; i <= n ; i ++){
        delete[] pos[i];
        pos[i] = NULL;
    }
    delete[] pos;
    pos = NULL;
}