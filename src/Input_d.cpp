#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <algorithm>
#include "Input_d.h"
using namespace std;

Input_d::Input_d(const char* FilePath){
    ifstream InputFile;
    stringstream InputStream;
    try{
        InputFile.open(FilePath);
        InputStream << InputFile.rdbuf();
        InputFile.close();
    }
    catch (ifstream::failure& e) { 
        cout << "ERROR::FILE_NOT_SUCCESSFULLY_READ:" << e.what() << endl;
    }
    string ASentence;
    while (getline(InputStream,ASentence)){
        int pos1 = ASentence.find(" ");
        string Premeter = ASentence.substr(0, pos1);
        int pos2 = ASentence.find(" ", pos1 + 1);
        string Operation = ASentence.substr(pos1 + 1, pos2 - pos1);
        if (Premeter == "cutoff") cutoff = stod(Operation);
        if (Premeter == "dr") dr = stod(Operation);
        if (Premeter == "mesh") mesh = stoi(Operation);
        if (Premeter == "f:"){
            f = new double[mesh];
            for (int i = 0 ; i < mesh ; i ++){
                InputStream >> f[i];
                while (((InputStream.peek() < '0' || InputStream.peek() > '9') && InputStream.peek() != '-') && i < mesh-1) InputStream.ignore();
            }
        }
    }
}

Input_d::~Input_d(){
    delete[] f;
    f = NULL;
}