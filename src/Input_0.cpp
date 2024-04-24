#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include "Input_0.h"
using namespace std;

Input_0::Input_0(const char* FilePath){
    ifstream InputFile;
    stringstream InputStream;
    //use ifstream to open file,put it into flow InputStream
    try {
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
        if (Premeter == "diago_lib") diago_lib = Operation;
        if (Premeter == "points_path") points_path = Operation;
        if (Premeter == "venergy_path") v_path = Operation;
        if (Premeter == "v_path") v_path = Operation;
        if (Premeter == "distribution_path") distribution_path = Operation;
        if (Premeter == "lx") lx = stod(Operation);
        if (Premeter == "ly") ly = stod(Operation);
        if (Premeter == "lz") lz = stod(Operation);
    }
    if (diago_lib[diago_lib.length()] == ' ') diago_lib.pop_back();
    if (v_path[v_path.length()] == ' ') v_path.pop_back();
    if (points_path[points_path.length()] == ' ') points_path.pop_back();
    if (distribution_path[distribution_path.length()] == ' ') distribution_path.pop_back();
}

Input_0::~Input_0(){}