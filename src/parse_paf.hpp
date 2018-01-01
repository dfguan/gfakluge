/*
 * =====================================================================================
 *
 *       Filename:  parse_paf.hpp
 *
 *    Description:  header of parse_paf.cpp
 *
 *        Version:  1.0
 *        Created:  21/12/2017 14:05:54
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#ifndef _PARSE_PAF_H
#define _PARSE_PAF_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdio.h>

#include "status_code.hpp"

//PS
#define BREAK 1
#define CONTINUE 2
#define END_Q 3


using namespace std;

typedef struct _aln_unit{
    string  seq_id;
    int     seq_len;
    bool    isConverted;
    int     seq_s;//data type might required to be extended 
    int     seq_e;
    int     ref_s;
    int     ref_e;
    _aln_unit& operator= (const _aln_unit& e){
        seq_id = e.seq_id;
        seq_len = e.seq_len;
        seq_s = e.seq_s;
        seq_e = e.seq_e;
        ref_s = e.ref_s;
        ref_e = e.ref_e;
        return *this;
    }
}aln_unit;

typedef struct _aln_block{
    string      ref_id;
    int         ref_len;
    int         aln_unit_limit;
    aln_unit    *alns;
    bool        isInt;
    int         aln_unit_index;
    _aln_block() {
        //allocate space for aln unit
        isInt = false;
        aln_unit_limit = 10; 
        aln_unit_index = 0;
        alns = new aln_unit[aln_unit_limit];
        if (alns != NULL) isInt = true; 
    }
    //double size of aln_unit
    aln_unit* reallocate() {
        int s = aln_unit_limit << 1;
        aln_unit *new_a = new aln_unit[s];
        if (new_a) {
            for (int i = 0; i < aln_unit_limit; ++i ) new_a[i] = alns[i];
            delete []alns;
            aln_unit_limit = s;
        }
        return new_a;
    }
    //push back aln_unit
    int push(aln_unit *e) {
        if (aln_unit_index >= aln_unit_limit) {
            alns = reallocate();
           if (alns) {
                alns[aln_unit_index++] = *e;                    
           } else 
               return RES_ERR;
        } else 
            alns[aln_unit_index++] = *e;                    
        return NORMAL; 
    }

    int pop() {
        if (aln_unit_index) {--aln_unit_limit; return NORMAL;}
        else 
            return END_Q;
    }
    int clear() {
        aln_unit_index = 0;
        return NORMAL;
    }
    ~_aln_block() {
        if (alns) delete []alns;
    }

}aln_block;

class paf_parser {

private:
    aln_block       aln_blk;//
    string          cur_ln; // current line
    string          pre_id;     // previous reference id
    ifstream        fl;
    vector<string>  tokens;
public:
    paf_parser() {}
    paf_parser(string fl_name) 
    {
        fl.open(fl_name);
        //read the first line;
        getline(fl, cur_ln);
        fl.seekg(0, ios::beg);
        tokens = split(cur_ln, '\t');
        pre_id = tokens[5];
    }    
    int open_file(string fl_name) 
    {
       if (fl.is_open()) fl.close();
       fl.open(fl_name);
       if (fl.is_open()) {
           //read the first line
            getline(fl, cur_ln);
            fl.seekg(0, ios::beg);
            tokens = split(cur_ln, '\t');
			pre_id = tokens[5];
            return NORMAL;
       } else 
           return IO_ERR;
    }
    int close_file() 
    {
        if (fl.is_open()) fl.close();
        return NORMAL;
    }
    vector<string> split(string s, char delim)
    {
        vector<string> ret;
        stringstream sstream(s);
        string temp;
        while(getline(sstream, temp, delim)){
            ret.push_back(temp);
        }        
        return ret;
    }
    
    int process(string& l)
    {
		//fprintf(stderr, "%s\n",l.c_str());
		if (l != "") {
            tokens = split(l, '\t');
            if (tokens[5] == pre_id) {
                aln_unit u;
                u.seq_id = tokens[0];
                u.seq_len = stoi(tokens[1]);
                u.seq_s = stoi(tokens[2]);
                u.seq_e = stoi(tokens[3]);
                aln_blk.ref_id = tokens[5];
                aln_blk.ref_len = stoi(tokens[6]);
                u.ref_s = stoi(tokens[7]);
                u.ref_e = stoi(tokens[8]);
                //coverted coordinate if reverse complementary
                if (tokens[4] == "-") {
                    u.isConverted = true;
                    int t = u.seq_s;
                    u.seq_s = u.seq_len - u.seq_e;
                    u.seq_e = u.seq_len - t; 
                } 
                aln_blk.push(&u);
                return CONTINUE;
            } else {
                pre_id = tokens[5];
                return BREAK; 
            }
        } else 
            return BREAK;
    }
    int read_next_block() 
    {
        aln_blk.clear();
        if (process(cur_ln) != BREAK) {
            while (getline(fl, cur_ln)) {
                if (process(cur_ln) == BREAK)
                    break;          
            }
        }
        return aln_blk.aln_unit_index;
    }
    aln_block * get_blk() 
    {
        return &aln_blk;    
    }
};
#endif

