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
#include "ds.hpp"
//PS
#define BREAK 1
#define CONTINUE 2
#define END_Q 3


using namespace std;
using ds::ary;

typedef struct _aln_unit{
    string  pat_id;
    int     pat_len;
    bool	isForward;
	bool	isCovered;
	int     qry_s;//data type might required to be extended 
    int     qry_e;
    int     pat_s;
    int     pat_e;
    _aln_unit& operator= (const _aln_unit& e){
        pat_id = e.pat_id;
        pat_len = e.pat_len;
		isForward = e.isForward;
		isCovered = e.isCovered;
		qry_s = e.qry_s;
        qry_e = e.qry_e;
        pat_s = e.pat_s;
        pat_e = e.pat_e;
        return *this;
    }
	//see if one query is covered in another
	bool operator<(const _aln_unit& r) const {
		if (qry_s < r.qry_s) return true;
		else if (qry_s == r.qry_s) return qry_e > r.qry_e;
		else return false;
	}
	int set_covered() { isCovered = true;return NORMAL;}
}aln_unit;

typedef struct _aln_block{
    string      qry_id;
    int         qry_len;
    int         aln_unit_limit;
    aln_unit    *alns;
    bool        isInit;
    int         aln_unit_index;
	bool		isSorted;
	bool		isSetCovered;
    _aln_block() {
        //allocate space for aln unit
        isInit = false;
        aln_unit_limit = 10; 
        aln_unit_index = 0;
        alns = new aln_unit[aln_unit_limit];
        if (alns != NULL) isInit = true; 
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
    // sort reads coordinates
	int sort_1() {
		if (aln_unit_index > 1) 
			sort(alns, alns + aln_unit_index - 1);
		isSorted = true;
		return NORMAL;
	}
	//set isCovered in aln_unit if qry_e > next.qry_e  
	int set_is_covered() {
		if (aln_unit_index > 1) {
			int pre_qry_e = alns[0].qry_e;	
			for (int i = 1;  i < aln_unit_index; ++i) {
				if (pre_qry_e >= alns[i].qry_e) 
					alns[i].set_covered();
				else
					pre_qry_e = alns[i].qry_e;
			}
		}  
		isSetCovered = true;
		return NORMAL;
	}
	//merge alignments requires to be sorted and set covered first
	int merge_alns(int thres, int thres_mapped_len) 
	{
		if (!isSorted) sort_1();
		if (!isSetCovered) set_is_covered();
		ary<int> c,g;
		int c_k,g_k;
		c_k = g_k = 0;
		for (int i = 0; i < aln_unit_index; ++i) if (!alns[i].isCovered) c[c_k++] = i; 
		
		g[g_k++] = 0;
		for (int i = 0; i < c_k - 1; ++i) {
			int p = c[i];
			int q = c[i+1];	
			if (!(alns[p].pat_id == alns[q].pat_id && alns[p].isForward == alns[q].isForward && (alns[p].isForward ? alns[q].pat_e > alns[p].pat_e : alns[p].pat_s > alns[q].pat_s) && (alns[q].qry_s - alns[p].qry_e) < thres && (alns[p].isForward ? (alns[q].pat_s - alns[p].pat_e) : (alns[p].pat_s - alns[q].pat_e)) < thres)) g[g_k++] = i + 1;
			//else {
				//cerr<<p<<"\t"<<q<<endl;
				//cerr<< alns[p].seq_id<<"\t"<<alns[p].ref_s<<"\t"<<alns[p].ref_e<<"\t"<<alns[p].seq_s<<"\t"<<alns[p].seq_e<<"\t"<<alns[p].isCovered<<"\t"<<alns[p].seq_e - alns[p].seq_s<<endl;
				//cerr<< alns[q].seq_id<<"\t"<<alns[q].ref_s<<"\t"<<alns[q].ref_e<<"\t"<<alns[q].seq_s<<"\t"<<alns[q].seq_e<<"\t"<<alns[q].isCovered<<"\t"<<alns[q].seq_e - alns[q].seq_s<<endl;
			//}
		}
		g[g_k] = c_k;	
		
		aln_unit t;
		int z = 0;
		for (int i = 0; i < g_k; ++i) {
			//fprintf(stderr, "%d\t%d\n", g[i], g[i+1]-1);
			int s = c[g[i]];
			int e = c[g[i+1] - 1];
			t.pat_s = alns[s].isForward ? alns[s].pat_s : alns[e].pat_s;
			t.pat_e = alns[s].isForward ? alns[e].pat_e : alns[s].pat_e;
			if (t.pat_e - t.pat_s < thres_mapped_len) 
				continue;
			t.qry_e = alns[e].qry_e;
			t.qry_s = alns[s].qry_s;
			t.pat_id = alns[s].pat_id;
			t.pat_len = alns[s].pat_len;	
			t.isCovered = false;
			t.isForward = alns[s].isForward;
			alns[z++] = t;
		}	
		//aln_unit_index = g_k;
		return z;	
	}	
	//release aln_block
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
	//bool			patIsQuery; //this occurs when one is trying to align reference to reads 
public:
    paf_parser() {}
    paf_parser(string fl_name) 
    {
        fl.open(fl_name);
        //read the first line;
        getline(fl, cur_ln);
        //fl.seekg(0, ios::beg);
        tokens = split(cur_ln, '\t');
        pre_id = tokens[0];
		//refIsQuery = r;

    }    
    //int set_refIQ(bool r) {
		//refIsQuery = r;
		//return NORMAL;	
	//}
	int open_file(string fl_name) 
    {
       if (fl.is_open()) fl.close();
       fl.open(fl_name);
       if (fl.is_open()) {
           //read the first line
            getline(fl, cur_ln);
            //fl.seekg(0, ios::beg);
            tokens = split(cur_ln, '\t');
			//pre_id = tokens[5];
			pre_id = tokens[0];
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
            if (tokens[0] == pre_id) {
                aln_unit u;
				//ref_id, ref_len, ref_start, ref_end,seq_id, seq_len, seq_s, seq_e 
				aln_blk.qry_id = tokens[0]; //this is query's id actually
				aln_blk.qry_len = stoi(tokens[1]); //query_length
				u.qry_s = stoi(tokens[2]);
				u.qry_e = stoi(tokens[3]);
				u.isForward = tokens[4] == "+"? true : false; 
				u.isCovered = false;
				u.pat_id = tokens[5];
				u.pat_len = stoi(tokens[6]);
				u.pat_s = stoi(tokens[7]);
				u.pat_e = stoi(tokens[8]);
                aln_blk.push(&u);
                return CONTINUE;
            } else {
                pre_id = tokens[0];
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
        } else 
			aln_blk.clear();
        return aln_blk.aln_unit_index;
    }
    aln_block * get_blk() 
    {
        return &aln_blk;    
    }
};
#endif

