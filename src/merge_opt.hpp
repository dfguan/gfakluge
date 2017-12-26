/*
 * =====================================================================================
 *
 *       Filename:  merge_opt.h
 *
 *    Description:  merge options 
 *
 *        Version:  1.0
 *        Created:  22/12/2017 10:52:48
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#ifndef _MERGE_OPT_H
#define _MERGE_OPT_H


#include <string>
using namespace std;

typedef struct _opt {
    bool         useful;
    int          bk_thres;   
    string       gfa_fn;
    string       map_fn;
}opts;


int help();
bool get_params(int argc ,char *argv[], opts *o);

#endif
