/*
 * =====================================================================================
 *
 *       Filename:  merge_utls.hpp
 *
 *    Description:  header file of merge utls
 *
 *        Version:  1.0
 *        Created:  22/12/2017 12:22:50
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#ifndef _MERGE_UTLS_H
#define _MERGE_UTLS_H

#include "gfakluge.hpp"
#include "parse_paf.hpp"

using namespace std;
using namespace gfak;


int proc_blk(aln_block *abk, int u_s, GFAKluge &g, int bk_thres, int *bk_count, int *jn_count);

#endif

