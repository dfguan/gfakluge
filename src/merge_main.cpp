/*
 * =====================================================================================
 *
 *       Filename:  gfa_merge.cpp
 *
 *    Description:  give mapping file of two gfa files, merge them 
 *
 *        Version:  1.0
 *        Created:  21/12/2017 14:00:14
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include <getopt.h>

//parse_paf.hpp and gfakluge.hpp included in merge utls.hpp
#include "merge_utls.hpp"

#include "merge_opt.hpp"

using namespace std;
using namespace gfak;

int main(int argc, char *argv[]) 
{
    opts o;
    //get the parameters we want
    if (!get_params(argc, argv, &o)) 
        exit(0);
    GFAKluge gg;
    gg.parse_gfa_file(o.gfa_fn);
    
    paf_parser pp;
    pp.open_file(o.map_fn);
    
    int unit_size;
    while ((unit_size = pp.read_next_block())) { 
        aln_block * a = pp.get_blk();
        proc_blk(a, unit_size, gg, o.bk_thres);     
    }

    pp.close_file();
    //write new gg out to stdout
     
    return NORMAL;
}



