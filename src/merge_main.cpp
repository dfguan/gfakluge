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
	//fprintf(stderr, "%s\n", o.gfa_fn.c_str());
    gg.parse_gfa_file(o.gfa_fn);
    gg.id_analyze();
	paf_parser pp;
    if (pp.open_file(o.map_fn)) {
		fprintf(stderr, "IO Error, Please Check your Paf file\n");
		return IO_ERR;
	}
	int bk_count, jn_count;
	bk_count = jn_count = 0;	
    int unit_size;
    while ((unit_size = pp.read_next_block())) { 
        aln_block * a = pp.get_blk();
		//for (int i = 0; i < unit_size; ++i) cerr<< a->alns[i].seq_id<<"\t"<<a->alns[i].ref_s<<"\t"<<a->alns[i].ref_e<<endl;
		//fprintf(stderr,"%d\n", unit_size);
		proc_blk(1, a, unit_size, gg, o.bk_thres,o.map_len_thres, &bk_count, &jn_count);     
    }

    pp.close_file();
    //print statistics
	fprintf(stderr, "[merge stats]: %d break point, %d joint point\n", bk_count, jn_count);
	
	//write new gg out to stdout
	cout<<gg.block_order_string_2(); 
    
	return NORMAL;
}



