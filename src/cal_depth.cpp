/*
 * =====================================================================================
 *
 *       Filename:  cal_depth.cpp
 *
 *    Description:  calculate depth of reference based on paf file
 *
 *        Version:  1.0
 *        Created:  13/03/2018 10:33:22
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include "parse_paf.hpp"
#include "gfakluge.hpp"

using namespace gfak;

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage : cal_dp <GFA_FILE> <PAF_FILE>\n");
		return LACK_PARAMS; 
	}	
	string gfa_fl_name = argv[1];
	string paf_fl_name = argv[2];
	GFAKluge  gg;
	gg.parse_gfa_file(gfa_fl_name);
	paf_parser pp;
    if (pp.open_file(gfa_fl_name)) {
		fprintf(stderr, "IO Error, Please Check your Paf file\n");
		return IO_ERR;
	}
	int unit_size;
	map<string, uint64_t> seq_to_depth;
	map<string, sequence_elem, custom_key> name_to_seq = gg.get_name_to_seq();
	while ((unit_size = pp.read_next_block())) {
		aln_block *a = pp.get_blk();	
		aln_unit *b;
		for (int i = 0; i < unit_size; ++i) {
			b = a->alns + i;
			seq_to_depth[b->pat_id] += b->pat_e - b->pat_s;	
		} 	
	}
	for (auto it : name_to_seq) {
		opt_elem o = {"DP","F",to_string(double(seq_to_depth[it.first])/it.second.length)};
		gg.add_tag(it.first, &o);
	}
	cout<<gg.block_order_string();
	
}
