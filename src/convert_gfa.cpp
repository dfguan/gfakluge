/*
 * =====================================================================================
 *
 *       Filename:  convert_gfa.cpp
 *
 *    Description:  convert GFA1 S->F, P->S(contigs) E with nothing 
 *
 *        Version:  1.0
 *        Created:  30/12/2017 09:41:42
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#include "gfakluge.hpp"
#include "status_code.hpp"
#include <string.h>
#include <getopt.h>

using namespace std;
using namespace gfak;


#define LIMITED_SIZE 100

char rc_tab[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 'T', 0, 'G', 0, 0, 0, 'C', 0, 0, 0, 0, 0, 0, 'N', 0, 
	0, 0, 0, 0, 'A', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 'T', 0, 'G', 0, 0, 0, 'C', 0, 0, 0, 0, 0, 0, 'N', 0, 
	0, 0, 0, 0, 'A', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};
int help() 
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: convert_gfa [options] <GFA_FILE>\n\n");
	fprintf(stderr, "options: -r FILE    required if a gfa file doesn't contain sequences\n");
	fprintf(stderr, "         -s         output non-overlapping string graph\n");
	fprintf(stderr, "         -h         print help information\n");
	return NORMAL;
}

string rc_dna_seq(string &s, int b, int e)
{
	string t = "";
	string temp = s.substr(b, e - b);
	for (int i = e - b - 1; i >= 0; --i) 
		t += rc_tab[temp[i]];
	return t;	
}

int check_seq(GFAKluge &g) 
{
	map<string, sequence_elem, custom_key> name_to_seq = g.get_name_to_seq();
	for (auto it : name_to_seq) {
		if (it.second.length != it.second.sequence.size())
			return LACK_KEY_INFO;
	}
	return NORMAL;
}
int output_seqs(vector<sequence_elem> &ss)
{
	for (auto s : ss) cout << s.to_string_2() << endl; 
	ss.clear();
	return NORMAL;
}


bool is_p_ctg(string g_id)
{
	for (size_t i = 0 ; i < g_id.length(); ++i) if (g_id[i] == '-') return false;
	return true;	
}
string zeroPadNum(int n)
{
	ostringstream ss;
	ss<<setw(7)<<setfill('0')<<n;
	return ss.str();
}
int conv_p_2_s_s(GFAKluge& gg)
{
	//convert to split sequence
	map<string, path_elem> pp = gg.get_name_to_path();		
	map<string, vector<edge_elem>> n_2_e = gg.get_seq_to_edges();
	map<string, sequence_elem, custom_key> n_2_s = gg.get_name_to_seq();
	int count = 0;
	for (auto it : pp) {
		if (is_p_ctg(it.first)) {
			vector<sequence_elem> ss;
			sequence_elem s;
			
			path_elem &p = it.second;
			string pre_seq_id = p.segment_names[0];
			string cur_seq_id;
			int cor_s, cor_e;
			bool isRc, isExist;

			s.name = "E"+zeroPadNum(count++); 
			s.length = n_2_s[pre_seq_id].length;
			s.sequence = n_2_s[pre_seq_id].sequence;
			//opt_elem o = {"LN","i",to_string(s.length)};
			//s.opt_fields.push_back(o);
			opt_elem o = {"sn","Z",pre_seq_id};
			s.opt_fields.push_back(o);
			ss.push_back(s); // the first one
		
			for (size_t i = 1 ; i < p.segment_names.size(); ++i) {
				isExist = false;
				for (auto a : n_2_e[pre_seq_id]) {
					if (a.sink_name == p.segment_names[i]) {
						cor_s = stoi(a.tags["ob"].val);
						cor_e = stoi(a.tags["oe"].val);
						isRc	=	false;
						isExist =	true;
					}
				}
				if (isExist) {
					if (cor_s > cor_e) {
						int t = cor_s;
						cor_s = cor_e;
						cor_e = t;
						isRc = true;	
					}	
				} else 
					return FL_FORMAT_ERR;
				sequence_elem t;
				t.name = "E"+zeroPadNum(count++); 
				t.sequence =  (isRc ? rc_dna_seq(n_2_s[p.segment_names[i]].sequence, cor_s, cor_e) : n_2_s[p.segment_names[i]].sequence.substr(cor_s, cor_e - cor_s));  
				t.length = cor_e - cor_s;
				pre_seq_id = p.segment_names[i];
				//o = {"LN","i",to_string(s.length)};
				//t.opt_fields.push_back(o); 
				o = {"sn","Z",pre_seq_id};
				t.opt_fields.push_back(o); 
				ss.push_back(t);
				if (ss.size() > LIMITED_SIZE) 
					output_seqs(ss);
			}
			if (ss.size()) 
				output_seqs(ss);
		}
	}	
	return NORMAL;

}
int conv_p_2_s(GFAKluge& gg)
{
	//get all path elements
	//if it is primary
	//segment id = g.id
	//contigs = first segment sequence
	//extract the first segment
	//pre = first
	//cur = current segment
	//extract edge(pre, cur);
	//if ob > oe cor_b = oe, cor_e = ob; rc = true
	//else cor_b = ob cor_e = oe, rc = false
	//temp_str = sub()
	//contigs += rc? rc():
	//length += 
	//output partof_contigs
	
	map<string, path_elem> pp = gg.get_name_to_path();		
	map<string, vector<edge_elem>> n_2_e = gg.get_seq_to_edges();
	map<string, sequence_elem, custom_key> n_2_s = gg.get_name_to_seq();
	vector<sequence_elem> ss;
	
	for (auto it : pp) {
		if (is_p_ctg(it.first)) {
			sequence_elem s;
			s.name = it.first;
			path_elem &p = it.second;
			//init first path element
			string pre_seq_id = p.segment_names[0];
			string cur_seq_id;
			int cor_s, cor_e;
			bool isRc, isExist;
			s.length = n_2_s[pre_seq_id].length;
			s.sequence = n_2_s[pre_seq_id].sequence;
			for (size_t i = 1 ; i < p.segment_names.size(); ++i) {
				isExist = false;
				for (auto a : n_2_e[pre_seq_id]) {
					if (a.sink_name == p.segment_names[i]) {
						cor_s = stoi(a.tags["ob"].val);
						cor_e = stoi(a.tags["oe"].val);
						isRc	=	false;
						isExist =	true;
					}
				}
				if (isExist) {
					if (cor_s > cor_e) {
						int t = cor_s;
						cor_s = cor_e;
						cor_e = t;
						isRc = true;	
					}	
				} else 
					return FL_FORMAT_ERR;
				s.sequence +=  (isRc ? rc_dna_seq(n_2_s[p.segment_names[i]].sequence, cor_s, cor_e) : n_2_s[p.segment_names[i]].sequence.substr(cor_s, cor_e - cor_s));  
				s.length += cor_e - cor_s;
				pre_seq_id = p.segment_names[i];
			}
			opt_elem o = {"LN","i",to_string(s.length)};
			s.opt_fields.push_back(o); 
			ss.push_back(s);
			if (ss.size() > LIMITED_SIZE) 
				output_seqs(ss);
		}
	}	
	if (ss.size()) 
		output_seqs(ss);
	return NORMAL;
}

string get_seq_id(string l)
{
	string t = l.substr(1, l.size()-1);
	int pos = t.find(' ');
		
	return t.substr(0, pos);//remove the first and last newline character 
}
int updated_gfa_sequences(GFAKluge &g, string f_n)
{
	ifstream fp(f_n);
	if (fp.is_open()) {
		//cout<<"file is open"<<endl;
		map<string, sequence_elem, custom_key> n_2_s = g.get_name_to_seq();
		map<string, sequence_elem, custom_key>::iterator it; 
		string line;
		getline(fp, line);	
		while (true) {
			if (fp.eof()) 
				break;
			if (line[0] == '>') {
				string seq_id =	get_seq_id(line);
				it = n_2_s.find(seq_id);
				if (it != n_2_s.end()) {
					string z = "";
					while (getline(fp, line)) {
						if (line[0] != '>')	
							z += line.substr(0, line.size());
						else 
							break;	
					}
					//cout<<z.size()<<endl;
					g.update_seq(seq_id, &z);
				} else {
					while (getline(fp, line)) 
						if (line[0] == '>')	
							break;	
				}	
			} 
		}
		fp.close();
		return NORMAL;
	} else {
		return IO_ERR;
	}	
}
int main(int argc, char *argv[])
{

	if (argc < 3) {
		help();
		return 1;
	}
	int c = 0;
	optind = 1;
	string ref_fl = 0; 
	bool s_c = false;
	while ((c = getopt(argc, argv, "r:sh")) < 0) {
		switch(c) {
			case 'r':
				ref_fl = optarg;
				break;
			case 's':
				s_c = true;
				break;
			case 'h':
				help();
				return 1;
			default:
				fprintf(stderr, "Undefined Params %c\n", c);
				help();
				return 1;
		}
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "[option_parse]: gfa file can't be ommited\n");
		help();
		return 1;
	}
	
	string	gfa_fl_name = argv[optind];


	GFAKluge gg;
	gg.parse_gfa_file(gfa_fl_name);
	
	if (ref_fl != "") {
		if (updated_gfa_sequences(gg, ref_fl)) {
			fprintf (stderr,"File open Error\n");
			return IO_ERR;
		};
	}	
	if (check_seq(gg)) {
		fprintf(stderr, "sequences error in gfa\n");	
		return LACK_KEY_INFO; 	
	}
	if (s_c) {
		if (conv_p_2_s_s(gg)) {
			fprintf(stderr,"File format Error\n");
			return FL_FORMAT_ERR;	
		}	
	}else {
		if (conv_p_2_s(gg)) {
			fprintf(stderr,"File format Error\n");
			return FL_FORMAT_ERR;	
		}
	}
	return NORMAL;	
}



