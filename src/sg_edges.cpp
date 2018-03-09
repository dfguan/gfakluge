/*
 * =====================================================================================
 *
 *       Filename:  gfaAfa2edges.cpp
 *
 *    Description:  combine edges in sg.gfa and sequences in pread4fasta 
 *
 *        Version:  1.0
 *        Created:  08/03/2018 12:17:21
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#include "sg_edges.hpp"
#include <fstream>

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

string get_seq_id(string l)
{
	return l.substr(1,l.size()-1);//remove the first and last newline character 
}

string rc_dna_seq(string &s, int b, int e)
{
	string t = "";
	string temp = s.substr(b, e - b);
	for (int i = e - b - 1; i >= 0; --i) 
		t += rc_tab[temp[i]];
	return t;	
}
int updated_gfa_sequences(GFAKluge &g, string f_n)
{
	ifstream fp(f_n);
	if (fp.is_open()) {
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
							z += line.substr(0, line.size()-1);
						else 
							break;	
					}
					it->second.sequence = z;
				} else {
					while (getline(fp, line)) 
						if (line[0] == '>')	
							break;	
				}	
			}
		}
		return NORMAL;
	} else {
		return IO_ERR;
	}	
}

int extract_gfa(GFAKluge &g)
{
		
	map<string, vector<edge_elem>> n_2_e = g.get_seq_to_edges();
	map<string, sequence_elem, custom_key> n_2_s = g.get_name_to_seq();
	map<string, path_elem> pp = g.get_name_to_path();
	int count ; 
	for (auto it : pp) {
		sequence_elem s;
		s.name = it.first;
		path_elem &p = it.second;
		int		ol ;
		//init first path element
		string pre_seq_id = p.segment_names[0];
		string cur_seq_id;
		//s.length = 0;
		//s.sequence = "";
		int cor_e, cor_s;	
		bool isRc;
		for (size_t i = 1 ; i < p.segment_names.size(); ++i) {
			for (auto a : n_2_e[pre_seq_id]) {
				if (a.sink_name == p.segment_names[i]) {
					int len = n_2_s[pre_seq_id].length;
					ol = stoi(a.tags["ol"].val);
					isRc = p.orientations[i-1];
					cor_s = isRc ? ol: 0;
					cor_e = isRc ? len : len - ol; 
				}
			}
			cout<<">E"<<setw(7)<<count<<endl;
			cout<<(isRc ? rc_dna_seq(n_2_s[pre_seq_id].sequence, cor_s, cor_e) : n_2_s[pre_seq_id].sequence.substr(cor_s, cor_e - cor_s))<<endl;  
			++count;
			pre_seq_id = p.segment_names[i];
		}
	}	
	return NORMAL;
}


int main(int argc, char *argv[])
{


	if (argc < 3) {
		cerr << "convert_falcon_gfa <GFA_FILE> <Preads.fasta>" << endl;
		return 1;
	}
	
	string	gfa_fl_name = argv[1];	
	string	fa_fl_name = argv[2];
	GFAKluge gg;
	gg.parse_gfa_file(gfa_fl_name);
	if (updated_gfa_sequences(gg, fa_fl_name)) {
		fprintf (stderr,"File open Error\n");
		return IO_ERR;
	};
	if (extract_gfa(gg)) {
		fprintf(stderr,"File format Error\n");
		return FL_FORMAT_ERR;	
	}
	return NORMAL;	
}




