/*
 * =====================================================================================
 *
 *       Filename:  merge_utls.cpp
 *
 *    Description:  merge utilities
 *
 *        Version:  1.0
 *        Created:  22/12/2017 12:22:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include "merge_utls.hpp"
#include "ds.hpp"

using ds::ary;

int proc_noneoverlap_edges(aln_unit *u, aln_unit *v, GFAKluge &g, int bk_thres)
{
		edge_elem edg;
        edg.id = g.get_new_id(2);//new id here got to figure out how to do 
        
        group_elem grp;
        grp.id = g.get_new_id(3); // new id here
        
        //add '$' check here;
        

        edg.source_name = u->seq_id;
        edg.source_orientation_forward = u->isForward;
        
        edg.sink_name = v->seq_id;
        edg.sink_orientation_forward = v->isForward;
        //! if less than bk_threshold here, should it be extended? 1400 8                                         
        edg.source_begin = edg.source_end =  u->isForward ? u->seq_e : u->seq_s;
        
        edg.sink_begin = edg.sink_end = v->isForward ? v->seq_s : v->seq_e;
        

        //set ends here and adjust coordinate if orientation is backward
       //int source_add_one = 1;
       edg.ends.set(0,0);
       edg.ends.set(1,0);
       //if (u->isConverted) {
            //int64_t t = edg.source_begin;
            //edg.source_end = edg.source_begin = u->seq_len - u->seq_s; //- source_add_one;
            //edg.source_end = u->seq_len - t; // - source_add_one;
       //} 
		if (edg.source_end + bk_thres >  (uint64_t) u->seq_len ) {
			edg.source_end = edg.source_begin = u->seq_len;
			edg.ends.set(0,1);
			edg.ends.set(1,1);
			
		} //source_add_one = 0;} 
       //int sink_add_one = 1;
	   if (edg.source_begin < (uint64_t) bk_thres) edg.source_end = edg.source_begin = 0;
       edg.ends.set(2,0);
       edg.ends.set(3,0);
       //if (v->isConverted) {
            //int64_t t = edg.sink_begin;
            //edg.sink_end = edg.sink_begin = v->seq_len - v->seq_e;// -sink_add_one;
            //edg.sink_end = v->seq_len;// - sink_add_one;
       //}  
		if (edg.sink_end + bk_thres > (uint64_t) v->seq_len) { // this is impossible
			edg.sink_begin = edg.sink_end = v->seq_len;
			edg.ends.set(2,1);
			edg.ends.set(3,1);
		}
       
	   if (edg.sink_begin < (uint64_t) bk_thres) edg.sink_end = edg.sink_begin = 0;
        int dist =  v->ref_e - v->ref_s + u->ref_e - u->ref_s;
        edg.alignment = to_string(dist) + "M";

         
        opt_elem o;
        o.key = "JN";
        o.type = "Z";
        o.val = grp.id;
        edg.tags[o.key] = o;
        g.add_edge(edg.source_name, edg);
            
        grp.ordered = true;
        grp.items.push_back(u->seq_id);
        grp.orientations.push_back(u->isForward);
        grp.items.push_back(edg.id);
        grp.orientations.push_back(true);
        grp.items.push_back(v->seq_id);
        grp.orientations.push_back(v->isForward);
        
        o.key = "SE";
        o.type = "B";
        o.val = "i";
        o.val += to_string(0) + "," + to_string(0);//problems
        grp.tags[o.key] = o;  
        o.key = "CG";
        o.type = "Z";
        o.val = to_string(dist) + "M";
        grp.tags[o.key] = o;  
        g.add_group(grp); 
        return NORMAL;
}

int proc_overlapped_edges(aln_unit *u, aln_unit *v, GFAKluge &g, int bk_thres)
{
    edge_elem edg;
    edg.id = g.get_new_id(2);//new id here got to figure out how to do 
    int s_d = v->ref_s - u->ref_s; // distance of start
    int e_d; //distance of ends 
    edg.source_name = u->seq_id;
    edg.source_orientation_forward = u->isForward;
    edg.sink_name = v->seq_id;
    edg.sink_orientation_forward = v->isForward;
    //if not contained 
    if (u->ref_e < v->ref_e)  {
        //fprintf(stderr, "N:%s\t%d\t%s\t%d\n", u->seq_id.c_str(), e_d, v->seq_id.c_str(), s_d);
        e_d = v->ref_e - u->ref_e;
        edg.source_begin = u->seq_s + s_d;//check if wrong to minus one here?  
        edg.source_end = u->seq_e;
         
        edg.sink_begin = v->seq_s; 
        edg.sink_end = v->seq_e - e_d;//!
        //int source_add_one = 1;
        edg.ends.set(0,0);
        edg.ends.set(1,0);
        
		edg.source_begin = u->isForward ? u->seq_s + s_d : u->seq_s;
		edg.source_end = u->isForward ? u->seq_e : u->seq_e - s_d;
		//if (u->isConverted) {
            //int64_t t = edg.source_begin;
            //edg.source_begin = u->seq_len - edg.source_end;// + source_add_one;
            //edg.source_end = u->seq_len - t; //+ source_add_one;
        //}

        if (edg.source_end + bk_thres > (uint64_t) u->seq_len) {
				edg.source_end = u->seq_len;
                edg.ends.set(1,1); 
        } 
		if (edg.source_begin < (uint64_t) bk_thres) edg.source_begin = 0;
        //int sink_add_one = 1;
        edg.ends.set(2,0);
        edg.ends.set(3,0);
		edg.sink_begin = v->isForward ? v->seq_s : v->seq_s + e_d;
		edg.sink_end = v->isForward ? v->seq_e - e_d : v->seq_e;
        //if (v->isConverted) {
            //int64_t t = edg.sink_begin;
            //edg.sink_begin = v->seq_len - edg.sink_end;// + sink_add_one;
            //edg.sink_end = v->seq_len - t;// + sink_add_one;
        //}
		if (edg.sink_end + bk_thres > (uint64_t) v->seq_len) {
			edg.sink_end = v->seq_len;
			edg.ends.set(3,1);
		}        
	   if (edg.sink_begin < (uint64_t) bk_thres) edg.sink_begin = 0;	
        edg.alignment = to_string(u->ref_e - v->ref_s + 1) + "M";  
        g.add_edge(edg.source_name, edg);
    } else {
		e_d = u->ref_e - v->ref_e;
		edg.source_begin = u->isForward ? u->seq_s + s_d : u->seq_e - e_d;
		edg.source_end = u->isForward ? u->seq_e - e_d : u->seq_e - s_d;
        //edg.source_begin = u->seq_s + s_d;//check if wrong  
        //edg.source_end = u->seq_e - e_d;
        //fprintf(stderr, "Y:%s\t%d\t%s\t%d\n", u->seq_id.c_str(), e_d, v->seq_id.c_str(), s_d);
        edg.sink_begin = v->seq_s; 
        edg.sink_end = v->seq_e;//!
        
        //int source_add_one = 1;
        edg.ends.set(0,0);
        edg.ends.set(1,0);
        //if (u->isConverted) {
                //source_add_one = 0;
            //int64_t t = edg.source_begin;
            //edg.source_begin = u->seq_len - edg.source_end; //+ source_add_one;
            //edg.source_end = u->seq_len - t;// + source_add_one;
        //}
        if (edg.source_end + bk_thres > (uint64_t) u->seq_len) {
			edg.source_end = u->seq_len;
			edg.ends.set(1,1);
		}
		if (edg.source_begin < (uint64_t) bk_thres) edg.source_begin = 0;
        //int sink_add_one = 1;
        edg.ends.set(2,0);
        edg.ends.set(3,0);
        //if (v->isConverted) {
            //int64_t t = edg.sink_begin;
            //edg.sink_begin = v->seq_len - edg.sink_end;
            //edg.sink_end = v->seq_len - t;
        //} 
		if (edg.sink_end + bk_thres > (uint64_t) v->seq_len) {
			edg.sink_end = v->seq_len;
			edg.ends.set(3, 1);
		}
        
		if (edg.sink_begin < (uint64_t) bk_thres) edg.sink_begin = 0;
		edg.alignment = to_string(v->ref_e - v->ref_s) + "M";//!!!whether should add one here according to coordinate  
        g.add_edge(edg.source_name, edg);
    }
    return NORMAL;
}

//if two alignments should be one? u with smaller reference start alignment position
bool being_connected(aln_unit *u, aln_unit *v, int thres)
{
	bool direct = u->isForward;
	//forward
	int d_r = v->ref_s - u->ref_e;
	int d_q;
	if (direct) {
		d_q = v->seq_s - u->seq_e;
	} else 
		d_q = u->seq_s - v->seq_e;
	if (d_r < thres && d_q < thres) return true;
	else return false;
}
int gen_brks(aln_unit *u, int b_s, GFAKluge &g, int bk_thres, int *bk_count)
{
	if (b_s > bk_thres) {
		opt_elem o;
		o.key = "BK";
		o.type = "B";
		o.val = "i";
		o.val += to_string(b_s) + "," + to_string(b_s -1);
		g.add_tag(u->seq_id, &o);
		//sequence_elem& s = ss[au[i].seq_id];//if key exist?
		//cout<<s.to_string_1()<<endl;
		++*bk_count;
	}
	int b_e = u->isForward ? u->seq_e : u->seq_s;
	if ((u->seq_len - b_e + 1) > bk_thres) {
		opt_elem o;
		o.key = "BK";
		o.type = "B";
		o.val = "i";
		o.val += to_string(b_e - 1) + "," + to_string(b_e);
		g.add_tag(u->seq_id, &o);
		//sequence_elem& s = ss[au[i].seq_id];//if key exist
		//s.opt_fields.push_back(o);
		++*bk_count;
	}

	return NORMAL;
}
int gen_jnts(aln_unit *u, aln_unit *v, GFAKluge &g, int bk_thres, int *jn_count)
{
	if (u->ref_e < v->ref_s) {
		proc_noneoverlap_edges(u, v, g, bk_thres);
		++*jn_count;
	} else {
		proc_overlapped_edges(u, v, g, bk_thres); 
	}
	return NORMAL;
}

//int group(aln_unit *a, int n, ary<int> &ind, int thres, ary<int> &grp_ind)
//{
	//grp_ind.push(0);

	//for (int i = 0; i < n ; ++i) {
		//int u = ind[i];
		//int v = ind[i+1];
		//if (!being_connected(a + u, a + v, thres)) {
			//grp_ind.push(i+1);
		//}		
	//}		
	//grp_ind.push(n);
//}
//process ref to ref alignment
int proc_blk_r2r(aln_block *abk, GFAKluge &g, int bk_thres, int map_thres, int* bk_count, int* jn_count)
{
    //int dist;
    abk->sort_1();
	abk->set_is_covered();
	int u_size = abk->merge_alns(bk_thres, map_thres);
	aln_unit *au = abk->alns;
	//ary<int> ind, grp_ind;
	//int k = 0;	
	//int k = 0;
	//cerr<<abk->ref_id<<"\t"<<abk->ref_len<<endl;
	for (int i = 0; i < u_size; ++i)  cerr<< au[i].seq_id<<"\t"<<au[i].ref_s<<"\t"<<au[i].ref_e<<"\t"<<au[i].seq_s<<"\t"<<au[i].seq_e<<"\t"<<au[i].isCovered<<"\t"<<au[i].seq_e - au[i].seq_s<<endl;
	int dist;
	for (int i = 0; i < u_size; ++i) {
        //check left side
        //if (!au[i].isCovered) {
			//ind[k++] = i;
			if (au[i].seq_s > bk_thres) {
				opt_elem o;
				o.key = "BK";
				o.type = "B";
				o.val = "i";
				o.val += to_string(au[i].seq_s) + "," + to_string(au[i].seq_s -1);
				g.add_tag(au[i].seq_id, &o);
				//sequence_elem& s = ss[au[i].seq_id];//if key exist?
				//cout<<s.to_string_1()<<endl;
				++*bk_count;
			}
			if ((dist = au[i].seq_len - au[i].seq_e + 1) > bk_thres) {
				opt_elem o;
				o.key = "BK";
				o.type = "B";
				o.val = "i";
				o.val += to_string(au[i].seq_e - 1) + "," + to_string(au[i].seq_e);
				g.add_tag(au[i].seq_id, &o);
				//sequence_elem& s = ss[au[i].seq_id];//if key exist
				//s.opt_fields.push_back(o);
				++*bk_count;
			}
		//}
    } 
	for (int i = 0; i < u_size - 1; ++i ) {
		gen_jnts(au+i, au+i+1, g, bk_thres, jn_count);
	}
	 //= group(au, ind , k, );

	//for (int i = 0; i < u_size; ++i) if (!au[i].isCovered) ind[k++] = i; 
	//int b_s = au[ind[0]].isForward ? au[ind[0]].seq_s : au[ind[0]].seq_e;
    //for (int i = 0; i < k - 1; ++i) {
		//search for uncovered aln
		//int u = ind[i]; 
		//int v = ind[i+1];
		//if (au[u].seq_id == au[v].seq_id) {
			//if (au[u].isForward == au[v].isForward) {
				//if (!being_connected(au+u, au+v, bk_thres)) {
					//gen_brks(au+u, b_s, g, bk_thres, bk_count);
					//gen_jnts(au+u, au+v, g, bk_thres, jn_count);// local repetitive pattern or misalign should be reversed alignment  
					//b_s = au[v].isForward ? au[v].seq_s : au[v].seq_e;
				//}		
			//} else  {
				//gen_brks(au+u, b_s, g, bk_thres, bk_count);
				//gen_jnts(au+u, au+v, g, bk_thres, jn_count);// local repetitive pattern or misalign should be reversed alignment  
				//b_s = au[v].isForward ? au[v].seq_s : au[v].seq_e;
			//}
		//} else {
			//gen_brks(au+u, b_s, g, bk_thres, bk_count);
			//if (!g.edge_exist(au[u].seq_id, au[v].seq_id)) {
				//gen_jnts(au+u, au+v, g, bk_thres, jn_count);
				//b_s = au[v].isForward ? au[v].seq_s : au[v].seq_e;
			//}
		//}
	//}	
	//gen_brks(au + ind[k-1], b_s, g, bk_thres, bk_count);
    return NORMAL;
}

int proc_blk_q2r(aln_block *abk, int u_size, GFAKluge &g, int bk_thres, int map_thres, int* bk_count, int* jn_count)
{
	return 0;	
}


int proc_blk(int type, aln_block *abk, int u_size, GFAKluge &g, int bk_thres,int maplen_thres, int* bk_count, int* jn_count)
{
	if (type == 1) {
		return	proc_blk_r2r(abk, g, bk_thres, maplen_thres, bk_count, jn_count);
	} else {
		return	proc_blk_q2r(abk, u_size, g, bk_thres, maplen_thres, bk_count, jn_count);
	} 
}
/*
int proc_blk(aln_block *abk, int u_size, GFAKluge &g, int bk_thres, int* bk_count, int* jn_count)
{
    int dist;
    aln_unit *au = abk->alns;
    //map<string, sequence_elem, custom_key> ss = g.get_name_to_seq();
    for (int i = 0; i < u_size; ++i) {
        //check left side
        if (au[i].seq_s > bk_thres) {
            opt_elem o;
            o.key = "BK";
            o.type = "B";
            o.val = "i";
            o.val += to_string(au[i].seq_s) + ":" + to_string(au[i].seq_s -1);
            g.add_tag(au[i].seq_id, &o);
			//sequence_elem& s = ss[au[i].seq_id];//if key exist?
			//cout<<s.to_string_1()<<endl;
			++*bk_count;
        }
        if ((dist = au[i].seq_len - au[i].seq_e + 1) > bk_thres) {
            opt_elem o;
            o.key = "BK";
            o.type = "B";
            o.val = "i";
            o.val += to_string(au[i].seq_e - 1) + ":" + to_string(au[i].seq_e);
            g.add_tag(au[i].seq_id, &o);
            //sequence_elem& s = ss[au[i].seq_id];//if key exist
            //s.opt_fields.push_back(o);
			++*bk_count;
        }
    } 
    for (int i = 0; i < u_size; ++i) {
        for (int j = i + 1; j < u_size; ++j) {
            if (!g.edge_exist(au[i].seq_id, au[j].seq_id)) {
                int u,v;// u is smaller one
                if (au[i].ref_s < au[j].ref_s) {
                    u = i;
                    v = j;
                } else {
                    u = j;
                    v = i;
                }
                
                if (au[u].ref_e < au[v].ref_s) {
                    proc_noneoverlap_edges(au + u, au + v, g, bk_thres);
					++*jn_count;
                } else {
                    proc_overlapped_edges(au + u, au + v, g, bk_thres); 
                }
            }
        }
    }
    return NORMAL;
}
*/
