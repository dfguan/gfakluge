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
int proc_noneoverlap_edges(aln_unit *u, aln_unit *v, GFAKluge &g)
{
        edge_elem edg;
        edg.id = g.get_new_id(2);//new id here got to figure out how to do 
        
        group_elem grp;
        grp.id = g.get_new_id(3); // new id here
        
        //add '$' check here;
        

        edg.source_name = u->seq_id;
        edg.source_orientation_forward = u->isConverted;
        
        edg.sink_name = v->seq_id;
        edg.sink_orientation_forward = v->isConverted;
        //! if less than bk_threshold here, should it be extended? 1400 8                                         
        edg.source_begin =  u->seq_e;
        edg.source_end =  u->seq_e;
        
        edg.sink_begin = v->seq_s;
        edg.sink_end = v->seq_s;
        

        //set ends here and adjust coordinate if orientation is backward
       //int source_add_one = 1;
       edg.ends.set(0,0);
       edg.ends.set(1,0);
       if (u->isConverted) {
            int64_t t = edg.source_begin;
            edg.source_begin = u->seq_len - edg.source_end; //- source_add_one;
            edg.source_end = u->seq_len - t; // - source_add_one;
       } 
		if (edg.source_end  == (uint64_t) u->seq_len) {edg.ends.set(0,1);edg.ends.set(1,1);} //source_add_one = 0;} 
       //int sink_add_one = 1;
       edg.ends.set(2,0);
       edg.ends.set(3,0);
       if (v->isConverted) {
            int64_t t = edg.sink_begin;
            edg.sink_begin = v->seq_len - edg.sink_end;// -sink_add_one;
            edg.sink_end = v->seq_len - t;// - sink_add_one;
       }  
		if (edg.sink_end == (uint64_t) v->seq_len) {
			edg.ends.set(2,1);
			edg.ends.set(3,1);
		}
       
        int dist =  edg.source_end - edg.source_begin + edg.sink_end - edg.sink_begin;
        edg.alignment = to_string(dist) + "M";

         
        opt_elem o;
        o.key = "JN";
        o.type = "Z";
        o.val = grp.id;
        edg.tags[o.key] = o;
        g.add_edge(edg.source_name, edg);
            
        grp.ordered = true;
        grp.items.push_back(u->seq_id);
        grp.orientations.push_back(u->isConverted);
        grp.items.push_back(edg.id);
        grp.orientations.push_back(true);
        grp.items.push_back(v->seq_id);
        grp.orientations.push_back(v->isConverted);
        
        o.key = "SE";
        o.type = "B";
        o.val = "i";
        o.val += to_string(u->seq_s) + "," + to_string(v->seq_e);
        grp.tags[o.key] = o;  
        o.key = "CG";
        o.type = "Z";
        o.val = to_string(dist) + "M";
        grp.tags[o.key] = o;  
        g.add_group(grp); 
        return NORMAL;
}

int proc_overlapped_edges(aln_unit *u, aln_unit *v, GFAKluge &g)
{
    edge_elem edg;
    edg.id = "*";//new id here got to figure out how to do 
    int s_d = v->ref_s - u->ref_s; // distance of start
    int e_d; //distance of ends 
    edg.source_name = u->seq_id;
    edg.source_orientation_forward = u->isConverted;
    edg.sink_name = v->seq_id;
    edg.sink_orientation_forward = v->isConverted;
    //if not contained 
    if (u->ref_e < v->ref_e)  {
        e_d = v->ref_e - u->ref_e;
        edg.source_begin = u->seq_e - s_d - 1;//check if wrong to minus one here?  
        edg.source_end = u->seq_e;
         
        edg.sink_begin = v->seq_s; 
        edg.sink_end = v->seq_e - e_d;//!
        //int source_add_one = 1;
        edg.ends.set(0,0);
        edg.ends.set(1,0);
        if (u->isConverted) {
            int64_t t = edg.source_begin;
            edg.source_begin = u->seq_len - edg.source_end;// + source_add_one;
            edg.source_end = u->seq_len - t; //+ source_add_one;
        }

        if (edg.source_end == (uint64_t) u->seq_len) {
                edg.ends.set(1,1); 
        } 
        //int sink_add_one = 1;
        edg.ends.set(2,0);
        edg.ends.set(3,0);
        if (v->isConverted) {
            int64_t t = edg.sink_begin;
            edg.sink_begin = v->seq_len - edg.sink_end;// + sink_add_one;
            edg.sink_end = v->seq_len - t;// + sink_add_one;
        }
		if (edg.sink_end == (uint64_t) v->seq_len) edg.ends.set(3,1);
         
        edg.alignment = to_string(u->ref_e - v->ref_e + 1) + "M";  
        g.add_edge(edg.source_name, edg);
    } else {
        e_d = u->ref_e - v->ref_e;
        edg.source_begin = u->seq_s;//check if wrong  
        edg.source_end = u->seq_e;
         
        edg.sink_begin = v->seq_s + s_d; 
        edg.sink_end = v->seq_e - e_d;//!
        
        //int source_add_one = 1;
        edg.ends.set(0,0);
        edg.ends.set(1,0);
        if (u->isConverted) {
                //source_add_one = 0;
            int64_t t = edg.source_begin;
            edg.source_begin = u->seq_len - edg.source_end; //+ source_add_one;
            edg.source_end = u->seq_len - t;// + source_add_one;
        }
        if (edg.source_end == (uint64_t) u->seq_len) edg.ends.set(1,1);

        //int sink_add_one = 1;
        edg.ends.set(2,0);
        edg.ends.set(3,0);
        if (v->isConverted) {
            int64_t t = edg.sink_begin;
            edg.sink_begin = v->seq_len - edg.sink_end;
            edg.sink_end = v->seq_len - t;
        } 
		if (edg.sink_end == (uint64_t) v->seq_len) edg.ends.set(3, 1);
        
		edg.alignment = to_string(v->ref_e - v->ref_s) + "M";//!!!whether should add one here according to coordinate  
        g.add_edge(edg.source_name, edg);
    }
    return NORMAL;
}

int proc_blk(aln_block *abk, int u_size, GFAKluge &g, int bk_thres)
{
    int dist;
    aln_unit *au = abk->alns;
    map<string, sequence_elem, custom_key> ss = g.get_name_to_seq();
    for (int i = 0; i < u_size; ++i) {
        //check left side
        if (au[i].seq_s > bk_thres) {
            opt_elem o;
            o.key = "BK";
            o.type = "B";
            o.val = "i";
            o.val += to_string(au[i].seq_s) + ":" + to_string(au[i].seq_s -1);
            sequence_elem& s = ss[au[i].seq_id];//not temporary?
            s.opt_fields.push_back(o);
        }
        if ((dist = au[i].seq_len - au[i].seq_e) > bk_thres) {
            opt_elem o;
            o.key = "BK";
            o.type = "B";
            o.val = "i";
            o.val += to_string(au[i].seq_e - 1) + ":" + to_string(au[i].seq_e);
            sequence_elem& s = ss[au[i].seq_id];//not temporary?
            s.opt_fields.push_back(o);
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
                    proc_noneoverlap_edges(au + u, au + v, g);
                } else {
                    proc_overlapped_edges(au + u, au + v, g); 
                }
            }
        }
    }
    return NORMAL;
}
