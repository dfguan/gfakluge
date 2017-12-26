/*
 * =====================================================================================
 *
 *       Filename:  merge_opt.cpp
 *
 *    Description:  parser of merge options
 *
 *        Version:  1.0
 *        Created:  22/12/2017 11:31:42
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "merge_opt.hpp"
#include "status_code.hpp"

struct option long_options[] =
{
    {"break_threshold", 1, NULL, 'b'},
    {"help", 0, NULL, 'h'},
    {0,0,0,0}
};

int help()
{
    fprintf(stderr, "gfa_merge [-b] <GFA_FILE> <PAF_FILE>\n");
    fprintf(stderr, "options:\n");
    fprintf(stderr, "       -b     break point threshold[100]\n");
    fprintf(stderr, "       -h     print help information\n");
    return NORMAL;
}



bool get_params(int argc ,char *argv[], opts *o) 
{
    int c,option_index = 0;
    optind = 1;
    while ((c = getopt_long(argc, argv, "hb:", long_options, &option_index)) > 0) {
        switch (c){
            case 'b':
                o->bk_thres = atoi(optarg);
                break;
            case 'h':
               help();
               return false;
            default:
               fprintf(stderr, "Undefined Parameters %c\n", c);
               help();
               return false; 
        }
    }
    if (optind + 2 > argc) {
        fprintf(stderr, "[option_parse]: arguments can't be omitted\n");
        
        help(); 
        return false;
    } else {
        ++optind;
        o->gfa_fn = argv[optind++];
        o->map_fn = argv[optind];
    }
    return true;
}

