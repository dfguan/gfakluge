/*
 * =====================================================================================
 *
 *       Filename:  convert_gfa.hpp
 *
 *    Description:  head file of convert_gfa.cpp 
 *
 *        Version:  1.0
 *        Created:  30/12/2017 11:05:20
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#ifndef _CONVERT_GFA_H
#define _CONVERT_GFA_H



#include "gfakluge.hpp"
#include "status_code.hpp"


using namespace std;
using namespace gfak;


int	convert_gfa(GFAKluge& gg);
#endif
