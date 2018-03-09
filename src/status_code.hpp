/*
 * =====================================================================================
 *
 *       Filename:  error.hpp
 *
 *    Description:  type of errors
 *
 *        Version:  1.0
 *        Created:  21/12/2017 14:34:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#ifndef _ERR_H
#define _ERR_H

//NORMAL STATUS
#define NORMAL 0
//ERRORs require to exit
#define RES_ERR 1 //resource error memory or cpu
#define IO_ERR 2 //file not found or some other ERR
#define FL_FORMAT_ERR 3
//WARNINGs 
#define EMPTY_FILE 4
#define NOT_ENOUGH 5
#define OTHER 6

#endif
