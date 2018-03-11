/*
 * =====================================================================================
 *
 *       Filename:  ds.hpp
 *
 *    Description:  personal data structure 
 *
 *        Version:  1.0
 *        Created:  12/02/2018 10:48:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#ifndef _DATA_STRUCTURE_H
#define _DATA_STRUCTURE_H
#include <iostream>
#include <string.h>
#include "status_code.hpp"
using namespace std;
namespace ds {
	//array 
	template <typename T> 
	struct ary{
		T *data;
		int size; // element limitation
		int max; 
		int ptr; //point to next empty element
		ary() {
			data = NULL;
			int n = 1024;
			while (!data) {
				data = new(std::nothrow) T[n];
				n >>= 1;
			}	
			size = n << 1; // will be zero if data is null
			ptr = 0;
			max = -1;
		}
		//expand
		int ary_extend() {
			int n;
			T *n_d;
			if (size < (1<<26))
			   n = size << 1; 	
			else 
				n = n + 1024;	
			try {
				n_d = new T[n];		
			} catch (std::bad_alloc &) {
				throw RES_ERR;
			}
			memcpy(n_d, data, sizeof(T) * size);
			size = n;
			delete []data;
			data = n_d;
		}
		T& operator[] (int i) {
			if (data) {
				if (i < size) {
					ptr = i;
					if (ptr > max) max = ptr;
					return data[ptr++];
				} else {
					try {
						ary_extend();
					} catch (int e) {
						fprintf(stderr, "Memory is not enough, exit...\n");
						delete []data; //release data before quit
						exit(1);
					}
					ptr = i;
					if (ptr > max) max = ptr;
					return data[ptr++];	
				} 
			} else {
				fprintf(stderr, "Memory is not enough, exit...\n");
				exit(1);
				//return data[i];// not working
			} 
		}
		int push(T e)
		{
			if (ptr > max) max = ptr;
			data[ptr++] = e;
			return NORMAL;
		}
		int size_1()
		{
			return max + 1;
		}
		~ary() {
			if (data) delete []data;
		}	
		
	};
	//struct dict {	};
	//struct 

}

#endif

