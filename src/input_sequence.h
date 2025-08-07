#ifndef INPUT_SEQUENCE_H
#define INPUT_SEQUENCE_H

#include <iostream>
#include <cstring>
#include <string>

struct Input_Sequence
{
	int seq_id;
    char* identifier;
    char* data;
};

struct Sequence_new
{
	int seq_id;
    const char* data;
	Sequence_new() {}
	Sequence_new(int id, const char* d) : seq_id(id), data(d) {}
//    char* identifier = null;
};
#endif
