/*
 *
 *
 */

#ifndef TRIGGERMAP_H
#define TRIGGERMAP_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <stack>

#include "LogiCalc.h"

class TriggerMap
{
public:
	TriggerMap(){};
	virtual ~TriggerMap(){};
	void MakeTable(std::string &);
	bool LookUp(uint32_t);
	void Dump();
protected:
private:
	bool ExtractBit(uint32_t, int);
	std::vector<bool> fLut;
	unsigned int fNsignal = 0;
	unsigned int fMapSize = 0;
};

#endif // TRIGGERMAP_H