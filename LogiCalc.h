/*
 *
 *
 */
#ifndef LOGICALC_H
#define LOGICALC_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <stack>
#include <bitset>

#include <algorithm>

#include "infixtorp.cxx"
// #include "Trigger.cxx"

class LogiCalc
{
public:
	LogiCalc(){};
	virtual ~LogiCalc(){};
	std::vector<std::string> & SetFormula(std::string);
	std::string& GetFormula();
	int GetSigMax(){return fSigMax;};
	bool Calc(uint32_t);
	void Dump(){};
protected:
	bool ExtractBit(uint32_t, int);
	std::vector<std::string> fCommands;
	std::stack<bool> fStack;

private:
	//std::vector<std::string> &Split(std::string);
	int ComPushBack(std::string &);

	
	
	int fNsigMax = 32;
	int fSigMax = 0;
	const std::vector<std::string> fFuncs = {"&", "|", "!", "x", "d"};
	const std::string fSeparator = std::string(" ");

};

#endif // LOGICALC_H