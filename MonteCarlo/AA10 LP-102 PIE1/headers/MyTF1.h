#ifndef MYTF1
#define MYTF1
#include "TF1.h"
//____________________________________________________________________

class MyTF1 : public TF1
{
public:
	using TF1::TF1;
	//using TF1::GetRandom;
	Double_t GetRandom();
};
//_________________________________________________________________


#endif
