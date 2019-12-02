#ifndef SOBOL
#define SOBOL
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF2.h"
//____________________________________________________________________

class MyTF1 : public TF1
{
public:
	// MyTF1() {};
	// MyTF1(const char *name, const char *formula,
	// 	Double_t from, Double_t to, const char *option) :
	// 	TF1(name, formula, from, to, option)
	// {};

	using TF1::TF1;
	using TF1::GetRandom;
	Double_t GetRandom(Double_t);
};
//_________________________________________________________________

class MyTH1F : public TH1F
{
public:

	using TH1F::TH1F;
	using TH1F::GetRandom;
	Double_t GetRandom(Double_t) const;
};
//_________________________________________________________________
class MyTH2F : public TH2F
{
public:

	using TH2F::TH2F;
	using TH2F::GetRandom2;
	void GetRandom2(Double_t&, Double_t&, Double_t);
};
//_________________________________________________________________
class MyTF2 : public TF2
{
public:

	using TF2::TF2;
	using TF2::GetRandom2;
	void GetRandom2(Double_t&, Double_t&, Double_t);
};
//_________________________________________________________________



#endif
