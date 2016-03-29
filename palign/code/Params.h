#ifndef _PARAMS_H
#define _PARAMS_H

#include "stdinc.h"
#include "Input.h"
//#include "PriorInfo.h"
//#include "ConstraintModule.h"
//#include "ScoreJoin.h"
//#include "PidThreshold.h"
//#include "ConstraintConservation.h"
//#include "ConstraintMarkov.h"
//#include "LocalGC.h"

class ConvergeCriterion {
public:
	ConvergeCriterion();
	virtual ~ConvergeCriterion();

	virtual bool isConverge(int numiters);

	//setExpirationTime(int seconds);
	virtual void setMaxIters(int iters);
private:
	int maxIters;
};

class Params {
public:
	Params(int argc, char **argv);
	virtual ~Params();
	static void printHelp();

	//virtual Nullset* getNullset(int index);
	virtual int getNumNullsets();

	//virtual ConstraintModule* getConstraintModule(string name);
	//virtual ConstraintModule* getConstraintModule(int index);
	//virtual int getNumConstraintModules();

	virtual void printParams();

	//variables
	Input *input;
	ConvergeCriterion *convergeCriterion;
	//ScoreJoin *scoreJoin; 

	//local-GC init
	//LocalGCInit *localGCInit;

	//list of indices that corresponds to index of ConstraintModule
	//bool normalModuleLst[ConstraintModule::NUM_TYPES];
	//bool initialModuleLst[ConstraintModule::NUM_TYPES];


	//weights
	//double verticalWeight;
	double weightV1;
	double weightV2;
	double weightH1;
	double weightH2;

	string vertProbMarkovianFilename;

protected:
	virtual void setDefaultVals();
	virtual void readArgv();
private:
	//vector<Nullset*> vectNullsets;
	//vector<ConstraintModule*> constraintModules;

	//input variables
	int markovOrder;
	unsigned int randomSeed;
	string fastaFilename;
	string bgFastaFilename;
	int nullsetSize;
	int maxIters;

	//pid
	//PidThreshold *pidThreshold;
	string pidFilename;
	double pidAlpha;

	int argc;
	char **argv;

};

#endif

