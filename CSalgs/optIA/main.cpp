#include <ecf/ECF.h>
#include "FunctionMinEvalOp.h"
/**
 *\brief Optimization Immune  Algorithm (opt-IA) 
 * this opt-IA implements:  - static cloning : all antibodies are cloned dup times, making the size of the clone population equal dup*spoplationSize
 *							- inversely proportional hypermutation : better antibodies are mutated less
 *							- static pure aging - if an antibody exceeds tauB number of trials, it is replaced with a new randomly created antibody
 *							- birthPhase: if the number of antibodies that survive the aging Phase is less than populationSize, new randomly created abs are added to the population
 *							- optional elitism
 * opt-IA algorithm accepts only a single FloatingPoint genotype
 * Additionally, opt-IA adds a FloatingPoint genotype (age) 
 */
class MyAlg : public Algorithm
{
protected:
        
		double ubound;
		double lbound;
		uint dimension;

		uint dup;		// number of clones per antibody
		double c;		// mutation parameter
		double tauB;	// maximum number of generations without improvement 
		string elitism;	// specifies whether to use elitism or not

		// sort vector of antibodies in regards to their fitness
		static bool sortPopulationByFitness (IndividualP ab1,IndividualP ab2) { return ( ab1->fitness->isBetterThan(ab2->fitness)); }
public:
        
        MyAlg()
        {		// define algorithm name
				name_ = "MyAlg";
        }


		//register any parameters
        void registerParameters(StateP state)
        {	
			registerParameter(state, "dup", (voidP) new uint(10), ECF::INT);
			registerParameter(state, "c", (voidP) new double(0.2), ECF::DOUBLE);
			registerParameter(state, "tauB", (voidP) new double(100), ECF::DOUBLE);
			registerParameter(state, "elitism", (voidP) new string("false"), ECF::STRING);
		}


		bool initialize(StateP state)
		{		
			voidP lBound = state->getGenotypes()[0]->getParameterValue(state, "lbound");
			lbound = *((double*) lBound.get());

			voidP uBound = state->getGenotypes()[0]->getParameterValue(state, "ubound");
			ubound = *((double*) uBound.get());

			voidP dimension_ = state->getGenotypes()[0]->getParameterValue(state, "dimension");
			dimension = *((uint*) dimension_.get());

			voidP dup_ = getParameterValue(state, "dup");
			dup = *((uint*) dup_.get());
			if( *((int*) dup_.get()) <= 0 ) {
				ECF_LOG(state, 1, "Error: opt-IA requires parameter 'dup' to be an integer greater than 0");
				throw "";}

			voidP c_ = getParameterValue(state, "c");
			c = *((double*) c_.get());
			if( c <= 0 ) {
				ECF_LOG(state, 1, "Error: opt-IA requires parameter 'c' to be a double greater than 0");
				throw "";}

			voidP tauB_ = getParameterValue(state, "tauB");
			tauB = *((double*) tauB_.get());
			if( tauB < 0 ) {
				ECF_LOG(state, 1, "Error: opt-IA requires parameter 'tauB' to be a nonnegative double value");
				throw "";}

			voidP elitism_ = getParameterValue(state, "elitism");
			elitism = *((string*) elitism_.get());
			if( elitism != "true" && elitism != "false"  ) {
				ECF_LOG(state, 1,  "Error: opt-IA requires parameter 'elitism' to be either 'true' or 'false'");
				throw "";}


			// algorithm accepts a single FloatingPoint Genotype
			FloatingPointP flp (new FloatingPoint::FloatingPoint);
			if(state->getGenotypes()[0]->getName() != flp->getName()) {
				ECF_LOG_ERROR(state, "Error: opt-IA algorithm accepts only a FloatingPoint genotype!");
				throw ("");}

			// algorithm adds another FloatingPoint genotype (age)
			FloatingPointP flpoint[2];
			for(uint iGen = 1; iGen < 2; iGen++) {
				flpoint[iGen] = (FloatingPointP) new FloatingPoint::FloatingPoint;
				state->setGenotype(flpoint[iGen]);

				flpoint[iGen]->setParameterValue(state, "dimension", (voidP) new uint(1));					

				// initial value of age parameter should be (or as close as possible to) 0				
				flpoint[iGen]->setParameterValue(state, "lbound", (voidP) new double(0));
				flpoint[iGen]->setParameterValue(state, "ubound", (voidP) new double(0.01));
				
			}
			ECF_LOG(state, 1, "opt-IA algorithm: added 1 FloatingPoint genotype (antibody age)");
			
            return true;
		}


		bool advanceGeneration(StateP state, DemeP deme)
		{	
			std::vector<IndividualP> clones;
			 
			cloningPhase(state, deme, clones);
			hypermutationPhase(state, deme, clones);
			agingPhase(state, deme, clones);
			selectionPhase(state, deme, clones);
            birthPhase(state, deme, clones);
			replacePopulation(state, deme, clones);

			return true;
		}


		bool cloningPhase(StateP state, DemeP deme, std::vector<IndividualP> &clones)
		{
			// storing all antibodies in a vector
			for( uint i = 0; i < deme->getSize(); i++ )  // for each antibody	
				clones.push_back(deme->at(i));

			for( uint i = 0; i < deme->getSize(); i++ ){ // for each antibody in clones vector
				IndividualP antibody = clones.at(i);
				
				// static cloning is fitness independent : : cloning each antibody dup times
				for (uint j = 0; j < dup; j++) 
					clones.push_back(copy(antibody));					
			}

			return true;
		}


		bool hypermutationPhase(StateP state, DemeP deme, std::vector<IndividualP> &clones)
		{	
			uint M;	// M - number of mutations of a single antibody 
			uint k;

			//sort 
			std::sort (clones.begin(), clones.end(), sortPopulationByFitness);

			for( uint i = 0; i < clones.size(); i++ ){ // for each antibody in vector clones
				IndividualP antibody = clones.at(i);
				
				FloatingPointP flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (antibody->getGenotype(0));
			    std::vector< double > &antibodyVars = flp->realValue;
				
				k = 1 + i/(dup+1);
				M =(int) ((1- 1/(double)(k)) * (c*dimension) + (c*dimension));
				
				// mutate M times
				for (uint j = 0; j < M; j++){
					uint param = state->getRandomizer()->getRandomInteger((int)antibodyVars.size());
					
					double randDouble1 = state->getRandomizer()->getRandomDouble();
					double randDouble2 = state->getRandomizer()->getRandomDouble();
					double value = antibodyVars[param] + (1-2*randDouble1)* 0.2 *  (ubound - lbound) * pow(2, -16*randDouble2 );
					
					if (value > ubound)
						value = ubound;
					else if (value <lbound)
						value = lbound;

					//produce a mutation on the antibody 
					antibodyVars[param] = value;
				}
				FitnessP parentFitness = antibody->fitness;
				evaluate(antibody);

				// if the clone is better than its parent, reset clone's age
				if(antibody-> fitness->isBetterThan(parentFitness)){					
					flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (antibody->getGenotype(1));
					double &age = flp->realValue[0];
					age = 0;
				} 
			}
			return true;
		}


		bool agingPhase(StateP state, DemeP deme,  std::vector<IndividualP> &clones)
		{	
			//sort 
			std::sort (clones.begin(), clones.end(), sortPopulationByFitness);

			std::vector<IndividualP> temp_clones;

			// if elitism = true , preserve the best antibody regardless of its age
			if (elitism == "true")
				temp_clones.push_back(clones.at(0));

			for (uint i = 0; i < clones.size(); i++){// for each antibody
				IndividualP antibody = clones.at(i);

				//age each antibody
				FloatingPointP flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (antibody->getGenotype(1));
				double &age = flp->realValue[0];
				age += 1;
				
				// static aging: if an antibody exceeds tauB number of trials, it is replaced with a new randomly created antibody
				if (age <=tauB)
					temp_clones.push_back(antibody);		  
			}
			clones = temp_clones;
			return true;
		}

		bool selectionPhase(StateP state, DemeP deme, std::vector<IndividualP> &clones)
		{	
			//sort 
			std::sort (clones.begin(), clones.end(), sortPopulationByFitness);
		
			//keep best populationSize antibodies ( or all if the number of clones is less than that ), erase the rest
			if(clones.size() > deme->getSize())
				clones.erase (clones.begin()+ deme->getSize(), clones.end());

			return true;
		}

		bool birthPhase(StateP state, DemeP deme, std::vector<IndividualP> &clones)
		{
			//number of new antibodies (randomly created)
			uint birthNumber = deme->getSize() - clones.size();

			//if no new antibodies are needed, return (this if part is optional, code works fine w/o it)
			if (birthNumber == 0) return true;

			IndividualP newAntibody = copy(deme->at(0));
			FloatingPointP flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (newAntibody->getGenotype(0));
	
			for (uint i = 0; i<birthNumber; i++){
				//create a random antibody
				flp->initialize(state);
				evaluate(newAntibody);
			
				//reset its age
				flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (newAntibody->getGenotype(1));
				double &age = flp->realValue[0];
				age = 0;

				//add it to the clones vector
				clones.push_back(copy(newAntibody));
			}
			return true;
		}

		bool replacePopulation(StateP state, DemeP deme, std::vector<IndividualP> &clones)
		{
			//replace population with the contents of the clones vector
			for( uint i = 0; i < clones.size(); i++ ) // for each antibody
				deme->replace(i, clones.at(i));
			
			clones.clear();
			
			return true;
		}



};
typedef boost::shared_ptr<MyAlg> MyAlgP;



//this main() function optimizes a single COCO function (Id set in config file)
//function Ids: noiseless 1-24, noisy 101-130


//int main(int argc, char **argv)
//{
//	StateP state (new State);
//	//set newAlg
//	MyAlgP alg = (MyAlgP) new MyAlg;
//	state->addAlgorithm(alg);
//
//	// set the evaluation operator
//	state->setEvalOp(new FunctionMinEvalOp);
//	state->initialize(argc, argv);
//	state->run();
//	int i;
//	cin>>i;
//
//	return 0;
//}



//
// this main() function iterates over multiple COCO functions and optimizes each one in turn
// function Ids: noiseless 1-24, noisy 101-130
//

int main(int argc, char **argv)
{
	// run for selected COCO functions
	for(uint function = 1; function < 25; function++) {

		// read XML config
		std::ifstream fin(argv[1]);
		if (!fin) {
			throw std::string("Error opening file! ");
		}

		std::string xmlFile, temp;
		while (!fin.eof()) {
			getline(fin, temp);
			xmlFile += "\n" + temp;
		}
		fin.close();

		// set log and stats parameters
		std::string funcName = uint2str(function);
		std::string logName = "log", statsName = "stats";
		if(function < 10) {
			logName += "0";
			statsName += "0";
		}
		logName += uint2str(function) + ".txt";
		statsName += uint2str(function) + ".txt";

		// update in XML
		XMLResults results;
		XMLNode xConfig = XMLNode::parseString(xmlFile.c_str(), "ECF", &results);
		XMLNode registry = xConfig.getChildNode("Registry");

		XMLNode func = registry.getChildNodeWithAttribute("Entry", "key", "coco.function");
		func.updateText(funcName.c_str());
		XMLNode log = registry.getChildNodeWithAttribute("Entry", "key", "log.filename");
		log.updateText(logName.c_str());
		XMLNode stats = registry.getChildNodeWithAttribute("Entry", "key", "batch.statsfile");
		stats.updateText(statsName.c_str());

		// write back
		std::ofstream fout(argv[1]);
		fout << xConfig.createXMLString(true);
		fout.close();


		// finally, run ECF on single function
		StateP state (new State);

		//set newAlg
		MyAlgP alg = (MyAlgP) new MyAlg;
		state->addAlgorithm(alg);
		// set the evaluation operator
		state->setEvalOp(new FunctionMinEvalOp);

		state->initialize(argc, argv);
		state->run();
	}

	return 0;
}
