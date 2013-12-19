#include <ecf/ECF.h>
#include "FunctionMinEvalOp.h"
/**
 * \brief Clonal Selection Algorithm (see e.g. http://en.wikipedia.org/wiki/Clonal_Selection_Algorithm)
 * 
 * this CLONALG implements: - static cloning :  n of the best antibodies are cloned beta time, making the size of the clones population  equal n*beta
 *							- proportional cloning: number of clones per antibody is proportional to that ab's fitness
 *							- inversely proportional hypermutation : better antibodies are mutated less
 *							- CLONALG1 - at new generation each antibody will be substituded by the best individual of its set of beta*population clones
 *							- CLONALG2 - new generation will be formed by the best (1-d)*populationSize clones ( or all if the number of clones is less than that )- birthPhase where d * populationSize of new antibodies are randomly created and added to the population
 *							
 * CLONALG algorithm accepts only a single FloatingPoint genotype (vector of real values).
 */
class MyAlg : public Algorithm
{
protected:
        
		double ubound;
		double lbound;
		uint dimension;

		uint n;					// number of antibodies cloned every generation
		double beta;			// parameter which determines the number of clones for every antibody
		double c;				// da uzmem iz postojeceg mutation parametra mutation parametar
		double d;				// fraction of population regenerated every generation
		string cloningVersion;	// specifies whether to use static or proportional cloning
		string selectionScheme;	// specifies which selection scheme to use CLONALG1 or CLONALG2

		// sort vector of antibodies in regards to their fitness
		static bool sortPopulationByFitness (IndividualP ab1,IndividualP ab2) { return ( ab1->fitness->isBetterThan(ab2->fitness)); }

		// sort vector of antibodies first by their antibody parents and then to their fitness
		static bool sortPopulationByParentAndFitness (IndividualP ab1,IndividualP ab2) 
		{ 
			FloatingPointP flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (ab1->getGenotype(1));
			double &parentAb1 = flp->realValue[0];
			flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (ab2->getGenotype(1));
			double &parentAb2 = flp->realValue[0];

			if (parentAb1 <parentAb2) return true;
			if (parentAb1 == parentAb2) return (ab1->fitness->isBetterThan(ab2->fitness));
	
			return false;
		}

public:
        
        MyAlg()
        {		// define algorithm name
				name_ = "MyAlg";
        }

        //register any parameters
        void registerParameters(StateP state)
        {	
			registerParameter(state, "n", (voidP) new uint(50), ECF::INT);
			registerParameter(state, "beta", (voidP) new double(0.1), ECF::DOUBLE);
			registerParameter(state, "c", (voidP) new double(0.8), ECF::DOUBLE);
			registerParameter(state, "d", (voidP) new double(0.0), ECF::DOUBLE);
			registerParameter(state, "cloningVersion", (voidP) new string("static"), ECF::STRING);
			registerParameter(state, "selectionScheme", (voidP) new string("CLONALG2"), ECF::STRING);
		}

        
        bool initialize(StateP state)
		{		

			voidP lBound = state->getGenotypes()[0]->getParameterValue(state, "lbound");
			lbound = *((double*) lBound.get());

			voidP uBound = state->getGenotypes()[0]->getParameterValue(state, "ubound");
			ubound = *((double*) uBound.get());

			voidP dimension_ = state->getGenotypes()[0]->getParameterValue(state, "dimension");
			dimension = *((uint*) dimension_.get());
			
			
			voidP populationSize_ = state->getRegistry()->getEntry("population.size");
			uint populationSize = *((uint*) populationSize_.get());

			voidP n_ = getParameterValue(state, "n");
			n = *((uint*) n_.get());
			if( n<1 || n>populationSize) {
				ECF_LOG(state, 1, "Error: CLONALG requires parameter 'n' to be an integer in range < 0, population.size] ");
				throw "";}


			voidP beta_ = getParameterValue(state, "beta");
			beta = *((double*) beta_.get());
			if( beta <= 0 ) {
				ECF_LOG(state, 1, "Error: CLONALG requires parameter 'beta' to be a double greater than 0");
				throw "";}


			voidP c_ = getParameterValue(state, "c");
			c = *((double*) c_.get());
			if( c <= 0 ) {
				ECF_LOG(state, 1, "Error: CLONALG requires parameter 'c' to be a double greater than 0");
				throw "";}


			voidP d_ = getParameterValue(state, "d");
			d = *((double*) d_.get());
			if( d<0 || d>1 ) {
				ECF_LOG(state, 1, "Error: CLONALG requires parameter 'd' to be a double in range [ 0, 1] ");
				throw "";}

			voidP cloning_ = getParameterValue(state, "cloningVersion");
			cloningVersion = *((string*) cloning_.get());
			if( cloningVersion != "static" && cloningVersion != "proportional"  ) {
				ECF_LOG(state, 1, "Error: CLONALG requires parameter 'cloningVersion' to be either 'static' or a 'proportional'");
				throw "";}

			voidP selection_ = getParameterValue(state, "selectionScheme");
			selectionScheme = *((string*) selection_.get());
			if( selectionScheme != "CLONALG1" && selectionScheme != "CLONALG2"  ) {
				ECF_LOG(state, 1, "Error: CLONALG requires parameter 'selectionScheme' to be either 'CLONALG1' or 'CLONALG2'");
				throw "";}
						

		    // algorithm accepts a single FloatingPoint Genotype
			FloatingPointP flp (new FloatingPoint::FloatingPoint);
			if(state->getGenotypes()[0]->getName() != flp->getName()) {
				ECF_LOG_ERROR(state, "Error: CLONALG algorithm accepts only a FloatingPoint genotype!");
				throw ("");
			}
			//add parentAntibody genotype
			if( selectionScheme == "CLONALG1"){
				FloatingPointP flpoint[2];
				for(uint iGen = 1; iGen < 2; iGen++) {
					flpoint[iGen] = (FloatingPointP) new FloatingPoint::FloatingPoint;
					state->setGenotype(flpoint[iGen]);

					flpoint[iGen]->setParameterValue(state, "dimension", (voidP) new uint(1));					

					// initial value of age parameter should be (or as close as possible to) 0				
					flpoint[iGen]->setParameterValue(state, "lbound", (voidP) new double(0));
					flpoint[iGen]->setParameterValue(state, "ubound", (voidP) new double(0.01));
					
				}
				ECF_LOG(state, 1, "CLONALG algorithm: added 1 FloatingPoint genotype (parentAntibody)");
			}

            return true;
        }

       
        bool advanceGeneration(StateP state, DemeP deme)
        {	
			  std::vector<IndividualP> clones;
			  if (selectionScheme == "CLONALG1" && state->getGenerationNo()== 0)
				 markAntibodies(deme);
			  cloningPhase(state, deme, clones);
			  hypermutationPhase(state, deme, clones);
			  selectionPhase(state, deme, clones);
              birthPhase(state, deme, clones);
			  replacePopulation(state, deme, clones);
			 
              return true;
        }
		
		
		bool markAntibodies(DemeP deme){
			//mark antibodies so the alg can know which clone belongs to which parent Antibody
			for( uint i = 0; i < deme->getSize(); i++ ) { // for each antibody
			
					FloatingPointP flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (deme->at(i)->getGenotype(1));
					double &parentAb = flp->realValue[0];
					parentAb = i;				
			}
			return true;
		}

		bool cloningPhase(StateP state, DemeP deme, std::vector<IndividualP> &clones)
		{	
			// calculate number of clones per antibody
			uint clonesPerAntibody = beta * deme->getSize();

			// storing all antibodies in a vector
			for( uint i = 0; i < deme->getSize(); i++ ) { // for each antibody	
				
				IndividualP antibody = deme->at(i);
				clones.push_back(antibody);
			}

			// sorting all antibodies
			std::sort (clones.begin(), clones.end(), sortPopulationByFitness);
			
			// leaving n best antibodies for cloning
			clones.erase (clones.begin()+n,clones.end());

			
			for( uint i = 0; i < n; i++ ){ // for each antibody in clones vector
				IndividualP antibody = clones.at(i);
			
				//static cloning : cloning each antibody beta*populationSize times
				if (cloningVersion == "static"){
					for (uint j = 0; j < clonesPerAntibody; j++) 
						clones.push_back(copy(antibody));
				}

				//proportional cloning 
				else{ 
					uint scaledNumberofClones = clonesPerAntibody/(i+1);
					for (uint j = 0; j < scaledNumberofClones ; j++) 
						clones.push_back(copy(antibody));
				}
		    }
			
			return true;
		}

		bool hypermutationPhase(StateP state, DemeP deme, std::vector<IndividualP> &clones)
		{	
			// calculate number of clones per antibody
			uint clonesPerAntibody = beta * deme->getSize() +1;
			
			uint M;	// M - number of mutations of a single antibody 
			uint k;

			for( uint i = 0; i < clones.size(); i++ ){ // for each antibody in vector clones
				IndividualP antibody = clones.at(i);
				
				FloatingPointP flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (antibody->getGenotype(0));
			    std::vector< double > &antibodyVars = flp->realValue;
				
				// inversely proportional hypermutation (better antibodies are mutated less)
				k = i/clonesPerAntibody +1;
				M = (int) ((1- 1/(double)(k)) * (c*dimension) + (c*dimension));
				
				
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
				
				FitnessP originalFitness = antibody->fitness;
				evaluate(antibody);

			}		
			return true;
		}
		
		bool selectionPhase(StateP state, DemeP deme, std::vector<IndividualP> &clones)
		{	
			if( selectionScheme == "CLONALG1") {
				uint selNumber = (uint)((1-d)*deme->getSize());
				std::sort (clones.begin(), clones.end(), sortPopulationByParentAndFitness);
			
				int j=0;
				std::vector<IndividualP> temp_clones;

				for (uint i = 0; i < clones.size(); i++){
					FloatingPointP flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (clones.at(i)->getGenotype(1));
					double &parentAb = flp->realValue[0];
					if (j == parentAb && j < selNumber){
						temp_clones.push_back(clones.at(i));
						j++;
					}
				}
				clones = temp_clones;				
			}


			else{
				std::sort (clones.begin(), clones.end(), sortPopulationByFitness);
				uint selNumber = (uint)((1-d)*deme->getSize());

				//keep best (1-d)*populationSize antibodies ( or all if the number of clones is less than that )
				if(selNumber < clones.size())
					clones.erase (clones.begin()+ selNumber, clones.end());
			}

			return true;
		}
		
		bool birthPhase(StateP state, DemeP deme, std::vector<IndividualP> &clones)
		{	
			//  birthNumber - number of new antibodies randomly created and added 
			uint birthNumber = deme->getSize() - clones.size();
			//mark the newAntibody's paretntAb 
			uint mark = clones.size();


			IndividualP newAntibody = copy(deme->at(0));
			FloatingPointP flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (newAntibody->getGenotype(0));

			for (uint i = 0; i<birthNumber; i++){
				//create a random antibody
				flp->initialize(state);
				evaluate(newAntibody);

				if (selectionScheme == "CLONALG1"){
					//mark ab's parentAb
					flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (newAntibody->getGenotype(1));
					double &parentAb = flp->realValue[0];
					parentAb = mark++;
				}

				//add it to the clones vector
				clones.push_back(copy(newAntibody));
			}
			
			return true;
		}

		bool replacePopulation(StateP state, DemeP deme, std::vector<IndividualP> &clones)
		{
			//replace population with the contents of clones vector
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
