#include <ecf/ECF.h>
#include "FunctionMinEvalOp.h"
/**
 * \brief Artificial Bee Colony algorithm (see e.g. http://www.scholarpedia.org/article/Artificial_bee_colony_algorithm)
 * 
 * ABC algorithm accepts only a single FloatingPoint genotype (vector of real values).
 * Additionally, it adds the following genotype for algorithm implementation:
 * 		- trial: floatingPoint genotype serving as a generation counter for each individual
 *		- probability: floatingPoint genotype for calculating the probability of getting chosen for each individual
 */

class MyAlg : public Algorithm
{
protected:
        // declare all available selection operators (not all get used)

        SelRandomOpP selRandomOp;
        SelBestOpP selBestOp;
		SelFitnessProportionalOpP selFitOp;
        
        uint limit;
		double ubound;
		double lbound;
public:
        
        MyAlg()
        {		// define algorithm name
				name_ = "MyAlg";

				// create selection operators needed
				// in this case, selRandomOp, SelBestOp and selFitOp
                selRandomOp = (SelRandomOpP) (new SelRandomOp);
                selBestOp = (SelBestOpP) (new SelBestOp);
				selFitOp = (SelFitnessProportionalOpP) (new SelFitnessProportionalOp);
        }

        //register any parameters
        void registerParameters(StateP state)
        {	
			// limit is a maximum number of cycles for each individual	
			registerParameter(state, "limit", (voidP) new uint(100), ECF::INT);              
        }

        
        bool initialize(StateP state)
		{		
			// initialize all operators
			selFitOp->initialize(state);
			selFitOp->setSelPressure(2);
			selBestOp->initialize(state);
			selRandomOp->initialize(state);
			
			voidP limit_ = getParameterValue(state, "limit");
			limit = *((uint*) limit_.get());

			voidP lBound = state->getGenotypes()[0]->getParameterValue(state, "lbound");
			lbound = *((double*) lBound.get());
			voidP uBound = state->getGenotypes()[0]->getParameterValue(state, "ubound");
			ubound = *((double*) uBound.get());

		// algorithm accepts a single FloatingPoint Genotype
			FloatingPointP flp (new FloatingPoint::FloatingPoint);
			if(state->getGenotypes()[0]->getName() != flp->getName()) {
				ECF_LOG_ERROR(state, "Error: ABC algorithm accepts only a FloatingPoint genotype!");
				throw ("");
			}

			FloatingPointP flpoint[3];
			for(uint iGen = 1; iGen < 3; iGen++) {

				flpoint[iGen] = (FloatingPointP) new FloatingPoint::FloatingPoint;
				state->setGenotype(flpoint[iGen]);

				flpoint[iGen]->setParameterValue(state, "dimension", (voidP) new uint(1));					

				// initial value of trial parameter should be (as close as possible to) 0				
				flpoint[iGen]->setParameterValue(state, "lbound", (voidP) new double(0));
				flpoint[iGen]->setParameterValue(state, "ubound", (voidP) new double(0.01));
				
			}
			ECF_LOG(state, 1, "ABC algorithm: added 2 FloatingPoint genotypes (trial, probability)");
 
            return true;
        }

       
        bool advanceGeneration(StateP state, DemeP deme)
        {		
//			In ABC, a colony of artificial bees search for artificial food sources (good solutions for a given problem):
//                , a colony of artificial bees contains three groups of bees: employed bees, onlooker bees, and scout bees 
//			- employed bees are associated with specific food sources
//			- onlooker bees watch the dance of employed bees within the hive to choose a food source
//			- scout bees search for food sources randomly
//
//			Initially, a randomly distributed initial population (food sources) is generated
//			REPEAT 
//			1)	Employed Bees Phase
//				1a)	for each food source createNewFoodSource() is called			
//
//			2)	Onlooker Bees Phase
//				2a)	each onlooker bee chooses a food source depending on their fitness values ( better individuals are more likely to be chosen)
//				2b) for each chosen food source createNewFoodSource() is called	
//	
//			3)	Scout Bees Phase
//				** as a rule, none or one food source gets abandoned per generation
//				3a) for each food source get trial variable
//				3b) remember the source if its trial exceeded limit 
//				3c) if there is a source that exceeded limit, replace it with a random food source
//
//			UNTIL(requirements are met)
//
//			*createNewFoodSource()
//				a)	for each food source find a neighbour (a random food source in the population) 
//				b)	produce a modification on the food source (discover a new food source)
//				c)	evaluate new food source 
//				d)	if the fitness value of the new one is better than that of the original source,
//						memorize the new source, forget the old one and set trial to 0
//					otherwise keep the old one and increment trial


              employedBeesPhase(state, deme);
			  onlookerBeesPhase(state, deme);	
			  scoutBeesPhase(state, deme);
              return true;
        }

		 bool employedBeesPhase(StateP state, DemeP deme)
        {	
			for( uint i = 0; i < deme->getSize(); i++ ) { // for each food source
				IndividualP food = deme->at(i);
				createNewFoodSource(food, state, deme);
			}
            return true;
        }

		 bool onlookerBeesPhase(StateP state, DemeP deme){
			calculateProbabilities(state, deme);
			int demeSize = deme->getSize();
			int i = 0;
			int n = 0;
			while( n < demeSize){
				int fact = i++ % demeSize;
				IndividualP food = deme->at(fact);
				FloatingPointP flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (food->getGenotype(2));
				double &probability = flp->realValue[0];
				
				if ( state->getRandomizer()->getRandomDouble() < probability){
					n++;
					createNewFoodSource(food, state, deme);
				}
				
					
			}
			 
			//for( uint i = 0; i < deme->getSize(); i++ ) { // for each food source
			//	//choose a food source depending on it's fitness value ( better individuals are more likely to be chosen)
			//	IndividualP food = selFitOp->select(*deme);
			//	createNewFoodSource(food, state, deme);
			//}
			return true;
		 }

		 bool scoutBeesPhase(StateP state, DemeP deme){
			IndividualP unimproved ;

			double maxTrial = 0;
			for( uint i = 0; i < deme->getSize(); i++ ) { // for each food source
				IndividualP food = deme->at(i);
				//get food source's trial variable
				FloatingPointP flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (food->getGenotype(1));
				double &trial = flp->realValue[0];
				
				//remember the source if its trial exceeded limit 
				if (trial > limit && trial >maxTrial){
					unimproved = food;
					maxTrial = trial;
				}					
			}

			//if there is a  food source that exceeded the limit, replace it with a random one
			if (unimproved != NULL){
					FloatingPointP flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (unimproved->getGenotype(1));
					double &trial = flp->realValue[0];
					trial = 0;
					flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (unimproved->getGenotype(0));
					flp->initialize(state);
					evaluate(unimproved);
			}

			return true;
		 }

		 bool createNewFoodSource(IndividualP food, StateP state, DemeP deme)
		 {  
			 //for each food source find a neighbour 
			IndividualP neighbour;
			do{
				neighbour = selRandomOp->select(*deme);
			}while(food->index == neighbour->index);

			//potential new food source
			IndividualP newFood = copy(food);

			FloatingPointP flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (food->getGenotype(0));
			std::vector< double > &foodVars = flp->realValue;
			flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (neighbour->getGenotype(0));
			std::vector< double > &neighbourVars = flp->realValue;
			flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (newFood->getGenotype(0));
			std::vector< double > &newFoodVars = flp->realValue;


			uint param = state->getRandomizer()->getRandomInteger((int)foodVars.size());
			double factor = state->getRandomizer()->getRandomDouble();
			double value = foodVars[param] * (1-2*factor)*(foodVars[param]-neighbourVars[param]);
			if (value > ubound)
				value = ubound;
			else if (value <lbound)
				value = lbound;

			//produce a modification on the food source (discover a new food source)
			newFoodVars[param] = value;
			evaluate(newFood);

			flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (food->getGenotype(1));
			double &foodTrial = flp->realValue[0];

//			d)	if the fitness value of the new food source is better than that of the original source,
//					memorize the new source, forget the old one and set trial to 0
//					otherwise keep the old one and increment trial
			if(newFood->fitness->isBetterThan( food->fitness) )
			{
				foodVars[param] = value;
				evaluate(food);
				foodTrial = 0;
			}
			else {
				foodTrial +=1;
			}
			return true;
		}
		bool calculateProbabilities(StateP state, DemeP deme){
			IndividualP bestFood =  selBestOp->select(*deme);
			double bestFitness = bestFood->fitness->getValue();

			for( uint i = 0; i < deme->getSize(); i++ ) { // for each food source
				IndividualP food = deme->at(i);
				double thisFitness = food->fitness->getValue();
				FloatingPointP flp = boost::dynamic_pointer_cast<FloatingPoint::FloatingPoint> (food->getGenotype(2));
				double &probability = flp->realValue[0];
				
				if (bestFitness == thisFitness)
					probability = 1.0;
				else if (thisFitness < bestFitness) //for maximization problems
						probability = 0.1 + 0.9 * thisFitness/bestFitness;
				else								//for minimization problems
						probability = 0.1 + 0.9 * bestFitness/thisFitness;
			}
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
