#include "ca.h"
#include <vector>

namespace cadv {
	//int no_births=0;
	//int no_deaths=0;
	CellAut *CellAut::instance = NULL;

	//FUNCTIONS FOR Cell
	void Cell::inicNeigh(int n, int type){
		switch(type){
			case 0:
				no_diff_neigh = n;
				diff_neigh_used = 0;
				diff_neigh = new Cell* [n] ; 
				break;
			case 1:
				no_met_neigh = n;
				met_neigh_used = 0;
				met_neigh = new Cell* [n] ; 
				break;
			case 2:
				no_repl_neigh = n;
				repl_neigh_used = 0;
				repl_neigh = new Cell* [n] ; 
				claims = new double [n+1];
				claims[0] = par_claimEmpty;
				break;
		}
	}

	void Cell::setNeigh(class Cell *np, const int type, int which_neigh){
		switch(type){
			case 0:
				if(which_neigh < 0) which_neigh = diff_neigh_used;
				if(which_neigh >= no_diff_neigh){
					std::cerr << "WARNING: Cell::setNeigh: cant set incorrect diff neighbour " << which_neigh << " when number of neighbours is " << no_diff_neigh << std::endl;
				} else {
					diff_neigh[which_neigh] = np;
					++diff_neigh_used;
				}
				break;
			case 1:
				if(which_neigh < 0) which_neigh = met_neigh_used;
				if(which_neigh >= no_met_neigh){
					std::cerr << "WARNING: Cell::setNeigh: cant set incorrect met neighbour " << which_neigh << " when number of neighbours is " << no_met_neigh << std::endl;
				} else {
					met_neigh[which_neigh] = np;
					++met_neigh_used;
				}
				break;
			case 2:
				if(which_neigh < 0) which_neigh = repl_neigh_used;
				if(which_neigh >= no_repl_neigh){
					std::cerr << "WARNING: Cell::setNeigh: cant set incorrect repl neighbour " << which_neigh << " when number of neighbours is " << no_repl_neigh << std::endl;
				} else {
					repl_neigh[which_neigh] = np;
					++repl_neigh_used;
				}
				break;
		}
	}

	void Cell::diff(){ // ONE Toffoli-Margoulus step
		static class rnarep::CellContent *temp_vals;

		temp_vals = vals; // 0 -> S
		if(gsl_rng_uniform(r) < 0.5) { //counterclockwise
			vals = diff_neigh[1]->vals; // 1 -> 0
			diff_neigh[1]->vals = diff_neigh[3]->vals; // 3 -> 1
			diff_neigh[3]->vals = diff_neigh[2]->vals; // 2 -> 3
			diff_neigh[2]->vals = temp_vals; // S -> 2

		} else { //cloclwise
			vals = diff_neigh[2]->vals; // 2 -> 0
			diff_neigh[2]->vals = diff_neigh[3]->vals; // 3 -> 2
			diff_neigh[3]->vals = diff_neigh[1]->vals; // 1 -> 3
			diff_neigh[1]->vals = temp_vals; // S -> 1
		}
	}

	double Cell::M(){
		double M = 1, akt = 0;

		for(int a = 0; a < par_noEA; a++){
			akt = 0;
			for(int met = 0; met < no_met_neigh; met++){
				// M(x) = prod(sum (a_i))
				akt += met_neigh[met]->vals->geta(a) ;
			}
			M *= akt;
		}

		return std::pow(M, reciproc_noEA);
	}

	void Cell::update(){
		double sum = par_claimEmpty;
		int decision = 0;

		if (vals->empty) { // if focal cell is empty (it can be filled with a copy)
			//REPLICATION
			for(int rep = 1; rep < no_repl_neigh; rep++) { // 0. neighbour is self, but it is empty
				if(!repl_neigh[rep]->vals->empty){
					sum += (claims[rep] = repl_neigh[rep]->vals->getR() * repl_neigh[rep]->M() );
				}
				else claims[rep] = 0.0;
			}
			//decision
			if(sum){
				decision = dvtools::brokenStickVals(claims, no_repl_neigh + 1, sum, gsl_rng_uniform(r)) ;
//				if(decision || (sum!=par_claimEmpty) ) {std::cout << "Replication: decision is: " << decision << " from claims:" << std::endl;
//									for(int op = 0; op < no_repl_neigh+1;op++ ) std::cout << claims[op]/sum << "\t";
//									std::cout << std::endl;}
				if(decision){ //claim 0 is claimEmpty NOTE that the probablity of staying empty is not fixed (e.g. 10%)! In case decision is negative see: brokenStickVals
						vals->replicate( *(repl_neigh[decision]->vals) ); //it is an empty cel, so no need to kill it first
						if(gsl_rng_uniform(r) < 0.5) { //havet to switch them at 50 percent
							switchit( *(repl_neigh[decision]) );
						}
						parent->cum_replications++;
//						no_births++;
//						std::cout << "Replication happend. The two molecules:" << std::endl << *(vals->get_seq()) << std::endl << *(repl_neigh[decision]->vals->get_seq()) << std::endl;  
				}
			}
		}
		else { //if focal cell is occupied (it can die)
			//DEGRADATION
			if(vals->Pdeg > gsl_rng_uniform(r) ) {
//				no_deaths++;
				parent->cum_deaths++;				
//				std::cout << "Degradation with Pdeg " << vals->Pdeg << std::endl;
				vals->die();
			}
		}
	}





	//FUNCTIONS FOR CellAut
	int CellAut::grid_init() {
		size = 0;
		
		//calculating sizes
		if(nrow <=0 || ncol <= 0) {
			size =0;
//			std::cout << "grid_init: size is 0" << std::endl;
			return(0);
		} 
		else if(nrow == 1 ) {
//			std::cout << "grid_init: nrow = 1" << std::endl;
			size = ncol;
		}
		else if(ncol == 1){
//			std::cout << "grid_init: ncol = 1" << std::endl;
			size = nrow;
		}
		else {
//			std::cout << "grid_init: matrix is more than 1D" << std::endl;
			switch (layout) { //the matrix is definitely more than 1D based on sizes
				case square:
				case hex:
					size = nrow * ncol;
					break;
				default:
					break;
			}
		}
//		std::cout << "grid_init: here2: ncol: " << ncol << ", nrow: " << nrow << ", size: " << size << std::endl;
		
		//initialising grid
		if(layout == square || layout == hex) {
			matrix = new Cell[size];
			if(! matrix ) {
				std::cerr << "ERROR: cadv::ca_init: matrix could not be initialized!" << std::endl;
				return(0);
			}
		}
//		std::cout << "grid_init: ncol: " << ncol << ", nrow: " << nrow << ", size " << size << std::endl;
		
		//setting parent
		for(int i=0; i < size; i++){
			matrix[i].parent = this;
		}

		return(size);
	}

	///gives back coordinates of nth cell
	std::vector<int> CellAut::getCoord(int n) {
		std::vector<int> coords;
		if (layout == square){
			coords.push_back(n % ncol); //x coord
			coords.push_back(n / ncol); //y coord
		}
		else if (layout == hex) {
			//x=n %/% ncol
			//y = n - x*ncol - x%/%2
			//z = 0-x-y
			coords.push_back( n / ncol); //x cubic coord
			coords.push_back( n - coords[0]*ncol - coords[0]/2 ); //y cubic coord
			coords.push_back(0 - coords[0] - coords[1]); //z cubic coord
		}

		return(coords);
	}

	///initialises matrix with predefined values, randomly
	void CellAut::init(std::string *pool, double* probs, int no_choices) {
		int i = 0;
		double sum = 0.0;
		
		for(i=0; i < no_choices; i++) {
			sum += probs[i];
		}
		for(i=0; i < size; i++) {
			*(matrix[i].vals) = pool[dvtools::brokenStickVals(probs, no_choices, sum, gsl_rng_uniform(r) )];
		}
	}	

	void CellAut::init_fromfile(char *infile) {
		std::string line, word;

		std::ifstream file(infile);
		Cell *cell = matrix;

		if(!file.is_open()) std::cerr << "ERROR: init_fromfile: file can not be opened!" << std::endl;

		rnarep::CellContent::no_replicators = 0; //clearing number of replicators 

		for (int cellnum = 0; cellnum < size && std::getline(file, line); cellnum++){
		//while (std::getline(file, line) ){
			std::istringstream linestream(line);
			linestream >> word;
//			std::cout << "init_fromfile assigning " << word << std::endl;			
			*(cell->vals) = word;
			cell++;
		}

		if(cell != (Cell *) matrix + size) std::cerr << "WARNING: file length is not equal to gridsize!" << std::endl;
//		std::cout << "Grid initialised with " << rnarep::CellContent::no_replicators << " replicators on a grid of " << size << " cells." << std::endl;
	}

	inline Cell* CellAut::get(int x, int y) {
		if(layout == square){
			return(matrix + y*ncol + x);
		}
		return(NULL);
	}
	///finds cell in pos [x,y] and gives back its number
	inline int CellAut::getN(int x, int y) {
		if(layout == square){
			return(y*ncol + x);
		}
		else if (layout == hex){
			return(y + x*ncol + x/2);
		}
		return(-1);
	}

	//a singel update step <- this is called by rUpdate and oUpdate
	int CellAut::updateStep(int cell){
		matrix[cell].update();
		return(0);
	}

	///Random update
	int CellAut::rUpdate(int gens){
		//setting up signal hadler
		signal(SIGINT, CellAut::static_signalHandler);
		//std::cerr << "Signal listener started" << std::endl;

		int iter=0, diff_until = dvtools::fracpart(diff * size * time);

		//check if output is open
		if(par_output_interval && !output) {
			std::cerr << "ERROR: output not open" << std::endl;
			return(-1);
		}

		//output/save in case of not init from start
		if(time){
			if(par_output_interval && (time % par_output_interval)) do_output();
			if(par_save_interval && (time % par_save_interval)) save();
		}

		for(mtime = time + gens ; time < mtime && rnarep::CellContent::no_replicators > 0; time++){ //updating generations
			// bubble extinction event
			if (par_bubble_interval && !(time % par_bubble_interval)) for(int repeta = par_no_bubi; repeta--;) bubble_sampling( std::abs(gsl_ran_gaussian(r, par_sd_bubblesize) + par_mean_bubblesize ) );

			//outputs
			if (par_output_interval && !(time % par_output_interval)) do_output();
			if (par_save_interval && !(time % par_save_interval)) save();

			for(iter = 0; iter < size; iter++){
				//UPDATING
				matrix[ gsl_rng_uniform_int(r, size) ].update();
				//DIFFUSION
				for(diff_until += diff; diff_until >= 1; diff_until--){
					matrix[ gsl_rng_uniform_int(r, size) ].diff();
				}
				// Secondry output
				if(cum_replications > (par_output_interval*size) ){
					do_output(time + static_cast<double>(iter+1)/static_cast<double>(size));
					cum_deaths = cum_replications = 0;
				}
			}
//			std::cout << "Cycle " << time << ": number of total deaths: " << no_deaths << ", number of total births: " << no_births << std::endl;
		}
		
		// saving/outputting
		if(par_output_interval) do_output();
		if(par_save_interval) if(save()) return -2;


		return rnarep::CellContent::no_replicators;
	}

	//Update according to a random order (in every generation all cells will be updated)
	int CellAut::oUpdate(int gens){
		int *order;
		int iter=0, temp = 0, target = 0, diff_until = dvtools::fracpart(diff * time);

		//check if output is open
		if(par_output_interval && !output) {
			std::cerr << "ERROR: output not open" << std::endl;
			return(-1);
		}

		//init order
		order = new int[size];
		for(iter = 0; iter < size; iter++){
			order[iter] = iter;
		}

		for(mtime = time + gens ; rnarep::CellContent::no_replicators > 0 && time < mtime ; time++){ //updating generations
			//outputs
			if (par_output_interval && !(time % par_output_interval)) do_output();
			if (par_save_interval && !(time % par_save_interval)) save();

			//UPDATING
			for(iter = 0; iter < size; iter++){
				target = gsl_rng_uniform_int(r, size - iter);
				if (target) {
					target += iter;
					temp = order[target];
					order[ target ] = order[ iter ];
					order[ iter ] = temp;
					//updateStep( temp );
					matrix[ temp ].update();
				}
				else {
					//updateStep( order[iter] );
					matrix[ order[iter] ].update();
				}
			}
			//DIFFUSION
			for(diff_until += diff; diff_until >= 1; diff_until--){
				for(iter = 0; iter < size; iter++){
					target = gsl_rng_uniform_int(r, size - iter);
					if (target) {
						target += iter;
						temp = order[target];
						order[ target ] = order[ iter ];
						order[ iter ] = temp;
						matrix[ temp ].diff();
					}
					else {
						matrix[ order[iter] ].diff();
					}
				}
			}
		}

		if(par_save_interval) save();
		if(std::strlen(par_output_filename)) do_output();

		delete [] (order);

		return(0);
	}

	//gives back coordinates of nth cell
	std::vector<int> CellAut::getCoord(int n, int type = 0) {
		std::vector<int> coords;
		if (layout == square){
			coords.push_back(n % ncol); //x coord
			coords.push_back(n / ncol); //y coord
		}
		else if (layout == hex) {
			if(type == 0){ //cube coordinates - axial coordinates are the choosen two from this three coords
				coords.push_back( n / ncol); //x cubic coord
				coords.push_back( n - coords[0]*ncol - coords[0]/2 ); //y cubic coord
				coords.push_back(0 - coords[0] - coords[1]); //z cubic coord
			}
			else if (type == 1){ //axial coords 
				coords.push_back( n / ncol); //x cubic coord
				coords.push_back( n - coords[0]*ncol - coords[0]/2 ); //y cubic coord
			}
			else if (type == 2){ //offset coords 
				coords.push_back(n % ncol); //x coord
				coords.push_back(n / ncol); //y coord
			}
		}

		return(coords);
	}

	void CellAut::blueprintNeigh(const double neigh_tipus, std::vector<int> & n_inic_x, std::vector<int> & n_inic_y){
		int maxDist = 0, x = 0, y = 0;

		n_inic_x.clear();
		n_inic_y.clear();

		//Create neighbourhood definition
		if(neigh_tipus == MARGOLUS_NEIGH) {
		    n_inic_x.push_back(0);
		    n_inic_y.push_back(0);

		    n_inic_x.push_back(1);
		    n_inic_y.push_back(0);

		    n_inic_x.push_back(0);
		    n_inic_y.push_back(1);

		    n_inic_x.push_back(1);
		    n_inic_y.push_back(1);
		}

		else{
		    if(CellAut::layout == square) {
			    if(2 <= neigh_tipus ) {
			    	//self
				n_inic_x.push_back(0);
				n_inic_y.push_back(0);
			    
				//other cells
				//maxDist = (int) std::log2((int) neigh_tipus - 1);
				maxDist = (int) std::sqrt(neigh_tipus)+1;
				for(x = -maxDist; x <= maxDist; x++) for(y = -maxDist; y <= maxDist; y++){ 
//						std::cout << "std::pow(2, " << x << ") + std::pow(2," << y << ") <= " << neigh_tipus << "\t" << std::pow(2, std::abs(x)) << " " << std::pow(2, std::abs(y)) << std::endl; 
						if( (x || y) && (x*x + y*y <= neigh_tipus ) ) {
							n_inic_x.push_back(x);
							n_inic_y.push_back(y);
						}
				} //end for x and y
			    }
		    }
		    else if(layout == hex){
			    if(3 <= neigh_tipus ) {
			    	//self
				n_inic_x.push_back(0);
				n_inic_y.push_back(0);

				//other cells
				maxDist = (int) std::log2((int) neigh_tipus - 2);
				for(x=-maxDist; x <= maxDist; x++) for(y=-maxDist; y <= maxDist; y++){ 
						if( (x || y) && (std::pow(2, std::abs(x)) + std::pow(2, std::abs(y)) + std::pow(2, std::abs(0-x-y)) <= neigh_tipus ) ) {
							n_inic_x.push_back(x);
							n_inic_y.push_back(y);
						}
				} //end for x and y
			    }
			    
		    }
		}
	}

	int CellAut::calcNeigh(int i, unsigned int n, std::vector<int> n_inic_x, std::vector<int> n_inic_y, const Ca_edge edge) const{
		if(edge == torus){
			if(layout == square){
					  return ( ( dvtools::Rmod( ((int)i/ncol + n_inic_y[n]) , nrow) ) * ncol + dvtools::Rmod( dvtools::Rmod(i , ncol) + n_inic_x[n] , ncol));
			}
			else if(layout == hex){
					   return(  dvtools::Rmod((int)i/ncol + n_inic_y[n] , nrow )  * ncol + dvtools::Rmod ( dvtools::Rmod( i, ncol) + n_inic_x[n] + ( n_inic_y[n] + (( (int)i / ncol)&1)  )/2 , ncol)); 
			}
		}
		else if (edge == wall){
				const int x = (int) i / ncol, y = i % ncol;

				if(x + n_inic_x[n] >= 0 && x + n_inic_x[n] < ncol && y + n_inic_y[n] >= 0 && y + n_inic_y[n] < nrow  ) {
					return (i + n_inic_x[n] * ncol + n_inic_y[n]); 	
				}  
		}
		else if (edge == mirror){ //does not work!!!
		}

		return -1;
	}

	unsigned int CellAut::countNoNeigh(const Ca_edge edge, std::vector<int> & n_inic_x, std::vector<int> & n_inic_y, const int i) const{
		if(edge == wall){
			const int x = (int) i / ncol, y = i % ncol;
			unsigned int noNei = 0;
			for(unsigned int n = 0; n < n_inic_x.size(); n++) {
				if(x + n_inic_x[n] >= 0 && x + n_inic_x[n] < ncol && y + n_inic_y[n] >= 0 && y + n_inic_y[n] < nrow  ) noNei++;
			}
			return noNei;
		} else {
		        return n_inic_x.size();
		}
		//return 0;
	}

	void CellAut::neighInic(const double neigh_tipus, const Ca_edge edge, int neigh_no = 1) {
		std::vector<int> n_inic_x;
		std::vector<int> n_inic_y;

		blueprintNeigh(neigh_tipus, n_inic_x, n_inic_y);

		for(int i=0; i < size; i++){ //iterate throught grid
			// count number of neighbours	
			unsigned int noNei = countNoNeigh(edge, n_inic_x, n_inic_y, i);

			// allocate space for neighbou pointers
			matrix[i].inicNeigh(noNei, neigh_no);
			
			//assign neighbours
			for(unsigned int n = 0; n < noNei; ++n) {	
				int myneigh = calcNeigh(i, n, n_inic_x, n_inic_y, edge);
				if( myneigh > -1 ) {
					matrix[i].setNeigh( matrix + myneigh, neigh_no); 	
				}
			}
		}
	} //end neighInic

	
	int CellAut::openOutputs(){
		std::string name, command;
		 
		name += par_outdir;
		
		//create directory output if it does not exist (Linux only!!)
		command += "mkdir -p ";
		command += par_outdir;
//		std::cout << command << std::endl;
		system(command.c_str());

		//create directory for output
		name += "/";
		name += par_ID;

		command.clear();
		command += "test -d ";
		command += name;
		
//		std::cout << command << std::endl;

		if(system(command.c_str())) { //directory does not exist
//			std::cout << "not exists" << std::endl;

			command.clear();
			command += "mkdir ";
			command += name;
			
//			std::cout << command << std::endl;
			system(command.c_str());
		}
		else{ //directory already exist
//			std::cerr << "exists" << std::endl;
			//check if files already exists
			command.clear();
			command += "test -f ";
			command += name;
			command += "/";
			command += par_output_filename;
//			std::cout << command << std::endl;
			if(!system(command.c_str())){ //exists already
				//try a new directory name
				name += "_";
				name += std::to_string(gsl_rng_uniform_int(r, 1000));
				command.clear();
				command += "test -d ";
				command += name;
//				std::cout << command << std::endl;
				if(!system(command.c_str())) { //already exzist
					std::cerr << "ERROR: conflicting simulations with identical IDs. Quitting..." << std::endl;
					return (1);
				}
				else{ // new dir does not exist
					//std::cerr << "WARNING: directory already existed with output files! Tried a new name: " << name << std::endl;
					command.clear();
					command += "mkdir ";
					command += name;
					
//					std::cout << command << std::endl;
					system(command.c_str());
				}
			} //files exist
		} //directory exists

		name += "/";
		//output directory exists and its path is in name

		//creating output file
		command.clear();
		command += name;
		command += par_output_filename;

		output.open(command);
		if(!output.is_open()) {
			std::cerr << "ERROR: output file (" << name << par_output_filename << ") cant be opened" << std::endl;
			return (2);
		}

		//adding header to output
		output << "time;replicators";
		output << ";no_par;mean_R_par;mean_length_par;mean_mfe_par" ;
		output << ";no_templ;mean_R_templ;mean_length_templ;mean_mfe_templ" ;
		for(int e = 0; e < par_noEA; e++){
			output << ";no_enz" << e << ";mean_R_enz" << e << ";mean_length_enz" << e << ";mean_mfe_enz" << e << ";mean_a_enz" << e ;
		}

		for(int e = 0; e < par_noEA; e++){
			output << ";no_Genz" << e << ";mean_R_Genz" << e << ";mean_length_Genz" << e << ";mean_mfe_Genz" << e << ";mean_a_Genz" << e ;
		}

		for(int a = 0; a <= par_noEA; a++) {
			output << ";no_A" << a ;
		}

		output << ";no_replications;no_deaths";

		output << std::endl;

		output.flush();

		//prepare vectors for output
		
		out_gen.reserve(par_noEA);
		out_spec.reserve(par_noEA);
		out_noA.reserve(par_noEA + 1);

		//creating SAVE directory
		command.clear();
		command += "mkdir -p ";
		
		name += par_savedir;

		command += name;
//		std::cout << command << std::endl;
		system(command.c_str());

		savedir = name;

		return 0;

	}

	void CellAut::bubble_sampling(const double bubblesize){
		// get blueprint for neigbourhood
		std::vector<int> x, y;
		blueprintNeigh(bubblesize, x, y);

		// count number of neighbours	
		unsigned int no_neigh = countNoNeigh(torus, x, y);

		if(no_neigh){
			// get filename for output
			std::string filename, command;
			unsigned int no_bubble=0;

			if(!savedir.length()) {
				std::cerr << "ERROR: No savedir inicialised! Please do run CellAut::openOutputs() before saving!" << std::endl;
				return;
			}
			
			filename = savedir; 
			filename += "/bubble_t"; 
			filename += std::to_string(time);
			filename += '_';
			
			//check if file already exists
			command.clear();
			command += "test -f ";
			command += filename;
			command += std::to_string(no_bubble);
			command += ".tsv";
			while(!system(command.c_str())){ //exists already
				command.clear();
				command += "test -f ";
				command += filename;
				command += std::to_string(++no_bubble);
				command += ".tsv";
			}
	
			filename += std::to_string(no_bubble);
			filename += ".tsv";

			// open output filestream
			std::ofstream out(filename);
			if(!out){
				std::cerr << "ERROR: Could not open file " << filename << " for saving bubble!" << std::endl;
				return;
			}

			// get random position
			const unsigned int middle = gsl_rng_uniform_int(r, size); 

			// get neighbours
			for(unsigned int n = 0; n < no_neigh; ++n){
				int neigh = calcNeigh(middle, n, x, y, torus);
				if( neigh > -1 ) {
					auto cellcont = matrix[neigh].vals;
					if(!cellcont->empty){
						// write it out
						out 	<< *(cellcont->get_seq())
							<< '\t' << cellcont->get_str()
							<< '\t' << cellcont->get_mfe() 
							<< '\t' << cellcont->getPfold()
							<< '\t' << cellcont->Pdeg 
							<< '\t' << cellcont->get_no_sites()
							<< '\t' << cellcont->getR()
							<< '\t' << 0 //cell->M()
							<< '\t' << cellcont->get_type() ; 
						for(double *a = cellcont->geta(), *a_until = cellcont->geta() + par_noEA; a != a_until; a++){
							out << '\t' << *a;
						}

						out 	<< '\t' << cellcont->get_prev_type()
							//<< '\t' << cellcont->get_type_rev() 
							<< std::endl;

						// kill it
						cellcont->die();
					}
				}
			}
			
			//close output filestream
			out.close(); // unneccessary
		}

	}

	void CellAut::do_output(){
		do_output(time);
	}

	void CellAut::do_output(const double otime){
		/* what i need:
			time, alive, [by akt: No, Rs mean, length mean, alpha mean, mfe mean], [by no akts: number]

		*/

		//clearing
		out_spec.assign(par_noEA, Outdata());
		out_gen.assign(par_noEA, Outdata());
		out_noA.assign(par_noEA + 1, 0);
		out_par = Outdata();
		out_templ = Outdata();
					      
		//calculating values
		for(Cell *cell = matrix, *end = (Cell *) matrix + size ; cell != end; cell++){
			if(!cell->vals->empty){ // if cell is not empty
				//how much activities does it have?
				int no_acts = cell->vals->get_no_acts();

				out_noA[no_acts]++;

				if(no_acts){ //it is not a parasite
					if(no_acts == 1){ // is a specialist
						// find which one is it (zero indexed: type==1 -> 0)
						unsigned int curract = 0;
						{
							auto n = cell->vals->get_type();
							while(n >>= 1) ++curract;
						}

						// store values
						out_spec[curract].add(cell->vals->getR(), cell->vals->get_length(), cell->vals->geta(curract), cell->vals->get_mfe());
					} else { // is a generalist
						for(int ea = 0; ea < par_noEA; ea++){
							double activity = cell->vals->geta(ea); //indexing of "out_" arrays starts at parazite, "a" starts with activity 0
							if(activity) { //if it has activity ea
								out_gen[ea].add(cell->vals->getR(), cell->vals->get_length(), activity, cell->vals->get_mfe());
							}
						}
					}
				}
				else { //it has no activity
					if(cell->vals->get_prev_type()){ //it is a template
						out_templ.add(cell->vals->getR(), cell->vals->get_length(), 0.0, cell->vals->get_mfe());
					} else { //it is a parazite
						out_par.add(cell->vals->getR(), cell->vals->get_length(), 0.0, cell->vals->get_mfe());
					}
				} 
			} // cell not empty
		} // tru cells in matrix

		//outputting
		// general data
		output << otime << ';' << rnarep::CellContent::no_replicators;

		//parasites
		if(out_par.no){
			output << ';' << out_par.no << ';' << out_par.avg_R() << ';' << out_par.avg_length() << ';' << out_par.avg_mfe();
		} else {
				output << ";0;0;0;0";
		}

		//templates
		if(out_templ.no){
			output << ';' << out_templ.no << ';' << out_templ.avg_R() << ';' << out_templ.avg_length() << ';' << out_templ.avg_mfe();
		} else {
				output << ";0;0;0;0";
		}

		//specialists
		unsigned long long int no;
		for(int ea = 0; ea < par_noEA; ea++) {
			Outdata *rep = &out_spec[ea]; 
			if((no = rep->no)){
				output << ';' << no << ';' << rep->avg_R() << ';' << rep->avg_length() << ';' << rep->avg_mfe() << ';' << rep->avg_a();
			}
			else {
				output << ";0;0;0;0;0";
			}

		}

		// generalists
		for(int ea = 0; ea < par_noEA; ea++) {
			Outdata *rep = &out_gen[ea]; 
			if((no = rep->no)){
				output << ';' << no << ';' << rep->avg_R() << ';' << rep->avg_length() << ';' << rep->avg_mfe() << ';' << rep->avg_a();
			}
			else {
				output << ";0;0;0;0;0";
			}

		}

		//by A
		for(int ea = 0; ea <= par_noEA; ea++) {
			output << ';' << out_noA[ea] ;
		}


		// cummulative data
		output << ';' << cum_replications << ';' << cum_deaths;

		output << std::endl;

		output.flush();

	}


	int CellAut::save(){
		/* Outputs:
		- text file, tab separated, each line represents a grid point from 0th to last
			values:
			seq str mfe Pfold Pdeg no_sites R M type [activities] prev_type
		- rng binary state file
		*/
		std::string emptystring("N\tN\t0\t-1\t-1\t-1\t-1\t-1\t0");
		std::string filename;

		//prepare string for empty cells
		for(int ea = 0; ea < par_noEA; ea++) emptystring += "\t0";
		emptystring += "\t-1\n";

		//open outputs
		if(!savedir.length()) {
			std::cerr << "ERROR: No savedir inicialised! Please do run CellAut::openOutputs() before saving!" << std::endl;
			return (1);
		}
		
		filename = savedir; 
		filename += '/'; 
		filename += std::to_string(time);
		filename += ".tsv";
		std::ofstream out(filename);
		
		if(!out){
			std::cerr << "ERROR: Could not open file for saving grid data!" << std::endl;
			return 2;
		}

		//going throught grid
		for(Cell *cell = matrix, *end = (Cell *) matrix + size ; cell != end; cell++){
			rnarep::CellContent *cellcont = cell->vals;
			//output values:
			// seq str mfe Pfold Pdeg no_sites R type [alphas]
			if(cellcont->empty) out << emptystring; 
			else {
				out 	<< *(cellcont->get_seq())
					<< '\t' << cellcont->get_str()
					<< '\t' << cellcont->get_mfe() 
					<< '\t' << cellcont->getPfold()
					<< '\t' << cellcont->Pdeg 
					<< '\t' << cellcont->get_no_sites()
					<< '\t' << cellcont->getR()
					<< '\t' << cell->M()
					<< '\t' << cellcont->get_type() ; 
				for(double *a = cellcont->geta(), *a_until = cellcont->geta() + par_noEA; a != a_until; a++){
					out << '\t' << *a;
				}

				out 	<< '\t' << cellcont->get_prev_type()
				 	//<< '\t' << cellcont->get_type_rev() 
					<< std::endl;
			}

			//out.flush();
		}

		out.close();

		//saving rng state
		filename.clear();
		filename = savedir;
		filename += "/rngstate";
		filename += std::to_string(time);
		filename += ".bin";
		//std::ofstream rngout(filename, std::ios::out | std::ios::binary);

		if(randomszam_mentes(filename.c_str(), r)) {
			std::cerr << "ERROR: Could not open file for saving random number generator state!" << std::endl;
			return 3;
		}

		return 0;
	}

	//signal handler
	void CellAut::signalHandler(int signal){
		mtime = time;
	}


}

