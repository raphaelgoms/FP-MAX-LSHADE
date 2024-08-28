/*
  L-SHADE implemented by C++ for Special Session & Competition on Real-Parameter Single Objective Optimization at CEC-2014
  See the details of L-SHADE in the following paper:

  * Ryoji Tanabe and Alex Fukunaga: Improving the Search Performance of SHADE Using Linear Population Size Reduction,  Proc. IEEE Congress on Evolutionary Computation (CEC-2014), Beijing, July, 2014.

  Version: 1.0   Date: 16/Apr/2014
  Written by Ryoji Tanabe (rt.ryoji.tanabe [at] gmail.com)
*/

#include"de.h"
// #include "cec14_test_func.cpp"
//#include "cec22_test_func.cpp"

void Tempo_CPU_Sistema(double *seg_CPU_total, double *seg_sistema_total)
{
  long seg_CPU, seg_sistema, mseg_CPU, mseg_sistema;
  struct rusage ptempo;

  getrusage(0,&ptempo);

  seg_CPU = ptempo.ru_utime.tv_sec;
  mseg_CPU = ptempo.ru_utime.tv_usec;
  seg_sistema = ptempo.ru_stime.tv_sec;
  mseg_sistema = ptempo.ru_stime.tv_usec;

 *seg_CPU_total     = (seg_CPU + 0.000001 * mseg_CPU);
 *seg_sistema_total = (seg_sistema + 0.000001 * mseg_sistema);
}

void cec14_test_func(double *,double *,int, int,int);
//void cec22_test_func(double *, double *,int,int,int);

int isDirectory(const char *path) {
   struct stat statbuf;
   if (stat(path, &statbuf) != 0)
       return 0;
   return S_ISDIR(statbuf.st_mode);
}

double difference(Individual i1, Individual i2) {
  double diff = 0;
  for (size_t i = 0; i < g_problem_size; i++) {
      diff += (i1[i] - i2[i]) * (i1[i] - i2[i]);   
  }
  return sqrt(diff);
}


double computeDiversity(vector<Individual> population) {
  double diff_acc = 0;
  for (Individual ind1 : population) {
    for (Individual ind2 : population) {
      diff_acc += difference(ind1, ind2);
    }
  }
  return diff_acc / population.size(); 
}


void fprintPopulation(vector <Individual> pop, int generation, string basepath) {
  
  string fullpath = basepath + "/population/g" + to_string(generation) + ".csv";
  std::ofstream pop_file;

  pop_file.open(fullpath);

  for (int j = 0; j < g_problem_size; j++) {
    pop_file << "x" << j+1;

    if (j < g_problem_size - 1)
      pop_file << ";";
  }
  pop_file << endl;

  // TODO: persist pop in file
  for (int i = 0; i < pop.size(); i++)
  { 
    for (int j = 0; j < g_problem_size; j++)
    {
      pop_file << pop[i][j];

      if (j < g_problem_size - 1)
        pop_file << ";";
    }
    pop_file << endl;
  }

  pop_file.close();
}

void fprintElite(vector<tuple<Individual, double>> & elite, int generation, string basepath) {
  string fullpath = basepath + "/elite/g" + to_string(generation) + ".csv";
  std::ofstream patterns_file;

  patterns_file.open(fullpath, std::ofstream::out);
  for (int j = 0; j < g_problem_size; j++) {
    patterns_file << "x" << j+1;

    if (j < g_problem_size - 1)
      patterns_file << ";";
  }
  patterns_file << endl;

  for (int i = 0; i < elite.size(); i++)
  { 
    for (int j = 0; j < g_problem_size; j++)
    {
      patterns_file << get<0>(elite[i])[j];

      if (j < g_problem_size - 1)
        patterns_file << ";";
    }
    patterns_file << endl;
  }

  patterns_file.close();
}



searchAlgorithm::searchAlgorithm() {

}

void searchAlgorithm::initializeParameters() {
  //function_number = g_function_number;
  problem_size = g_problem_size;
  max_num_evaluations = g_max_num_evaluations;
  pop_size = g_pop_size;
  initializeFitnessFunctionParameters();
}

void searchAlgorithm::evaluatePopulation(const vector<Individual> &pop, vector<Fitness> &fitness) {

  for (int i = 0; i < pop_size; i++) {
    if (function_number > -1) {
      fitness[i] = evaluateIndividual(pop[i]);
      //cout << function_number << fitness[i] << " ";
    }  
  }
}

Fitness searchAlgorithm::evaluateIndividual(Individual individual) {
    
    Fitness fitness;

    // if (coco_exp_func) {
    //   coco_exp_func(individual, fitness);
    //   return fitness[0];
    // }
    cec14_test_func(individual, &fitness, problem_size, 1, function_number);
    return fitness;
}

void searchAlgorithm::initializeFitnessFunctionParameters() {
  //epsilon is an acceptable error value.
  epsilon = pow(10.0, -8);
  max_region = g_max_region;
  min_region = g_min_region;

  lower_bounds = vector<variable>(problem_size);
  upper_bounds = vector<variable>(problem_size); 

  for (int i = 0; i < problem_size; i++)
  {
    lower_bounds[i] = min_region;
    upper_bounds[i] = max_region;
  }

  //optimum = 0; // CEC2013
  optimum = g_function_number * 100; // -> 2014
  // optimum =  g_optimum[g_function_number - 1]; // -> 2022
  // optimum = f->getMinValue();
}

//set best solution (bsf_solution) and its fitness value (bsf_fitness) in the initial population
void searchAlgorithm::setBestSolution(const vector<Individual> &pop, const vector<Fitness> &fitness, Individual &bsf_solution, Fitness &bsf_fitness) {
  int current_best_individual = 0;

  for (int i = 1; i < pop_size; i++) {
    if (fitness[current_best_individual] > fitness[i]) {
      current_best_individual = i;
    }
  }

  bsf_fitness = fitness[current_best_individual];
  for (int i = 0; i < problem_size; i++) {
    bsf_solution[i] = pop[current_best_individual][i];
  }
}

// make new individual randomly
Individual searchAlgorithm::makeNewIndividual() {
  Individual individual = (variable*)malloc(sizeof(variable) * problem_size);

  for (int i = 0; i < problem_size; i++) {
    individual[i] = ((upper_bounds[i] - lower_bounds[i]) * randDouble()) + lower_bounds[i];
  }

  return individual;
}

/*
  For each dimension j, if the mutant vector element v_j is outside the boundaries [x_min , x_max], we applied this bound handling method
  If you'd like to know that precisely, please read:
  J. Zhang and A. C. Sanderson, "JADE: Adaptive differential evolution with optional external archive,"
  IEEE Tran. Evol. Comput., vol. 13, no. 5, pp. 945–958, 2009.
 */
void searchAlgorithm::modifySolutionWithParentMedium(Individual child, Individual parent) {
  int l_problem_size = problem_size;
  variable l_min_region = min_region;
  variable l_max_region = max_region;

  for (int j = 0; j < l_problem_size; j++) {
    if (child[j] < l_min_region) {
      child[j]= (l_min_region + parent[j]) / 2.0;
    }
    else if (child[j] > l_max_region) {
      child[j]= (l_max_region + parent[j]) / 2.0;
    }
  }
}

tuple<double, double> searchAlgorithm::lineSearch(Individual current_position, double h, int axis, bool stop_in_first_improve) {
    int j, start, end, z;
    double ax_pos, f_curr, s_domain, e_domain;

    Individual x_curr = current_position;
    
    ax_pos = x_curr[axis];
    s_domain = lower_bounds[axis];
    e_domain = upper_bounds[axis];

    double f_star = evaluateIndividual(x_curr);
    double pos = lower_bounds[axis], best_pos = ax_pos;

    start = (int)floor((ax_pos - s_domain) / h);
    double step;

    int i=0;
    for (j = 0; j < start; j++)
    { 
        step = h * pow(2, i); i++;
        x_curr[axis] = ax_pos - step;
        if (x_curr[axis] < s_domain)
          break;
        
        f_curr = evaluateIndividual(x_curr);

        if (f_curr < f_star)
        {
            f_star = f_curr;
            best_pos = x_curr[axis];

            if (stop_in_first_improve) 
              return { best_pos, f_star };
        }
    }

    end = (int)floor((e_domain - ax_pos) / h);
    i = 0;
    for (j = 0; j <= end; j++)
    {
        step = h * pow(2, i); i++;
        x_curr[axis] = ax_pos + step;
        if (x_curr[axis] > e_domain)
          break;
        f_curr = evaluateIndividual(x_curr);

        if (f_curr < f_star)
        {
            f_star = f_curr;
            best_pos = x_curr[axis];

            if (stop_in_first_improve) 
              return { best_pos, f_star };;
        }
    }

    return { best_pos, f_star };
}

/**
 * TODO : contruct a solution from the hibrid of patterns and p-best(momentum),
 *  in a greed-randomized way. Stop method when find the fist cost improve.
 * 
 * 
 * H1: If I start with random values in unfixed positions, find improve will
 *  probably be easy, but the solution quality low. 
 *  -> Possibly this will improve the diversity
 *  
 * 
 * 
 * To remember: 
 * - First improve, can has two levels:
 *    1. Solution level
 *    2. Axis level
*/
bool searchAlgorithm::constructGreedyRandomized(Individual& current_position, set<int> unfixed_positions, double h, double alpha)
{
    if (!alpha) 
      alpha = randDouble();

    set<int> rcl; // Restricted canditade list
    bool reuse = false, imprc = false;

    double gplus, gminus;
    Individual z = (variable*)malloc(sizeof(variable) * problem_size);
    Individual _x = (variable*)malloc(sizeof(variable) * problem_size);
    vector<float> g(problem_size);
    
    double curr_cost = evaluateIndividual(current_position);
    //cout << "=========================>>>> constructGreedyRandomized <<<<==============================="<< endl;
    while (unfixed_positions.size()) {
        gminus = INFINITY;
        gplus = -INFINITY;

        for (auto& i : unfixed_positions) {
            if (!reuse) {
                
                tuple<double, double> result = lineSearch(current_position, h, i, true);
                // cout << " z_i = "<< z[i] << endl; 
                z[i] = get<0>(result);
                g[i] = get<1>(result);

                // _x = current_position; _x[i] = z[i];
                // g[i] = evaluateIndividual(_x);
                if (g[i] < curr_cost)
                  return true;
                // cout << " g_i = "<< g[i] << endl; 
            }
            if (gminus > g[i]) gminus = g[i];
            if (gplus < g[i]) gplus = g[i];
        }

        rcl.clear();
        double threshold = gminus + alpha * (gplus - gminus);
        for (auto& i : unfixed_positions) {
            if (g[i] <= threshold) {
                rcl.insert(i);
            }
        }

        int j = radomlySelectElement(rcl);
        if (current_position[j] == z[j]) {
            reuse = true;
        }
        else {
            current_position[j] = z[j];
            reuse = false;
            imprc = true;
            return imprc;
        }
        unfixed_positions.erase(j); //j
        //cout << "i = " << unfixed_positions.size() << endl;
        // return false;
    }

    return imprc;
}