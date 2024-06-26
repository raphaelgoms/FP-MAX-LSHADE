/*
  L-SHADE implemented by C++ for Special Session & Competition on Real-Parameter Single Objective Optimization at CEC-2014

  Version: 1.0   Date: 16/Apr/2014
  Written by Ryoji Tanabe (rt.ryoji.tanabe [at] gmail.com)
*/

#include "de.h"
#include <iostream>

#include <boost/math/statistics/linear_regression.hpp>
#include <cstring>
#include <mlpack/methods/kmeans/kmeans.hpp>
#include <mlpack/core/cv/metrics/silhouette_score.hpp>
#include <mlpack/methods/kmeans/allow_empty_clusters.hpp>

using namespace mlpack::kmeans;
using namespace mlpack::cv;

using namespace std;
using boost::math::statistics::simple_ordinary_least_squares_with_R_squared;

using namespace std;

KMEANS_LSHADE::KMEANS_LSHADE(int cec_function_number, PatternSelectionStrategy pattern_sel_strategy, SolutionFillingStrategy filling_strategy)
{
  // dataset = new Dataset;

  this->function_number = cec_function_number;
  this->pattern_sel_strategy = pattern_sel_strategy;
  this->filling_strategy = filling_strategy;

  //m_generator.seed(static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count()));
}

Fitness KMEANS_LSHADE::run()
{   
  m_generator.seed(this->seed);
  this->function_calls = 0;
  this->used_gen_count = 0;

  // cout << scientific << setprecision(8);
  initializeParameters();
  setSHADEParameters();

  // cout << pop_size << endl;
  // cout << arc_size << endl;
  // cout << p_best_rate << endl;
  // cout << memory_size << endl;

  bool print = false;
  int generation = 1;

  clusters_count_min = INFINITY;
  clusters_count_max = 0;
  clusters_count_avg = 0;

  vector<Individual> pop;
  vector<Fitness> fitness(pop_size, 0);

  vector<Individual> children;
  vector<Fitness> children_fitness(pop_size, 0);

  // initialize population
  for (int i = 0; i < pop_size; i++)
  {
    pop.push_back(makeNewIndividual());
    children.push_back((variable *)malloc(sizeof(variable) * problem_size));
  }

  // evaluate the initial population's fitness values
  evaluatePopulation(pop, fitness);

  int initial_pop_size = pop_size;
  std::string results_root_path = "";
  std::ofstream pop_costs_file;

  if (debug_mode) {
    while (!isDirectory((results_root_path + "results").c_str()))
      results_root_path += "../";
    results_root_path += "results/generations/KMEANS_LSHADE/";

    string fullpath = results_root_path + "pop-costs.csv";

    int check = mkdir("results/generations",0777);
    check = mkdir("results/generations/KMEANS_LSHADE",0777);
    
   
    pop_costs_file.open(fullpath, std::ofstream::out | std::ofstream::app);

    for (int j = 0; j < pop_size; j++) {
      pop_costs_file << "i" << j+1;

      if (j < pop_size - 1)
        pop_costs_file << ";";
    }
    pop_costs_file << endl;
    
    fprintPopulation(pop, generation, results_root_path + "/generations/KMEANS_LSHADE/");

    for (int j = 0; j < pop_size; j++) {
      pop_costs_file << fitness[j]-optimum;

      if (j < pop_size - 1)
        pop_costs_file << ";";
    }
    pop_costs_file << endl;
  }
  

  Individual bsf_solution = (variable *)malloc(sizeof(variable) * problem_size);
  Fitness bsf_fitness;
  int nfes = 0;

  if ((fitness[0] - optimum) < epsilon)
    fitness[0] = optimum;
  bsf_fitness = fitness[0];
  for (int j = 0; j < problem_size; j++)
    bsf_solution[j] = pop[0][j];
  /////////////////////////////////////////////////////////////////////////
  for (int i = 0; i < pop_size; i++)
  {
    nfes++;

    if ((fitness[i] - optimum) < epsilon)
      fitness[i] = optimum;

    if (fitness[i] < bsf_fitness)
    {
      bsf_fitness = fitness[i];
      for (int j = 0; j < problem_size; j++)
        bsf_solution[j] = pop[i][j];
    }

    // if (nfes % 1000 == 0) {
    //   //      cout << nfes << " " << bsf_fitness - optimum << endl;
    //   cout << bsf_fitness - optimum << endl;
    // }

    if (nfes >= max_num_evaluations)
      break;
  }
  ////////////////////////////////////////////////////////////////////////////

  // for external archive
  int arc_ind_count = 0;
  int random_selected_arc_ind;
  vector<Individual> archive;
  for (int i = 0; i < arc_size; i++)
    archive.push_back((variable *)malloc(sizeof(variable) * problem_size));

  int num_success_params;
  vector<variable> success_sf;
  vector<variable> success_cr;
  vector<variable> dif_fitness;

  // the contents of M_f and M_cr are all initialiezed 0.5
  vector<variable> memory_sf(memory_size, 0.5);
  vector<variable> memory_cr(memory_size, 0.5);

  variable temp_sum_sf;
  variable temp_sum_cr;
  variable sum;
  variable weight;

  // memory index counter
  int memory_pos = 0;

  // for new parameters sampling
  variable mu_sf, mu_cr;
  int random_selected_period;
  variable *pop_sf = (variable *)malloc(sizeof(variable) * pop_size);
  variable *pop_cr = (variable *)malloc(sizeof(variable) * pop_size);

  // for current-to-pbest/1
  int p_best_ind;
  int p_num = round(pop_size * p_best_rate);
  int *sorted_array = (int *)malloc(sizeof(int) * pop_size);
  Fitness *temp_fit = (Fitness *)malloc(sizeof(Fitness) * pop_size);

  // for linear population size reduction
  int max_pop_size = pop_size;
  int min_pop_size = 4;
  int plan_pop_size;

  int cmp_count = 0, succ_count = 0;
  double mean_gen_cost = 0.0, best_gen_cost = 1.0e+10, costs_sum = 0.0;

  std::ofstream success, clusters_stats, gen_costs, conv_speed;

  // success.open("stats/success_rate", std::ofstream::out | std::ofstream::app);
  // gen_costs.open("stats/gen_costs", std::ofstream::out | std::ofstream::app);
  // clusters_stats.open("stats/clusters_stats", std::ofstream::out | std::ofstream::app);
  // conv_speed.open("stats/conv_speed", std::ofstream::out | std::ofstream::app);

  // data mining aux. structures
  map<int, double> _pattern;
  vector<map<int, double>> patterns;
  //vector<tuple<Individual, double>> elite;

  string basepath = results_root_path;
  int currentPatternIndex = 0;

  //elite_max_size = 30;
  
  // main loop
  elite.clear();
  int restart_gen = 10;
  int decrement = 300;
  int restarts_count = 0;
  double x = 1;

  patterns_count = 0;
  patterns_usage_count = 0;
  vector<interval_pattern> ipatterns;

#ifdef PRINT_STEP_MARKS
cout << "before while" << endl;
#endif

  while (nfes < max_num_evaluations)
  {
    // if (generation > dm_gen_step) 
    //   break;

    if (bsf_fitness - optimum < 10e-8) {
      bsf_fitness = optimum;
      break;
    }

    this->used_gen_count++;

#ifdef PRINT_STEP_MARKS
cout << "elite handling" << endl;
#endif

    for (int i = 0; i < pop_size; i++) sorted_array[i] = i;
    for (int i = 0; i < pop_size; i++) temp_fit[i] = fitness[i];
    sortIndexWithQuickSort(&temp_fit[0], 0, pop_size - 1, sorted_array);
    patterns.clear();

    vector<vector<double>> dataset;
    if (elite_max_size > 0) {
      // update elite set: ========================
      if (elite_type == CROSS_GENERATION) {
        updateElite(pop, sorted_array, temp_fit);
      } else { // BY_GENERATION
        elite.clear();
        for (int i=0; i<p_num; i++) {
          elite.push_back({ pop[sorted_array[i]], fitness[sorted_array[i]] });
        }
      }
      // ==========================================
    }
    else {
      // elite_transactions.clear();
      // for (int i = 0; i < pop_size; i++) {
      //   Transaction_ t;

      //   for (int j = 0; j < problem_size; j++)
      //     t.push_back(Item{ to_string(j) + "-" + to_string( (int)((pop[i][j] - min_region) / discretization_step) ) });

      //   elite_transactions.push_back(t);
      // }
    }
      
#ifdef PRINT_STEP_MARKS
cout << "pattern mining" << endl;
#endif

    double rho_nfes = 0.0;
    bool print_debug = false;
    if (elite.size()) {
      if (generation >= 1 && generation % dm_gen_step == 0) {


        arma::mat data(g_problem_size, elite.size()); // n_rows, n_cols
        // The assignments will be stored in this vector.
        for(int i=0;i<elite.size();++i){
            for(int j=0;j<g_problem_size;++j){
                data(j,i) = get<0>(elite[i])[j];//(i,j) at the i-th row and j-th column
            }
        }

        // cout << data.size() << endl;

        arma::Row<size_t> assignments;
        arma::Row<size_t> bestAssignments;
        arma::mat centroids;
        arma::mat bestCentroids;

        KMeans<> kmeans(10);
        kmeans.Cluster(data, number_of_patterns, assignments, centroids);
        
        if (centroids.size() > 0)
        {
          // cout << centroids.n_rows << endl;

          /* code */
          int pc = min((int)centroids.n_cols, pop_size);
          for (int k = 0; k < pc; k++)
          {
            int idx = sorted_array[pop_size - 1 - k];
            for (size_t j = 0; j < g_problem_size; j++)
              pop[idx][j] = centroids(j, k);
          }
        }
      }
    }

    if (print_debug) cout << "BEST COST = " << bsf_fitness - optimum << endl; 
    if (print_debug) cout << "----" << endl;
    // ========================

#ifdef PRINT_STEP_MARKS
cout << "children generation" << endl;
#endif

    for (int target = 0; target < pop_size; target++)
    {
      // In each generation, CR_i and F_i used by each individual ax_pos are generated by first selecting an index r_i randomly from [1, H]
      random_selected_period = rand() % memory_size;
      mu_sf = memory_sf[random_selected_period];
      mu_cr = memory_cr[random_selected_period];

      // generate CR_i and repair its value
      if (mu_cr == -1) {
        pop_cr[target] = 0;
      }
      else
      {
        pop_cr[target] = gauss(mu_cr, 0.1);
        if (pop_cr[target] > 1) pop_cr[target] = 1;
        else if (pop_cr[target] < 0) pop_cr[target] = 0;
      }

      // generate F_i and repair its value
      do
      {
        pop_sf[target] = cauchy_g(mu_sf, 0.1);
      } while (pop_sf[target] <= 0);

      if (pop_sf[target] > 1) pop_sf[target] = 1;

      /**
       * Strategies:
       *  1. pattern and pbest with same id
       *  2. pattern getted from queue
       *  3. pattern and pbest with random(different) ids
      */

      // cout << 1 << endl;
      unsigned r = rand(); // LEMBRAR: a distribuição de random numbers influencia
      p_best_ind = sorted_array[r % p_num];

      operateCurrentToPBest1BinWithArchive(pop, &children[target][0], target, p_best_ind, pop_sf[target], pop_cr[target], archive, arc_ind_count);
    }

#ifdef PRINT_STEP_MARKS
cout << "generation alternation" << endl;
#endif

    generation++;
    if (debug_mode) fprintPopulation(pop, generation, results_root_path + "/generations/KMEANS_LSHADE/");
    // evaluate the children's fitness values
    evaluatePopulation(children, children_fitness);

    /////////////////////////////////////////////////////////////////////////
    // update the bsf-solution and check the current number of fitness evaluations
    //  if the current number of fitness evaluations over the max number of fitness evaluations, the search is terminated
    //  So, this program is unconcerned about L-SHADE algorithm directly
    for (int i = 0; i < pop_size; i++)
    {
      nfes++;
      mean_gen_cost += fitness[i];

      // following the rules of CEC 2014 real parameter competition,
      // if the gap between the error values of the best solution found and the optimal solution was 10^{−8} or smaller,
      // the error was treated as 0
      if ((children_fitness[i] - optimum) < epsilon) children_fitness[i] = optimum;

      if (children_fitness[i] < bsf_fitness) {
        bsf_fitness = children_fitness[i];
        for (int j = 0; j < problem_size; j++) bsf_solution[j] = children[i][j];
      }

      if (nfes >= max_num_evaluations) break;
    }

    mean_gen_cost /= pop_size;
    //gen_costs << mean_gen_cost << "," << bsf_fitness << endl;
    ////////////////////////////////////////////////////////////////////////////

    // generation alternation
    int chosed_pattern_idx=0;
    for (int i = 0; i < pop_size; i++) { 
      cmp_count++;
      if (children_fitness[i] == fitness[i]) {
        fitness[i] = children_fitness[i];
        for (int j = 0; j < problem_size; j++) pop[i][j] = children[i][j];
      }
      else if (children_fitness[i] < fitness[i]) { 
        succ_count++;

        dif_fitness.push_back(fabs(fitness[i] - children_fitness[i]));
        fitness[i] = children_fitness[i];
        for (int j = 0; j < problem_size; j++) pop[i][j] = children[i][j];
        
        // successful parameters are preserved in S_F and S_CR
        success_sf.push_back(pop_sf[i]);
        success_cr.push_back(pop_cr[i]);

        // parent vectors ax_pos which were worse than the trial vectors u_i are preserved
        if (arc_size > 1) {
          if (arc_ind_count < arc_size) {
            for (int j = 0; j < problem_size; j++) archive[arc_ind_count][j] = pop[i][j];
            arc_ind_count++;
          }
          // Whenever the size of the archive exceeds, randomly selected elements are deleted to make space for the newly inserted elements
          else
          {
            random_selected_arc_ind = rand() % arc_size;
            for (int j = 0; j < problem_size; j++) archive[random_selected_arc_ind][j] = pop[i][j];
          }
        }
        
      }
    }

#ifdef PRINT_STEP_MARKS
cout << "parameter adaptation" << endl;
#endif

    num_success_params = success_sf.size();
    //success << num_success_params << endl;

    // if numeber of successful parameters > 0, historical memories are updated
    if (num_success_params > 0)
    {
      memory_sf[memory_pos] = 0;
      memory_cr[memory_pos] = 0;
      temp_sum_sf = 0;
      temp_sum_cr = 0;
      sum = 0;

      for (int i = 0; i < num_success_params; i++)
        sum += dif_fitness[i];

      // weighted lehmer mean
      for (int i = 0; i < num_success_params; i++)
      {
        weight = dif_fitness[i] / sum;

        memory_sf[memory_pos] += weight * success_sf[i] * success_sf[i];
        temp_sum_sf += weight * success_sf[i];

        memory_cr[memory_pos] += weight * success_cr[i] * success_cr[i];
        temp_sum_cr += weight * success_cr[i];
      }

      memory_sf[memory_pos] /= temp_sum_sf;

      if (temp_sum_cr == 0 || memory_cr[memory_pos] == -1)
        memory_cr[memory_pos] = -1;
      else
        memory_cr[memory_pos] /= temp_sum_cr;

      // increment the counter
      memory_pos++;
      if (memory_pos >= memory_size)
        memory_pos = 0;

      // clear out the S_F, S_CR and delta fitness
      success_sf.clear();
      success_cr.clear();
      dif_fitness.clear();
    }

    // calculate the population size in the next generation
    plan_pop_size = round((((min_pop_size - max_pop_size) / (double)max_num_evaluations) * nfes) + max_pop_size);

    if (pop_size > plan_pop_size)
    {
      reduction_ind_num = pop_size - plan_pop_size;
      if (pop_size - reduction_ind_num < min_pop_size)
        reduction_ind_num = pop_size - min_pop_size;

      reducePopulationWithSort(pop, fitness);

      // resize the archive size
      arc_size = pop_size * g_arc_rate;
      if (arc_ind_count > arc_size) arc_ind_count = arc_size;

      // resize the number of p-best individuals
      p_num = round(pop_size * p_best_rate);
      if (p_num <= 1)
        p_num = 2;
    }

    costs_sum += bsf_fitness - optimum;
    if (debug_mode) {
      for (int j = 0; j < initial_pop_size; j++) {
        if (j < pop_size) {
          pop_costs_file << fitness[j] - optimum;
        } else {
          pop_costs_file << -1;
        }

        if (j < initial_pop_size - 1)
            pop_costs_file << ";";
      }
      pop_costs_file << endl;
    }
    
  }

#ifdef PRINT_STEP_MARKS
cout << "after while" << endl;
#endif

  if (debug_mode)
    pop_costs_file.close();

  //success.close();
  //gen_costs.close();
  // conv_speed << endl;
  // conv_speed.close();

  this->clusters_count_avg /= used_gen_count;

  if (debug_mode) {
      cout << "min. k = " << this->clusters_count_min << endl;
      cout << "avg. k = " << this->clusters_count_avg << endl;
      cout << "max. k = " << this->clusters_count_max << endl;
      cout << "G. count = " << generation << endl;
      cout << "cost = " << bsf_fitness - optimum << endl;
  }
   
  //cout << "Qtde de padrões = " << patterns_count/patterns_usage_count << endl;
  if (patterns_usage_count)
    patterns_count_avg = patterns_count/patterns_usage_count;

  this->success_rate = (double)succ_count / cmp_count;
  this->run_avg_cost = (double)costs_sum / generation;
  return bsf_fitness - optimum;
}

void KMEANS_LSHADE::updateElite(vector<Individual> & curr_pop, int* sorted_indexes, double* fitness)
{
    bool has_update = false;
    int max = elite_max_size > curr_pop.size() ? curr_pop.size() : elite_max_size;

    if (elite.size() == 0) {
      for (int i = 0; i < max; i++) {
        Individual ind = curr_pop[sorted_indexes[i]];
        elite.push_back({ ind, fitness[i] });
      }
    } 
    else 
    {
      for (int i = 0; i < max; i++)
      {
        std::vector<tuple<Individual, double>>::iterator elite_member = elite.end();
        std::vector<set<int>>::iterator elite_t_member = elite_transactions.end();

        int pos = elite.size();
        while (elite_member != elite.begin() && fitness[i] < get<1>(*(elite_member-1))) {
            elite_member--;
            elite_t_member--;
            pos--;
        }

        if (pos < elite.size()) {
          int index = sorted_indexes[i];
          Individual ind = curr_pop[index];
          Fitness fit = fitness[i];

          elite.insert(elite_member, { ind , fit });

          has_update = true;
        } 
        else if (pos < elite_max_size)
        {
          Individual ind = curr_pop[sorted_indexes[i]];
          elite.push_back({ ind , fitness[i] });
        }
     }

    }

   
    if (has_update && elite.size() > elite_max_size) {
        elite.resize(elite_max_size);
        elite_transactions.resize(elite_max_size);
    }
}


void KMEANS_LSHADE::operateCurrentToPBest1BinWithArchive(const vector<Individual> &pop, Individual child, int &target, int &p_best_individual, variable &scaling_factor, variable &cross_rate, const vector<Individual> &archive, int &arc_ind_count)
{
  int r1, r2;

  do
  {
    r1 = rand() % pop_size;
  } while (r1 == target);
  do
  {
    r2 = rand() % (pop_size + arc_ind_count);
  } while ((r2 == target) || (r2 == r1));

  int random_variable = rand() % problem_size;

  if (r2 >= pop_size)
  {
    r2 -= pop_size;
    for (int i = 0; i < problem_size; i++)
    {
      if ((randDouble() < cross_rate) || (i == random_variable))
      {
        child[i] = pop[target][i] + scaling_factor * (pop[p_best_individual][i] - pop[target][i]) + scaling_factor * (pop[r1][i] - archive[r2][i]);
      }
      else
      {
        child[i] = pop[target][i];
      }
    }
  }
  else
  {
    for (int i = 0; i < problem_size; i++)
    {
      if ((randDouble() < cross_rate) || (i == random_variable))
      {
        child[i] = pop[target][i] + scaling_factor * (pop[p_best_individual][i] - pop[target][i]) + scaling_factor * (pop[r1][i] - pop[r2][i]);
      }
      else
      {
        child[i] = pop[target][i];
      }
    }
  }

  // If the mutant vector violates bounds, the bound handling method is applied
  modifySolutionWithParentMedium(child, pop[target]);
}


void KMEANS_LSHADE::operateCurrentToPBest1BinWithArchiveAndXPatternCross(const vector<Individual>& pop, Individual child, int& target, int& p_best_individual,
    variable& scaling_factor, variable& cross_rate, const vector<Individual>& archive, int& arc_ind_count, map<int, double> pattern) {
    int r1, r2;

    set<int> unfixed_positions;

    do {
        r1 = rand() % pop_size;
    } while (r1 == target);
    do {
        r2 = rand() % (pop_size + arc_ind_count);
    } while ((r2 == target) || (r2 == r1));

    int random_variable = rand() % problem_size;

    if (r2 >= pop_size) 
    {
        r2 -= pop_size;
        for (int i = 0; i < problem_size; i++) {
            if ((randDouble() < cross_rate) || (i == random_variable)) {
              child[i] = pop[target][i] + scaling_factor * (pop[p_best_individual][i] - pop[target][i]) + scaling_factor * (pop[r1][i] - archive[r2][i]);
            }
            else {
              child[i] = pop[target][i];
            }
        }
    }
    else 
    {
        for (int i = 0; i < problem_size; i++) {
            if ((randDouble() < cross_rate) || (i == random_variable)) {
              child[i] = pop[target][i] + scaling_factor * (pop[p_best_individual][i] - pop[target][i]) + scaling_factor * (pop[r1][i] - pop[r2][i]);
            }
            else {
              child[i] = pattern[i];
            }
        }
    }

    //If the mutant vector violates bounds, the bound handling method is applied
    modifySolutionWithParentMedium(child, pop[target]);
}


void KMEANS_LSHADE::reducePopulationWithSort(vector<Individual> &pop, vector<Fitness> &fitness)
{
  int worst_ind;

  for (int i = 0; i < reduction_ind_num; i++)
  {
    worst_ind = 0;
    for (int j = 1; j < pop_size; j++)
    {
      if (fitness[j] > fitness[worst_ind])
        worst_ind = j;
    }

    pop.erase(pop.begin() + worst_ind);
    fitness.erase(fitness.begin() + worst_ind);
    pop_size--;
  }
}

void KMEANS_LSHADE::setSHADEParameters()
{
  arc_rate = g_arc_rate;
  arc_size = (int)round(pop_size * arc_rate);
  p_best_rate = g_p_best_rate;
  memory_size = g_memory_size;
}



