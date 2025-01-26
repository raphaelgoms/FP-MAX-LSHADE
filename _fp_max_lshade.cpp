/*
  L-SHADE implemented by C++ for Special Session & Competition on Real-Parameter Single Objective Optimization at CEC-2014

  Version: 1.0   Date: 16/Apr/2014
  Written by Ryoji Tanabe (rt.ryoji.tanabe [at] gmail.com)
*/

#include "de.h"
#include<algorithm>
using namespace std;

FP_MAX_LSHADE::FP_MAX_LSHADE(int cec_function_number, 
  PatternSelectionStrategy pattern_sel_strategy, 
  SolutionFillingStrategy filling_strategy)
{
  this->function_number = cec_function_number;
  this->pattern_sel_strategy = pattern_sel_strategy;
  this->filling_strategy = filling_strategy;
}


Fitness FP_MAX_LSHADE::run()
{
  //cout << this->seed << endl;
  srand(this->seed);
  debug_mode = false;
  double s_CPU_inicial, s_CPU_final;
  double s_total_inicial, s_total_final;

  Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
  
  new_elite = vector<tuple<Individual, double>>(elite_max_size);
  points = vector<Individual>(elite_max_size);
  clock_t time = clock();

  this->used_gen_count = 0;
  initializeParameters();
  setSHADEParameters();

  int generation = 1;

  vector<Individual> pop;
  vector<Fitness> fitness(pop_size, 0);

  vector<Individual> children;
  vector<Fitness> children_fitness(pop_size, 0);

  int initial_pop_size = pop_size;

  // initialize population
  for (int i = 0; i < pop_size; i++)
  {
    pop.push_back(makeNewIndividual());
    children.push_back((variable *)malloc(sizeof(variable) * problem_size));
  }

  // evaluate the initial population's fitness values
  evaluatePopulation(pop, fitness);

  elite.clear();
  elite_transactions.clear();


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
  double mean_gen_cost = 0, best_gen_cost = 1.0e+10, costs_sum = 0.0;

  std::ofstream success, clusters_stats, gen_costs, conv_speed;

  std::string results_root_path = "";
  std::ofstream pop_costs_file;
  std::ofstream dm_effect_data_file;

  if (debug_mode) {
    results_root_path = "/home/raphael/Code/FP-MAX-LSHADE/results/generations/FP-MAX-LSHADE/"
                      + config 
                      + "/D" + to_string(g_problem_size) 
                      + "/cec-f" + to_string(function_number);

    int check = mkdir("/home/raphael/Code/FP-MAX-LSHADE/results/generations",0777);
    check = mkdir("/home/raphael/Code/FP-MAX-LSHADE/results/generations/FP-MAX-LSHADE",0777);
    check = mkdir(("/home/raphael/Code/FP-MAX-LSHADE/results/generations/FP-MAX-LSHADE/"+ config).c_str(),0777);
    check = mkdir(("/home/raphael/Code/FP-MAX-LSHADE/results/generations/FP-MAX-LSHADE/"+ config + "/D" + to_string(g_problem_size)).c_str(),0777);
    
    check = mkdir(results_root_path.c_str(),0777);

    results_root_path += "/s" + to_string(seed);
    check = mkdir(results_root_path.c_str(),0777);
    check = mkdir((results_root_path +"/elite").c_str(),0777); // elite
    check = mkdir((results_root_path +"/population" ).c_str(),0777); // pop

    string fullpath = results_root_path + "/pop-costs.csv";
    pop_costs_file.open(fullpath);

    for (int j = 0; j < pop_size; j++) {
        pop_costs_file << "i" << j+1;

        if (j < pop_size - 1)
          pop_costs_file << ";";
    }

    pop_costs_file << endl;

    fullpath = results_root_path + "/dm-effect.csv";
    dm_effect_data_file.open(fullpath);
    dm_effect_data_file << "G;PD;NP;APS;MSR;PSR" << endl; 
  }

    int j=0;
  // main loop
  while (nfes < max_num_evaluations)
  //while (true)
  {
    // if ( (double)(clock()-time)/CLOCKS_PER_SEC >= 1) {
    //   break;
    // }

    if (debug_mode)
    {
      //fprintPopulation(pop, generation, results_root_path);
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

      dm_effect_data_file << generation << ";" << computeDiversity(pop) << ";";
    }

    if (bsf_fitness - optimum < 10e-8) {
      bsf_fitness = optimum;
      break;
    }

    this->used_gen_count++;

    for (int i = 0; i < pop_size; i++)
      sorted_array[i] = i;
    for (int i = 0; i < pop_size; i++)
      temp_fit[i] = fitness[i];
    sortIndexWithQuickSort(&temp_fit[0], 0, pop_size - 1, sorted_array);

    updateElite(pop, sorted_array, temp_fit);

    // MINING //////////////////////
    bool hasPatterns = false;
    if (pattern_usage_strategy == INSERT_IN_POP 
        && (generation > 0 && generation % dm_gen_step == 0)) {  

      elite_transactions.clear();
      for (int i = 0; i < elite.size(); i++) {
        elite_transactions.push_back(transformToTransaction(get<0>(elite[i])));
      }
      //cout << elite_transactions.size() << endl;

      std::set<Pattern> patterns = computeFrequentIntervalSets(support, discretization_step);
      vector<interval_pattern> ipatterns = computePatternBounds(patterns, discretization_step);   
      
      int insertionsCount = 0;
      double avgPatternSize = 0;
      succ_count = 0;
      cmp_count = 0;

      // std::sort(ipatterns.begin(), ipatterns.end(), 
      //   [](interval_pattern &p1, interval_pattern &p2){ 
      //       return p1.size() > p2.size();
      //   });

      for (auto ipatt : ipatterns) {

        // if (insertionsCount>g_pop_size)
        //   break;
        // cout << "pcr = " << patterns_count_rate << endl;
        
        if (insertionsCount >= std::round(pop_size * patterns_count_rate))
          break;

        if (ipatt.size() < std::round(problem_size * patterns_size_rate))
          continue;

        //cout << "psr = " << patterns_size_rate << endl; 
        // avgPatternSize += ipatt.size();

        int idx = sorted_array[g_pop_size - 1 - insertionsCount++];
        //int idx = sorted_array[rand() % p_num]; insertionsCount++;
        Individual new_ind = makeNewIndividualFromPattern(ipatt);

        // cmp_count++;
        // if(debug_mode && evaluateIndividual(new_ind) < evaluateIndividual(pop[idx]))
        //   succ_count++;

        pop[idx] = new_ind;

        // if (insertionsCount > g_pop_size * 0.8)
        //   break;
      }

      if (debug_mode && ipatterns.size()) {
        hasPatterns = true;

        dm_effect_data_file << cmp_count << ";";
        dm_effect_data_file << avgPatternSize / ipatterns.size() << ";";
        dm_effect_data_file << (double)succ_count / cmp_count << ";";
      }

    }

    ////////////////////////////////

    if (debug_mode && !hasPatterns) dm_effect_data_file << ";;;";

    for (int target = 0; target < pop_size; target++)
    {
      // In each generation, CR_i and F_i used by each individual x_i are generated by first selecting an index r_i randomly from [1, H]
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
        if (pop_cr[target] > 1)
          pop_cr[target] = 1;
        else if (pop_cr[target] < 0)
          pop_cr[target] = 0;
      }

      // generate F_i and repair its value
      do
      {
        pop_sf[target] = cauchy_g(mu_sf, 0.1);
      } while (pop_sf[target] <= 0);

      if (pop_sf[target] > 1)
        pop_sf[target] = 1;

      // p-best individual is randomly selected from the top pop_size *  p_i members
      p_best_ind = sorted_array[rand() % p_num];
      operateCurrentToPBest1BinWithArchive(pop, &children[target][0], target, p_best_ind, pop_sf[target], pop_cr[target], archive, arc_ind_count);
    }

    generation++;
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
      if ((children_fitness[i] - optimum) < epsilon)
        children_fitness[i] = optimum;

      if (children_fitness[i] < bsf_fitness)
      {
        bsf_fitness = children_fitness[i];
        for (int j = 0; j < problem_size; j++)
          bsf_solution[j] = children[i][j];
      }

      if (nfes >= max_num_evaluations)
        break;
    }

    mean_gen_cost /= pop_size;
    //gen_costs << mean_gen_cost << "," << bsf_fitness << endl;
    ////////////////////////////////////////////////////////////////////////////

    cmp_count = 0;
    succ_count = 0;
    // generation alternation
    for (int i = 0; i < pop_size; i++)
    {
      cmp_count++;
      if (children_fitness[i] == fitness[i])
      {
        fitness[i] = children_fitness[i];
        for (int j = 0; j < problem_size; j++)
          pop[i][j] = children[i][j];
      }
      else if (children_fitness[i] < fitness[i])
      {
        succ_count++;
        dif_fitness.push_back(fabs(fitness[i] - children_fitness[i]));
        fitness[i] = children_fitness[i];
        for (int j = 0; j < problem_size; j++)
          pop[i][j] = children[i][j];

        // successful parameters are preserved in S_F and S_CR
        success_sf.push_back(pop_sf[i]);
        success_cr.push_back(pop_cr[i]);

        // parent vectors x_i which were worse than the trial vectors u_i are preserved
        if (arc_size > 1)
        {
          if (arc_ind_count < arc_size)
          {
            for (int j = 0; j < problem_size; j++)
              archive[arc_ind_count][j] = pop[i][j];
            arc_ind_count++;
          }
          // Whenever the size of the archive exceeds, randomly selected elements are deleted to make space for the newly inserted elements
          else
          {
            random_selected_arc_ind = rand() % arc_size;
            for (int j = 0; j < problem_size; j++)
              archive[random_selected_arc_ind][j] = pop[i][j];
          }
        }
      }
    }


    if (debug_mode) {
      dm_effect_data_file <<  (double)succ_count / cmp_count << endl;
    }

    num_success_params = success_sf.size();

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
      if (arc_ind_count > arc_size)
        arc_ind_count = arc_size;

      // resize the number of p-best individuals
      p_num = round(pop_size * p_best_rate);
      if (p_num <= 1)
        p_num = 2;
    }

    costs_sum += bsf_fitness - optimum;

  }

  if (debug_mode) { 
    pop_costs_file.close();
    dm_effect_data_file.close();
  }

  this->success_rate = (double)succ_count / cmp_count;
  this->run_avg_cost = (double)costs_sum / generation;

  return bsf_fitness - optimum;
}

void FP_MAX_LSHADE::operateCurrentToPBest1BinWithArchive(const vector<Individual> &pop, Individual child, int &target, int &p_best_individual, variable &scaling_factor, variable &cross_rate, const vector<Individual> &archive, int &arc_ind_count)
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

void FP_MAX_LSHADE::reducePopulationWithSort(vector<Individual> &pop, vector<Fitness> &fitness)
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

void FP_MAX_LSHADE::setSHADEParameters()
{
  arc_rate = g_arc_rate;
  arc_size = (int)round(pop_size * arc_rate);
  p_best_rate = g_p_best_rate;
  memory_size = g_memory_size;
}

void FP_MAX_LSHADE::updateElite(vector<Individual> & curr_pop, int* sorted_indexes, double* fitness)
{
    // bool has_update = false;
    // int max = elite_max_size > curr_pop.size() ? curr_pop.size() : elite_max_size;

    // if (elite.size() == 0) {
    //   for (int i = 0; i < max; i++) {
    //     Individual ind = curr_pop[sorted_indexes[i]];
    //     // elite.push_back({ ind, fitness[i] });
    //     elite_transactions.push_back(transformToTransaction(ind));
    //   }
    // } 
    // else 
    // {
    //   for (int i = 0; i < max; i++)
    //   {
    //     std::vector<tuple<Individual, double>>::iterator elite_member = elite.end();
    //     std::vector<set<int>>::iterator elite_t_member = elite_transactions.end();

    //     int pos = elite_transactions.size();
    //     while (elite_member != elite.begin() && fitness[i] < get<1>(*(elite_member-1))) {
    //         // elite_member--;
    //         elite_t_member--;
    //         pos--;
    //     }

    //     if (pos < elite_transactions.size()) {

    //       int index = sorted_indexes[i];
    //       Individual ind = curr_pop[index];
    //       Fitness fit = fitness[i];

    //       // elite.insert(elite_member, { ind , fit });
    //       elite_transactions.insert(elite_t_member, transformToTransaction(ind));

    //       has_update = true;
    //     } 
    //     else if (pos < elite_max_size)
    //     {
    //       Individual ind = curr_pop[sorted_indexes[i]];
    //       // elite.push_back({ ind , fitness[i] });
    //       elite_transactions.push_back(transformToTransaction(ind));
    //     }
    //  }

    // }
  
    // if (has_update && elite_transactions.size() > elite_max_size) {
    //     // elite.resize(elite_max_size);
    //     elite_transactions.resize(elite_max_size);
    // }

    int i=0, j=0, k=0;
    //cout << elite_max_size << endl;
    while (true) {
      if (j >= elite_max_size)
        break;

      if (k >= elite_max_size)
        break; 

      // cout<< i << " " << j << " "<< k << endl;
      while (i < elite.size() && k < elite_max_size && get<1>(elite[i]) <= fitness[sorted_indexes[j]]) {
        //points[k] = get<0>(elite[i]);
        new_elite[k] = elite[i];
        k++; i++;
      }

      if (k >= elite_max_size)
        break; 

      //points[k] = curr_pop[sorted_indexes[j]];
      new_elite[k] = { curr_pop[sorted_indexes[j]], fitness[sorted_indexes[j]] };
      k++; j++;
    }

    // cout << i << " " << j << " "<< k << endl;
    elite = new_elite;
}

set<int> FP_MAX_LSHADE::transformToTransaction(Individual ind) {
  set<int> t;

  for (int j = 0; j < problem_size; j++) {
    Item itm = to_string(j) + "-" + to_string( (int)((ind[j] - min_region) / discretization_step) );
    
    if (!mapIntervalToItemID.count(itm)) {
      mapIntervalToItemID[itm] = itemID;
      mapItemIDToInterval[itemID] = itm;
      itemID++;
    }

    t.insert(mapIntervalToItemID[itm]);
  }

  return t;
}

std::set<Pattern> FP_MAX_LSHADE::computeFrequentIntervalSets(int support, double discretizationStep)
{
  std::set<Pattern> patterns;
  Dataset* dataset = new Dataset;
	for(set<int> transaction: elite_transactions)
    dataset->push_back(transaction);

	FISet* freq_itemsets = fpmax(dataset, support);  
  for (FISet::iterator it=freq_itemsets->begin(); it!=freq_itemsets->end(); ++it) {
    set<Item> its;
    for (set<int>::iterator it2=it->begin(); it2!=it->end(); ++it2)
      its.insert(mapItemIDToInterval[*it2]);
    patterns.insert({ its, it->support() });
  }

  delete dataset;
  delete freq_itemsets;

  return patterns;
}


vector<interval_pattern> FP_MAX_LSHADE::computePatternBounds(const std::set<Pattern>& patterns, double discretizationStep)
{
    double lb, ub;
    vector<interval_pattern> ipatterns;

    for (auto p : patterns)
    {
        map<int, tuple<double, double>> intervalPattern;
        set<Item> itemset = get<0>(p);
        
        for (Item i : itemset) {
            std::size_t pos = i.find("-");  
            
            string index = i.substr(0, pos);
            string _interval = i.substr(pos+1);

            int attribute = stoi(index);
            int interval = stoi(_interval);
            
            lb = min_region + (interval) * discretizationStep;
            ub = min_region + (interval + 1) * discretizationStep;

            intervalPattern.insert({ attribute, { lb, ub } });
        }

        if (intervalPattern.size())
          ipatterns.push_back(intervalPattern);
    }

    return ipatterns;
}

Individual FP_MAX_LSHADE::makeNewIndividualFromPattern(interval_pattern pattern, Individual base)
{
  Individual newOne = new double[g_problem_size];
  for (size_t k = 0; k < g_problem_size; k++) {
    if (pattern.find(k) != pattern.end()) {
      double lower = get<0>(pattern[k]);
      double upper = get<1>(pattern[k]);
      newOne[k] = ((upper - lower) * randDouble()) + lower;
    } else {
      if (base)
        newOne[k] = base[k];
      else
        newOne[k] = ((upper_bounds[k] - lower_bounds[k]) * randDouble()) + lower_bounds[k];
    }  
  }
  return newOne;
}
