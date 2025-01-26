// HybridContinuousGrasp.cpp : Este arquivo contém a função 'main'. A execução do programa começa e termina ali.
//
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "de.h"
#include "nl_shade_lbc.h"
#include "nl_shade_rsp_mid.h"

#include <sys/stat.h>
int exist(const char *name)
{
  struct stat   buffer;
  return (stat (name, &buffer) == 0);
}

using namespace std;

double *OShift,*M,*y,*z,*x_bound;
int *Rn;
int ini_flag=0,n_flag,func_flag,*SS;
FILE *fpt;

double g_optimum[12] = {300, 400, 600, 800, 900, 1800, 2000, 2200, 2300, 2400, 2600, 2700};

// void Tempo_CPU_Sistema(double *seg_CPU_total, double *seg_sistema_total)
// {
//   long seg_CPU, seg_sistema, mseg_CPU, mseg_sistema;
//   struct rusage ptempo;

//   getrusage(0,&ptempo);

//   seg_CPU = ptempo.ru_utime.tv_sec;
//   mseg_CPU = ptempo.ru_utime.tv_usec;
//   seg_sistema = ptempo.ru_stime.tv_sec;
//   mseg_sistema = ptempo.ru_stime.tv_usec;

//  *seg_CPU_total     = (seg_CPU + 0.000001 * mseg_CPU);
//  *seg_sistema_total = (seg_sistema + 0.000001 * mseg_sistema);
// }

bool isEliteTypeCode(const char *str, searchAlgorithmWithMining::EliteType &elite_type)
{
	if (!strcmp("--bygen", str)) {
		elite_type = searchAlgorithmWithMining::EliteType::BY_GENERATION;
		return true;
	} else if (!strcmp("--crossgen", str)) {
		elite_type = searchAlgorithmWithMining::EliteType::CROSS_GENERATION;
		return true;
	} 

	return false;
}

bool isPatternUsageStrategyCode(const char *str, searchAlgorithmWithMining::PatternUsageStrategy &selection_usage_code)
{
	if (!strcmp("--mutation", str)) {
		selection_usage_code = searchAlgorithmWithMining::PatternUsageStrategy::MUTATION_SCHEME;
		return true;
	} else if (!strcmp("--insert", str)) {
		//cout << "INSERT" << endl;
		selection_usage_code = searchAlgorithmWithMining::PatternUsageStrategy::INSERT_IN_POP;
		return true;
	} 

	return false;
}

bool isPatternSelStrategyCode(const char *str, searchAlgorithmWithMining::PatternSelectionStrategy &selection_strategy_code)
{
	if (!strcmp("--random", str)) {
		selection_strategy_code = searchAlgorithmWithMining::PatternSelectionStrategy::RANDOM;
		return true;
	} else if (!strcmp("--queue", str)) {
		selection_strategy_code = searchAlgorithmWithMining::PatternSelectionStrategy::QUEUE;
		return true;
	} 

	return false;
}

bool isFillingStrategyCode(const char *str, searchAlgorithmWithMining::SolutionFillingStrategy &filling_strategy_code)
{
	if (!strcmp("--pbest", str)) {
		filling_strategy_code = searchAlgorithmWithMining::SolutionFillingStrategy::P_BEST;
		return true;
	} else if (!strcmp("--momentum", str)) {
		filling_strategy_code = searchAlgorithmWithMining::SolutionFillingStrategy::MOMENTUM;
		return true;
	} else if (!strcmp("--linesearch", str)) {
		filling_strategy_code = searchAlgorithmWithMining::SolutionFillingStrategy::LINE_SEARCH;
		return true;
	} else if (!strcmp("--rand", str)) {
		filling_strategy_code = searchAlgorithmWithMining::SolutionFillingStrategy::RAND;
		return true;
	} 

	return false;
}

bool isAlgorithmCode(const char *str, string &algorithmName)
{
	if (!strcmp("--dmc", str)) {
		algorithmName = "dmc";
		return true;
	} else if (!strcmp("--xdmc", str)) {
		algorithmName = "xdmc";
		return true;
	} else if (!strcmp("--mxdmc", str)) {
		algorithmName = "mxdmc";
		return true;
	}  else if (!strcmp("--cgrasp", str)) {
		algorithmName = "cgrasp";
		return true;
	} else if (!strcmp("--rmxdmc", str)) {
		algorithmName = "rmxdmc";
		return true;
	} else if (!strcmp("--amxdmc", str)) {
		algorithmName = "amxdmc";
		return true;
	} else if (!strcmp("--imcgrasp", str)) {
		algorithmName = "imcgrasp";
		return true;
	} else if (!strcmp("--lshade", str)) {
		algorithmName = "lshade";
		return true;
	} 
	else if (!strcmp("--dmlshade", str)) {
		algorithmName = "dmlshade";
		return true;
	}
	else if (!strcmp("--fimlshade", str)) {
		algorithmName = "fimlshade";
		return true;
	}
	else if (!strcmp("--fpglshade", str)) {
		algorithmName = "fpglshade";
		return true;
	}
	else if (!strcmp("--afimlshade", str)) {
		algorithmName = "afimlshade";
		return true;
	}
	else if (!strcmp("--fpmaxlshade", str)) {
		algorithmName = "fpmaxlshade";
		return true;
	}
	else if (!strcmp("--kmeanslshade", str)) {
		algorithmName = "kmeanslshade";
		return true;
	}
	else if (!strcmp("--nlshadelbc", str)) {
		algorithmName = "nlshadelbc";
		return true;
	}
	else if (!strcmp("--nlshaderspmid", str)) {
		algorithmName = "nlshaderspmid";
		return true;
	}

	algorithmName = "";
	return false;
}

int g_pop_size;
double g_arc_rate;
int g_memory_size;
double g_p_best_rate;
int g_function_number;
int g_problem_size;
int g_restart_generation;
unsigned int g_max_num_evaluations;
int g_number_of_used_vars = 0;
double g_min_region;
double g_max_region;
double g_elite_rate;

int main(int argc, char **argv)
{

	double elite_rate = 0.0;
	double clusters_rate = 0.0;
	double patterns_size_rate = 0.0;
	double patterns_count_rate = 1.0;
	string seed_type = "timed";

	string algCode;
	string funcName;
	string function_prefix = "CEC2014";
	int iFuncNumb;
	int problemDim = 10;
	int seed = -1;
	int cec_year = 2014;

	string mining_algorithm = "";
	int number_of_patterns_to_mining = 0;
	int eliteSize = 0;
	int dmStartMoment = 1;
	int dmGenStep = 1;
	int dmFreqStrategy; 
	double patternPercentUsed = 0.0;
	double standardDeviation = 0.0;

	double maxtime = 0.0;

	searchAlgorithmWithMining::EliteType elite_type;
	searchAlgorithmWithMining::PatternUsageStrategy pattern_use_strategy;
	searchAlgorithmWithMining::PatternSelectionStrategy pattern_sel_strategy;
  	searchAlgorithmWithMining::SolutionFillingStrategy filling_strategy;

	int support = 10;
	double discretization_step = 0.1;

	bool irace_mode = false;
	bool analysis_mode = false;
	int max_evals;

	int cec_function_number = -1;
	g_restart_generation = 0;
	
	string arguments("");
	bool print = false;

	// g_min_region = -80.0;
	// g_max_region = 80.0;

	g_min_region = -100.0;
	g_max_region = 100.0;

	std::string data_root_path = "";
	while (!isDirectory((data_root_path + "input_data_cec22").c_str()))
		data_root_path += "../";
	data_root_path += "input_data_cec22/Rand_Seeds.txt";
	
	fpt = fopen(data_root_path.c_str(),"r");
	Rn =(int*)malloc(1000*sizeof(int));
	float f;
	for(int i=0;i<1000;i++)
	{
		fscanf(fpt,"%f",&f);
		Rn[i]=f;
	}
	fclose(fpt);

	int i = 0;
	bool alg_found = false;

	while (++i < argc) 
	{
		if (!alg_found && isAlgorithmCode(argv[i], algCode)) {
			alg_found = true;
			if (print) cout << "0: " << algCode << endl;
			arguments +=  string(argv[i]) + " ";
		}
		else
		if (isEliteTypeCode(argv[i], elite_type)) {
			string et=argv[i];
			if (print) cout << "Elite Type: " << et << endl;
			arguments += et + " ";
		}
		else
		if (isPatternUsageStrategyCode(argv[i], pattern_use_strategy)) {
			string puc=argv[i];
			if (print) cout << "Pattern usage: " << puc << endl;
			arguments += puc + " ";
		}
		else
		if (isPatternSelStrategyCode(argv[i], pattern_sel_strategy)) {
			string psc=argv[i];
			if (print) cout << "Pattern seletion: " << psc << endl;
			arguments += psc + " ";
		}
		else
		if (isFillingStrategyCode(argv[i], filling_strategy)) {
			string fsc=argv[i];
			if (print) cout << "filling_strategy: " << fsc << endl;
			arguments += fsc + " ";
		}
		else
		if (!strcmp("-i", argv[i])) {
			funcName = argv[++i];
			std::transform(funcName.begin(), funcName.end(), funcName.begin(), ::toupper);
			//iFuncNumb = getFuncNumb(funcName.c_str());
			if (print) cout << "Instance: " << iFuncNumb << endl;
			//arguments += "-i " + string(argv[i])  + " ";
		}
		else
		if (!strcmp("--nvar", argv[i])) {
			problemDim = atoi(argv[++i]);
			if (print) cout << "Number of variables: " << problemDim << endl;
			//arguments += "--nvar " + string(argv[i])  + " ";
		} 
		else
		if (!strcmp("--seed", argv[i])) {
			seed = atoi(argv[++i]);
			if (print) cout << "Seed: " << seed << endl;
		} 
		else
		if (!strcmp("--elsz", argv[i])) {
			eliteSize = atoi(argv[++i]);
			if (print) cout << "Elite size: " << eliteSize << endl;
			arguments += "--elsz " +  string(argv[i])  + " ";
		}
		else
		if (!strcmp("--gs", argv[i])) {
			dmGenStep = atoi( argv[++i]);
			if (print) cout << "DM generations step: " << dmGenStep << endl;
			arguments += "--gs " +  string(argv[i])  + " ";
		}
		else
		if (!strcmp("--dmstart", argv[i])) {
			dmStartMoment = atoi( argv[++i]);
			if (print) cout << "DM start moment: " << dmStartMoment << endl;
			arguments += "--dmstart " +  string(argv[i])  + " ";
		}
		else
		if (!strcmp("--dmfreq", argv[i])) {
			
			if (!strcmp("jot", argv[++i])) {
				dmFreqStrategy = 0; // just one time
			} else {
				dmFreqStrategy = 1; // continuous
			}

			if (print) cout << "DM frequency: " << argv[i] << endl;
			arguments += "--dmfreq " + string(argv[i])  + " ";
		}
		else
		if (!strcmp("--ptsz", argv[i])) {
			patternPercentUsed = atof( argv[++i]);
			if (print) cout << "Percent of pattern used: " << patternPercentUsed << endl;
			arguments += "--ptsz " + string(argv[i]) + " ";
		}
		else
		if (!strcmp("--sd", argv[i])) {
			
			standardDeviation = atof( argv[++i]);
			if (print) cout << "Standard Deviation used: " << standardDeviation << endl;
			arguments += "--sd " +  string(argv[i])  + " ";
		}
		else
		if (!strcmp("--irace", argv[i])) {
			irace_mode = true;
			arguments += " --irace ";
		}
		else
		if (!strcmp("--a", argv[i])) {
			analysis_mode = true;
			arguments += " --analysis ";
		}
		else
		if (!strcmp("--max_evals", argv[i])) {
			
			max_evals = atoi( argv[++i]);
			if (print) cout << "Max num. of evals: " << standardDeviation << endl;
			arguments += "--max_evals " +  string(argv[i])  + " ";
		}
		else
		if (!strcmp("--cec", argv[i])) {
			cec_function_number = atoi( argv[++i]);
			if (print) cout << "CEC Function: " << cec_function_number << endl;
			//arguments += "--cec " +  string(argv[i])  + " ";
			funcName = function_prefix + "-f" + to_string(cec_function_number) + "-d";
		}
		else
		if (!strcmp("--rg", argv[i])) {
			g_restart_generation = atoi( argv[++i]);
			if (print) cout << "Restart Generation Step: " << g_restart_generation << endl;
			arguments += "--rg " +  string(argv[i])  + " ";
		}
		else
		if (!strcmp("--k", argv[i])) {
			number_of_patterns_to_mining = atoi( argv[++i]);
			if (print) cout << "Number of patterns to mining: " << number_of_patterns_to_mining << endl;
			arguments += "--k " +  string(argv[i])  + " ";
		}
		else
		if (!strcmp("--x", argv[i])) {
			mining_algorithm = "xmeans";
			arguments += "--k " +  string(argv[i])  + " ";
		}
		else
		if (!strcmp("--uv", argv[i])) {
			g_number_of_used_vars = atoi( argv[++i]);
			//arguments += "--uv " +  string(argv[i])  + " ";
		}
		else
		if (!strcmp("--ps", argv[i])) {
			g_number_of_used_vars = atoi( argv[++i]);
			//arguments += "--uv " +  string(argv[i])  + " ";
		}
		else
		if (!strcmp("--s", argv[i])) {
			support = atof( argv[++i]);
			if (print) cout << "Support: " << support << endl;
			arguments += "--s " +  string(argv[i])  + " ";
		}
		else
		if (!strcmp("--h", argv[i])) {
			discretization_step = atof( argv[++i]);
			if (print) cout << "Discretization Step: " << discretization_step << endl;
			arguments += "--h " +  string(argv[i])  + " ";
		}
		else
		if (!strcmp("--y", argv[i])) {
			cec_year = atoi( argv[++i]);
			cout << cec_year << endl;
			if (print) cout << "Discretization Step: " << discretization_step << endl;
			arguments += "--y " +  string(argv[i])  + " ";
		}
		else 
		if (!strcmp("--t", argv[i])) {
			cout << "time getted" << endl;
			maxtime = atof( argv[++i]);
			if (print) cout << "Max time: " << maxtime << endl;
			arguments += "--t ";
		}
		else
		if (!strcmp("--er", argv[i])) {
			elite_rate = atof( argv[++i]);
			if (print) cout << "Elite Rate: " << elite_rate << endl;
			arguments += "--er " +  string(argv[i])  + " ";
		}
		else 
		if (!strcmp("--cr", argv[i])) {
			clusters_rate = atof( argv[++i]);
			if (print) cout << "Clusters Rate: " << clusters_rate << endl;
			arguments += "--cr " +  string(argv[i])  + " ";
		}
		else
		if (!strcmp("--psr", argv[i])) {
			patterns_size_rate = atof( argv[++i]);
			if (print) cout << "Pattern Size Rate: " << patterns_size_rate << endl;
			arguments += "--psr " +  string(argv[i])  + " ";
		}
		else 
		if (!strcmp("--pcr", argv[i])) {
			patterns_count_rate = atof( argv[++i]);
			if (print) cout << "Patterns Count Rate: " << patterns_count_rate << endl;
			arguments += "--pcr " +  string(argv[i])  + " ";
		}  
		else if (!strcmp("--sf", argv[i])) {
			seed_type = "fixed";
			if (print) cout << "Seeds type: " << seed_type << endl;
			arguments += "--sf ";
		} 
		else{
			arguments += string(argv[i])  + " ";
		}
	}

	if (!g_number_of_used_vars) {
		g_number_of_used_vars = problemDim;
	}

	max_evals = problemDim * 10000;
	

	if (!alg_found){ 
		algCode = "cgrasp";
		arguments += algCode + " ";
	}

	// Construct The correct names of function to show in the costs:
	string funcCode = funcName;
	

	int rand_idx;
	int numOfIterations = 100;

	tuple<vector<double>, double> best_solution;
	g_problem_size = problemDim;

	g_pop_size = (int)round(g_problem_size * 18);
	g_memory_size = 6;
	g_arc_rate = 2.6;
	g_p_best_rate = 0.11;
	g_max_num_evaluations = max_evals; //g_problem_size * 1200;
	g_elite_rate = elite_rate;

	if (cec_function_number < 0)
		g_function_number = iFuncNumb;
	else {
		g_function_number = cec_function_number;
		g_min_region = -100;
		g_max_region = 100;
	}

	searchAlgorithm *lshade, *dmlshade, *fimlshade, *fp_growth_lshade, *afim_lshade, *fp_max_lshade, 
		*kmeans_lshade;
	Optimizer OptZ;
	
	
	if (cec_function_number > -1) {
		lshade = new LSHADE(cec_function_number);
		fp_max_lshade = new FP_MAX_LSHADE(cec_function_number);
		kmeans_lshade = new KMEANS_LSHADE(cec_function_number);
	}

	double cost;
	double patterns_count;
	double avg_cost = 0.0;
	double avg_succ_rate = 0.0;
	double avg_time = 0.0;
	double avg_cfos = 0.0;
	double avg_patterns_count = 0.0;
	double avg_used_gens = 0.0;
	double avg_dm_atvs = 0.0;

	double min = 1.0e+30;

	double s_CPU_inicial, s_CPU_final;
  	double s_total_inicial, s_total_final;

	double avg_cgrasp = 0.0, avg_imcgrasp = 0.0, avg_lshade = 0.0, avg_dmlshade = 0.0;

	if (irace_mode)
	{	
		srand(seed);
		if (algCode == "lshade") {
			
			//Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
			((LSHADE*)lshade)->maxtime = maxtime;
			cost = lshade->run();

			//Tempo_CPU_Sistema(&s_CPU_final, & s_total_final);   
		}
		else if (algCode == "fpmaxlshade") {
			//cout << "fpmaxlshade" << endl;
			((FP_MAX_LSHADE*)fp_max_lshade)->number_of_patterns = number_of_patterns_to_mining;
			((FP_MAX_LSHADE*)fp_max_lshade)->mining_algorithm = mining_algorithm;
			((FP_MAX_LSHADE*)fp_max_lshade)->pattern_sel_strategy = pattern_sel_strategy;
			((FP_MAX_LSHADE*)fp_max_lshade)->pattern_usage_strategy = pattern_use_strategy;
			((FP_MAX_LSHADE*)fp_max_lshade)->filling_strategy = filling_strategy;
			((FP_MAX_LSHADE*)fp_max_lshade)->elite_type = elite_type;
			((FP_MAX_LSHADE*)fp_max_lshade)->dm_start_gen = dmStartMoment;
			((FP_MAX_LSHADE*)fp_max_lshade)->dm_gen_step = dmGenStep;
			((FP_MAX_LSHADE*)fp_max_lshade)->config = dmGenStep;

			((FP_MAX_LSHADE*)fp_max_lshade)->patterns_count_rate = patterns_count_rate;
			((FP_MAX_LSHADE*)fp_max_lshade)->patterns_size_rate = patterns_size_rate;
		
			if(eliteSize==0)
				eliteSize = (int)std::round(g_elite_rate * g_pop_size);

			//cout <<" elite size = " << eliteSize << endl;
			((FP_MAX_LSHADE*)fp_max_lshade)->elite_max_size = eliteSize;
			((FP_MAX_LSHADE*)fp_max_lshade)->seed = seed; //Rn[rand_idx];

			float minsup = ((float)support * eliteSize) / 100;
			//cout << "support = " << minsup << endl;

			((FP_MAX_LSHADE*)fp_max_lshade)->support = minsup;
			((FP_MAX_LSHADE*)fp_max_lshade)->discretization_step = discretization_step; //Rn[rand_idx];
			((FP_MAX_LSHADE*)fp_max_lshade)->cec_year = cec_year;
			
			

			//Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
			cost = fp_max_lshade->run();
			//Tempo_CPU_Sistema(&s_CPU_final, & s_total_final);

			//avg_patterns_count += afim_lshade->
		}

		cout << cost;
		return 0;
	}
	else
	{


	int nruns=51; // 30;
	if (analysis_mode) { 
		nruns = 1;
	}

	std::string results_root_path = "";

	while (!isDirectory((results_root_path + "results").c_str()))
		results_root_path += "../";
	results_root_path += "results";

	std::ofstream stats_test_data;
	string dirname = results_root_path + "/stats_test_data/";
	if (algCode == "dmlshade") 
		dirname += "uv=" + to_string(g_number_of_used_vars);

	int check = mkdir(dirname.c_str(),0777);
	string sdfilepath(dirname +  "/out(" + arguments + ").csv");

	if (!exist(sdfilepath.c_str())) {
		stats_test_data.open(sdfilepath);
		
		stats_test_data << "f";
		for (int i = 0; i < nruns; i++) {
			stats_test_data << ";s" << i+1;
		}
		stats_test_data << endl;

	} else {
		stats_test_data.open(sdfilepath, std::ofstream::out | std::ofstream::app);
	}

	std::ofstream run_avg_cost_data;
	dirname = results_root_path + "/run_avg_cost_data/";
	check = mkdir(dirname.c_str(),0777);
	if (algCode != "fimlshade") 
		dirname += "uv=" + to_string(g_number_of_used_vars);

	check = mkdir(dirname.c_str(),0777);
	sdfilepath = dirname +  "/out(" + arguments + ").csv";

	if (!exist(sdfilepath.c_str())) {
		run_avg_cost_data.open(sdfilepath);
		
		run_avg_cost_data << "f";
		for (int i = 0; i < nruns; i++) {
			run_avg_cost_data << ";s" << i+1;
		}
		run_avg_cost_data << endl;

	} else {
		run_avg_cost_data.open(sdfilepath, std::ofstream::out | std::ofstream::app);
	}

	std::ofstream clusters_stats;

	dirname = "stats/"+ funcCode;
	check = mkdir(dirname.c_str(),0777);

	string cs_filepath(dirname + "/(" + arguments + ").csv");
	
	if (!exist(cs_filepath.c_str())) {
		clusters_stats.open(cs_filepath);
		clusters_stats << "min K;avg. K;max K" << endl;

	} else {
		clusters_stats.open(cs_filepath, std::ofstream::out | std::ofstream::app);
	}
	
	stats_test_data << cec_function_number << "-" << g_problem_size;
	run_avg_cost_data << cec_function_number << "-" << g_problem_size;

	int best_solution_exec = -1;
	double best_cost = 1e+50;

	for (int i = 0; i < nruns; i++)
	{	 
		//cout << "--------------" << endl;
		//cout << "RUN " << (i + 1) << ":" << endl;
		//cout << endl;

		//if (i+1==5){
		//	continue;
		//}

		if (seed_type == "fixed") 
		{
			int run_id = i;
			int func_no = g_function_number+1;
			
			rand_idx = (g_problem_size/10*func_no*nruns+run_id)-nruns;
			rand_idx = rand_idx%1000+1;

			seed = Rn[rand_idx];
			srand(Rn[rand_idx]);
		} 
		else { 
			srand(time(NULL));
		}

		if (algCode == "lshade") {

			lshade->seed = seed;
			lshade->cec_year = cec_year;
			lshade->maxtime = maxtime;
			
			if (analysis_mode){ 
				((LSHADE*)lshade)->debug_mode = true;
			}

			Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
			cost = lshade->run();
			Tempo_CPU_Sistema(&s_CPU_final, & s_total_final);   

			avg_succ_rate += lshade->success_rate;
			avg_used_gens += lshade->used_gen_count;
			avg_cfos += lshade->function_calls;
			
			run_avg_cost_data << ";" << lshade->run_avg_cost;  
		}
		else if (algCode == "fpmaxlshade") {
			//cout << seed << endl;
			((FP_MAX_LSHADE*)fp_max_lshade)->seed = seed;
			((FP_MAX_LSHADE*)fp_max_lshade)->number_of_patterns = number_of_patterns_to_mining;
			((FP_MAX_LSHADE*)fp_max_lshade)->mining_algorithm = mining_algorithm;
			((FP_MAX_LSHADE*)fp_max_lshade)->pattern_sel_strategy = pattern_sel_strategy;
			((FP_MAX_LSHADE*)fp_max_lshade)->pattern_usage_strategy = pattern_use_strategy;
			((FP_MAX_LSHADE*)fp_max_lshade)->filling_strategy = filling_strategy;
			((FP_MAX_LSHADE*)fp_max_lshade)->elite_type = elite_type;
			((FP_MAX_LSHADE*)fp_max_lshade)->dm_start_gen = dmStartMoment;
			((FP_MAX_LSHADE*)fp_max_lshade)->dm_gen_step = dmGenStep;

			if(eliteSize==0)
				eliteSize = (int)std::round(g_elite_rate * g_pop_size);

			// cout <<" elite size = " << eliteSize << endl;

			((FP_MAX_LSHADE*)fp_max_lshade)->elite_max_size = eliteSize;
			((FP_MAX_LSHADE*)fp_max_lshade)->seed = seed; //Rn[rand_idx];

			float minsup = ((float)support * eliteSize) / 100;
			//cout << "support = " << minsup << endl;

			((FP_MAX_LSHADE*)fp_max_lshade)->support = minsup;
			((FP_MAX_LSHADE*)fp_max_lshade)->discretization_step = discretization_step; //Rn[rand_idx];
			((FP_MAX_LSHADE*)fp_max_lshade)->cec_year = cec_year;

			if (analysis_mode){ 
				((FP_MAX_LSHADE*)fp_max_lshade)->debug_mode = true;
				((FP_MAX_LSHADE*)fp_max_lshade)->config = arguments;
			}

			((FP_MAX_LSHADE*)fp_max_lshade)->patterns_count_rate = patterns_count_rate;
			((FP_MAX_LSHADE*)fp_max_lshade)->patterns_size_rate = patterns_size_rate;

			Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
			cost = fp_max_lshade->run();
			Tempo_CPU_Sistema(&s_CPU_final, & s_total_final);

			//avg_patterns_count += afim_lshade->
		} 
		else if (algCode == "kmeanslshade") {

			((KMEANS_LSHADE*)kmeans_lshade)->number_of_patterns = number_of_patterns_to_mining;
			((KMEANS_LSHADE*)kmeans_lshade)->mining_algorithm = mining_algorithm;
			((KMEANS_LSHADE*)kmeans_lshade)->pattern_sel_strategy = pattern_sel_strategy;
			((KMEANS_LSHADE*)kmeans_lshade)->pattern_usage_strategy = pattern_use_strategy;
			((KMEANS_LSHADE*)kmeans_lshade)->filling_strategy = filling_strategy;
			((KMEANS_LSHADE*)kmeans_lshade)->elite_type = elite_type;
			((KMEANS_LSHADE*)kmeans_lshade)->dm_start_moment = dmStartMoment;
			((KMEANS_LSHADE*)kmeans_lshade)->dm_gen_step = dmGenStep;

			((KMEANS_LSHADE*)kmeans_lshade)->elite_max_size = eliteSize;
			((KMEANS_LSHADE*)kmeans_lshade)->seed = seed;

			Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
			cost = kmeans_lshade->run();
			Tempo_CPU_Sistema(&s_CPU_final, & s_total_final);

		} 
		else if (algCode == "nlshadelbc") {

			////////////////NInds     	  NVars      func              Run  memory    arch size maxFEvals
			OptZ.Initialize(problemDim*23,problemDim,g_function_number,i,20*problemDim,1,max_evals,seed);
			
			Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
			cost = OptZ.MainCycle();
			Tempo_CPU_Sistema(&s_CPU_final, & s_total_final);

			OptZ.Clean();
		} 
		// else if (algCode == "nlshaderspmid") {
		// 	int dimensionality = problemDim;
		// 	int POP_SIZE=dimensionality*5;
		// 	int oldPopSize=POP_SIZE;

		// 	RSPOptimizer RSPOptZ;
			
		// 	////////////////NInds     	  NVars      func              Run  memory    arch size maxFEvals
		// 	RSPOptZ.Initialize(POP_SIZE,dimensionality,g_function_number,0,20*dimensionality,2.1, max_evals,10);
			
		// 	Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);

		// 	cost = RSPOptZ.Run();

		// 	Tempo_CPU_Sistema(&s_CPU_final, & s_total_final);

		// 	RSPOptZ.Clean();
			
		// }

		avg_cost += cost;
		stats_test_data << ";" << cost;
		avg_time += (s_CPU_final - s_CPU_inicial);
		avg_patterns_count += patterns_count;
		
		//cout << "cost = " << cost << endl;
		//cout << "run time = " << (s_CPU_final - s_CPU_inicial) << endl;
		//cout << "patterns count = " << patterns_count << endl;

		if (cost < min) min = cost;
		
		//cout << "--------------" << endl;
		// cout << endl;
	}
	cout << endl;
	stats_test_data << endl;
	stats_test_data.close();

	run_avg_cost_data << endl;
	run_avg_cost_data.close();

	avg_cost /= nruns;
	avg_time /= nruns;
	avg_cfos /= nruns;
	avg_succ_rate /= nruns;
	avg_used_gens /= nruns;
	avg_patterns_count /= nruns;
	
	if (algCode != "lshade") {
		avg_dm_atvs /= nruns;
	}

	clusters_stats.close();

	cout << "avg cost = " << avg_cost << endl;
	cout << "best solution run = " << best_solution_exec << endl;
	cout << endl;
	
	// Save average cost for a function
	dirname = results_root_path;
	if (algCode != "lshade") 
		dirname += "/uv=" + to_string(g_number_of_used_vars);

	check = mkdir(dirname.c_str(),0777);
	string filepath(dirname +  "/out(" + arguments + ").csv");

	std::ofstream outfile;
	if (!exist(filepath.c_str())) {
		outfile.open(filepath);
		
		if (cec_function_number < 0)
			outfile << "Func;Dim;BestFO;AvgFO;CFOs;Time" << endl;
		else {
			outfile << "Func;Dim;BestFO;AvgFO;CFOs;Time;SR;GC";
			if (algCode != "lshade")
				outfile << ";DMAtv;BR;PC";
			outfile << endl;
		}

	} else {
		outfile.open(filepath, std::ofstream::out | std::ofstream::app);
	}

	if (cec_function_number < 0)
		outfile << to_string(cec_function_number) << ";" << problemDim << ";" << min << ";" << avg_cost << ";" << avg_cfos << ";" << avg_time << endl;
	else {
		outfile << to_string(cec_function_number) << ";" << problemDim << ";" << min << ";" << avg_cost << ";" << avg_cfos << ";" << avg_time << ";" << avg_succ_rate << ";" << avg_used_gens;
		if (algCode != "lshade")
			outfile << ";" << avg_dm_atvs << ";" << best_solution_exec << ";" << avg_patterns_count;
		outfile << endl;
	}
	outfile.close();
	}
}
