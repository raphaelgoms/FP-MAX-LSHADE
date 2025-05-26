#ifndef MINT_H_
#define MINT_H_

#include <vector>
#include <tuple>
#include <map>

using std::vector;
using std::tuple;
using std::map;

typedef int interval_id;
typedef vector<interval_id> transaction;
typedef vector<tuple<double, double>> mint_pattern;


typedef struct candidate_join {
    mint_pattern h_j;
    mint_pattern h_k;
    double delta_L;
};

class Mint {
 public:
    explicit Mint(vector<vector<double>> dataSet, 
                    vector<double> lb,
                    vector<double> ub, 
                    double h);
    vector<mint_pattern> run();

 private:
    double m_h;
    vector<double> m_lower_bound;
    vector<double> m_upper_bound;
    vector<vector<double>> m_dataSet;

    map<interval_id, tuple<double, double>> id_to_interval_map;

    vector<candidate_join> createCandidates();
    vector<transaction> discretizeData();
    vector<mint_pattern> getElementaryHyperRectangles(vector<mint_pattern> D);
    tuple<candidate_join, double> popLargestGain(vector<candidate_join> C);
    double computeTotalDescriptionLength(vector<transaction> H);
};
#endif  // MINT_H_
