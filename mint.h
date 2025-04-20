#ifndef MINT_H_
#define MINT_H_

#include <vector>
#include <tuple>

using std::vector;
using std::tuple;

typedef vector<tuple<double, double>> mint_pattern;

typedef struct candidate_join {
    mint_pattern h_j;
    mint_pattern h_k;
    double delta_L;
};

class Mint {
 public:
    explicit Mint(vector<vector<double>> dataSet);
    vector<mint_pattern> run();

 private:
    vector<candidate_join> createCandidates();
    vector<mint_pattern> discretizeData();
    vector<mint_pattern> getElementaryHyperRectangles(vector<mint_pattern> D);
    tuple<candidate_join, double> popLargestGain(vector<candidate_join> C);
    double computeTotalDescriptionLength(vector<mint_pattern> H);
};
#endif  // MINT_H_
