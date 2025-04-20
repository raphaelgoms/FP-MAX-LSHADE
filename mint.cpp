#include "./mint.h"

#include <algorithm>
#include <vector>
#include <tuple>

Mint::Mint(vector<vector<double>> dataSet) { }

bool findIn(mint_pattern h, vector<mint_pattern> H) {
    return std::find(H.begin(), H.end(), h) != H.end();
}

vector<mint_pattern> join(vector<mint_pattern> h,
    mint_pattern h_j, mint_pattern h_k) {
    return vector<mint_pattern>();
}

vector<mint_pattern> Mint::run() {
    vector<mint_pattern> D = discretizeData();
    vector<mint_pattern> H = getElementaryHyperRectangles(D);
    double L_total = computeTotalDescriptionLength(H);

    vector< candidate_join > C = createCandidates();  // candidates set
    while (C.size()) {
        auto C_new = vector<candidate_join>();
        tuple<candidate_join, double> candidate_pattern = popLargestGain(C);

        double delta_L = std::get<1>(candidate_pattern);
        while (delta_L > 0) {
            candidate_join cj = std::get<0>(candidate_pattern);
            if (findIn(cj.h_j, H) && findIn(cj.h_k, H)) {
                auto H_ = join(H, cj.h_j, cj.h_k);
                auto L_new = computeTotalDescriptionLength(H_);
                if (L_new <= L_total) {
                    auto new_candidates =
                        vector<candidate_join>();  // TODO(raphael): construct new candidates set
                    C_new.insert(C_new.end(),
                        new_candidates.begin(), new_candidates.end());
                    L_total = L_new;
                }
            }

            tuple<candidate_join, double> candidate_pattern = popLargestGain(C);
        }
        C = C_new;
    }



    return vector<mint_pattern>();
}

vector<candidate_join> Mint::createCandidates() {
    return vector<candidate_join>();
}

vector<mint_pattern> Mint::discretizeData() {
    return vector<mint_pattern>();
}

vector<mint_pattern>
    Mint::getElementaryHyperRectangles(vector<mint_pattern> D) {
    return vector<mint_pattern>();
}

tuple<candidate_join, double> Mint::popLargestGain(vector< candidate_join > C) {
    return tuple<candidate_join, double>();
}

double Mint::computeTotalDescriptionLength(vector<mint_pattern> H) {
    return 0.0;
}
