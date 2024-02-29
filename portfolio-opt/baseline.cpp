// Baseline with random weights

#include <omp.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

double simulatePortfolioReturn(const std::vector<double>& weights, const std::vector<double>& assetReturns) {
    double portfolioReturn = 0.0;
    for (size_t i = 0; i < weights.size(); ++i) {
        portfolioReturn += weights[i] * assetReturns[i];
    }
    return portfolioReturn;
}

int main() {
    const int numPortfolios = 10000;
    const int numAssets = 5; 
    std::vector<double> assetReturns = {0.05, 0.07, 0.06, 0.04, 0.03}; 
    std::vector<double> portfolioReturns(numPortfolios);

    srand(time(nullptr));

    #pragma omp parallel for
    for (int i = 0; i < numPortfolios; ++i) {
        std::vector<double> weights(numAssets);
        double weightSum = 0.0;

        for (int j = 0; j < numAssets; ++j) {
            weights[j] = rand() / (double)RAND_MAX; // random weights
            weightSum += weights[j];
        }

        for (int j = 0; j < numAssets; ++j) {
            weights[j] /= weightSum;
        }

        portfolioReturns[i] = simulatePortfolioReturn(weights, assetReturns);
    }

    double maxReturn = 0.0;
    for (int i = 0; i < numPortfolios; ++i) {
        if (portfolioReturns[i] > maxReturn) {
            maxReturn = portfolioReturns[i];
        }
    }

    std::cout << "Best return with baseline optimization: " << maxReturn << std::endl;

    return 0;
}
