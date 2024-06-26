#include <fstream>
#include <iostream>

class RSPOptimizer
{
public:
    bool FitNotCalculated;
    int Int_ArchiveSizeParam;
    int MemorySize;
    int MemoryIter;
    int SuccessFilled;
    int MemoryCurrentIndex;
    int NVars;
    int popSize;
    int NIndsMax;
    int NIndsMin;
    int besti;
    int func_num;
    int runNum;
    int Generation;
    int ArchiveSize;
    int CurrentArchiveSize;
    double F;
    double Cr;
    double bestfit;
    double ArchiveSizeParam;
    double Right;
    double Left;
    int* Rands;
    int* Indexes;
    int* BackIndexes;
    double* Weights;
    double* Donor;
    double* Trial;
    double* Fitmass;
    double* popFitTmp;
    double* FitmassCopy;
    double* BestInd;
    double* tempSuccessCr;
    double* tempSuccessF;
    double* FGenerated;
    double* CrGenerated;
    double* MemoryCr;
    double* MemoryF;
    double* FitDelta;
    double* ArchUsages;
    double** Popul;
    double** populTemp;
    double** Archive;
    std::ofstream infoLog;

    void Initialize(int newNInds, int newNVars, int new_func_num, int newrunNum, int NewMemSize, double NewArchSizeParam, int newMaxFEval, uint seed=0);
    void restart(int newNInds, int newNVars, int new_func_num, int newrunNum, int NewMemSize, double NewArchSizeParam);
    void Clean();
    void MainCycle(double);
    void FindNSaveBest(bool init, int ChosenOne);
    inline double GetValue(const int index, const int popSize, const int j);
    void CopyToArchive(double* RefusedParent);
    void SaveSuccessCrF(double Cr, double F, double FitD);
    void UpdateMemoryCrF();
    double MeanWL_general(double* Vector, double* TempWeights, int Size, double g_p, double g_m);
    void RemoveWorst(int popSize, int NewNInds);
    void RemoveTooNear(int popSize, int NewNInds);

    double Run();
};
