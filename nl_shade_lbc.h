
class Optimizer
{
public:
    bool FitNotCalculated;
    int Int_ArchiveSizeParam;
    int MemorySize;
    int MemoryIter;
    int SuccessFilled;
    int MemoryCurrentIndex;
    int NVars;
    int NInds;
    int NIndsMax;
    int NIndsMin;
    int besti;
    int func_num;
    int RunN;
    int Generation;
    int ArchiveSize;
    int CurrentArchiveSize;
    double MWLp1;
    double MWLp2;
    double MWLm;
    double LBC_fin;
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
    double* FitMass;
    double* FitMassTemp;
    double* FitMassCopy;
    double* BestInd;
    double* tempSuccessCr;
    double* tempSuccessF;
    double* FGenerated;
    double* CrGenerated;
    double* MemoryCr;
    double* MemoryF;
    double* FitDelta;
    double* FitMassArch;
    double** Popul;
    double** PopulTemp;
    double** Archive;

    void Initialize(int newNInds, int newNVars, int new_func_num, int newRunN, int NewMemSize, double NewArchSizeParam, int newMaxFEval, uint seed=0);
    void Clean();
    double MainCycle();
    void FindNSaveBest(bool init, int ChosenOne);
    inline double GetValue(const int index, const int NInds, const int j);
    void CopyToArchive(double* RefusedParent,double RefusedFitness);
    void SaveSuccessCrF(double Cr, double F, double FitD);
    void UpdateMemoryCrF();
    double MeanWL_general(double* Vector, double* TempWeights, int Size, double g_p, double g_m);
    void RemoveWorst(int NInds, int NewNInds);
};