int PoissonSolver(int level,struct PARAM *param, struct CPU *cpu);
void Prolongation(int level, struct PARAM *param, struct CPU *cpu);
void FillDens(unsigned int level, struct CPU *cpu, struct PARAM *param);
void PoissonForce(unsigned int level, struct CPU *cpu, struct PARAM *param);


