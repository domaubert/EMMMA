void cic(unsigned int level, struct CPU *cpu, struct PARAM *param);
void L_resetispart(unsigned int level, struct CPU *cpu, struct PARAM *param);
void L_accelpart(unsigned int level, struct CPU *cpu, struct PARAM *param, REAL* adt, int is);
void L_movepart(unsigned int level, struct CPU *cpu, struct PARAM *param, REAL* adt, int is);
void L_levpart(unsigned int level, struct CPU *cpu, struct PARAM *param,int is);
REAL L_comptstep(unsigned int level, struct CPU *cpu, struct PARAM *param);
