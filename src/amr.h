int crit_phy(struct CELL *cell);
void mark_phy(unsigned int level, struct CPU *cpu, struct PARAM *param,int ismooth);
void mark_child(unsigned int level, struct CPU *cpu, struct PARAM *param,int ismooth);
void mark_nei(unsigned int level, struct CPU *cpu, struct PARAM *param,int ismooth);
void create_cell(unsigned int level, struct CPU *cpu, struct PARAM *param);
void destroy_cell(unsigned int level, struct CPU *cpu, struct PARAM *param);
