unsigned long bitscat(unsigned int a);
unsigned int bitgat(unsigned long x);
void C2M(unsigned long *c, unsigned int x,unsigned int y, unsigned int z);
void LC2M(unsigned long *c, unsigned int x,unsigned int y, unsigned int z, unsigned int level);
int get_level(unsigned long c);
void M2C(unsigned long c, unsigned int *x,unsigned int *y, unsigned int *z, unsigned int *level);
void assign_level(int level, unsigned long *c);
void key2pos(unsigned long key, REAL *x, unsigned int *level);
void key2cen(unsigned long key, REAL *x, unsigned int *level);
void reorg(struct CPU *cpu, struct PARAM *param);
void nei6(struct CELL *cell, unsigned long neikey[]);
struct CELL * getcell(unsigned long *key,int level,struct CPU *cpu);
int comp(const void * a, const void * b);
void nei27(struct CELL *cell, unsigned long neikey[]);

void reorgpart(struct CPU *cpu, struct PARAM *param);
struct PART * getpart(unsigned long *key,int level,struct CPU *cpu);
void amr_update_key_part(unsigned int level, struct CPU *cpu, struct PARAM *param);
unsigned long pos2key(REAL *pos, unsigned int level);
void update_key_part(unsigned int level, struct CPU *cpu, struct PARAM *param);
