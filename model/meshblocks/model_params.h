#ifndef MODEL_PARAMS_H
#define MODEL_PARAMS_H

#include <stdio.h>
#include <stddef.h>

#define SLOW_LIGHT (0)

extern double sigma_cut;

double get_dump_t(char *fnam, int dumpidx);
size_t get_athenak_datastart(char *fname);
void load_athenak_data(int n, char *fnam, int index, size_t datastart);

void strim(char *s);
void get_line_token(char *line, char *delimiter, int index, char **result);
int read_line(FILE *fp, char **line, char *message);
int read_int(FILE *fp, size_t bytes);
double read_double(FILE *fp, size_t bytes);

int is_within_meshblock(size_t mb, double x1, double x2, double x3);
int get_meshblock(double x1, double x2, double x3);

#endif // MODEL_PARAMS_H
