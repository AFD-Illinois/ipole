#ifndef DICT_H
#define DICT_H

#define DICT_STRLEN (128)

struct of_dictel {
  char key[DICT_STRLEN];
  char value[DICT_STRLEN];
  struct of_dictel *next;
};

typedef struct of_dict {
  int size;
  struct of_dictel *head;
} dict;

dict *dict_new();
void dict_del(dict *d);
void dict_add(dict *d, char *key, char *value);
char *dict_get(dict *d, char *key, char *fallback);
void dict_print(dict *d);

#endif // DICT_H
