#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "dict.h"

// this is terribly inefficient, but it's simple and easy-to-implement
// and therefore "okay" given how infrequently we should use it.

dict *dict_new()
{
  dict *d = malloc(sizeof(dict));
  return d;
}

void dict_del(dict *d)
{
  if (d == NULL) return;
  struct of_dictel *de = d->head;
  while (de != NULL) {
    struct of_dictel *t = de;
    de = de->next;
    free(t);
  }
  free(d);
}

void dict_add(dict *d, char *key, char *value)
{
  if (d == NULL) return;
  if (strlen(key) > DICT_STRLEN || strlen(value) > DICT_STRLEN) {
    fprintf(stderr, "! length of key or value in dict_add > DICT_STRLEN. ignoring.\n");
    return;
  }
  // this implementation (along with dict_get) will "overwrite" old values with new ones
  struct of_dictel *de = malloc(sizeof(struct of_dictel));
  strncpy(de->key, key, DICT_STRLEN-1);
  strncpy(de->value, value, DICT_STRLEN-1);
  de->next = d->head;
  d->head = de;
}

char *dict_get(dict *d, char *key, char *fallback)
{
  if (d == NULL) return NULL;
  if (strlen(key) > DICT_STRLEN) {
    fprintf(stderr, "! length of key in dict_get > DICT_STRLEN. returning fallback.\n");
    return fallback;
  }
  struct of_dictel *de = d->head;
  while (de != NULL) {
    if (strcmp(de->key, key) == 0) {
      return de->value;
    }
    de = de->next;
  }

  return fallback; 
}

void dict_print(dict *d)
{
  if (d == NULL) return;
  struct of_dictel *de = d->head;
  while (de != NULL) {
    fprintf(stderr, "key/value: \"%s\" -> \"%s\"\n", de->key, de->value);
    de = de->next;
  }
}

