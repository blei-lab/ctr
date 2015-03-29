#include "corpus.h"
#include <assert.h>
#include <stdio.h>

c_corpus::c_corpus() {
  m_num_docs = 0;
  m_size_vocab = 0;
  m_num_total_words = 0;
}

c_corpus::~c_corpus() {
  for (int i = 0; i < m_num_docs; i ++) {
    c_document * doc = m_docs[i];
    delete doc;
  }
  m_docs.clear();

  m_num_docs = 0;
  m_size_vocab = 0;
  m_num_total_words = 0;
}

void c_corpus::read_data(const char * data_filename, int OFFSET) {

  int length = 0, count = 0, word = 0, n = 0, nd = 0, nw = 0;

  FILE * fileptr;
  fileptr = fopen(data_filename, "r");
  nd = 0;
  nw = 0;

  printf("reading data from %s\n", data_filename);
  while ((fscanf(fileptr, "%10d", &length) != EOF)) {
    c_document * doc = new c_document(length);
    for (n = 0; n < length; n++) {
      fscanf(fileptr, "%10d:%10d", &word, &count);
      word = word - OFFSET;
      doc->m_words[n] = word;
      doc->m_counts[n] = count;
      doc->m_total += count;
      if (word >= nw) 
        nw = word + 1;
    }
    m_num_total_words += doc->m_total;
    m_docs.push_back(doc);
    nd++;
  }
  fclose(fileptr);
  m_num_docs = nd;
  m_size_vocab = nw;
  printf("number of docs  : %d\n", nd);
  printf("number of terms : %d\n", nw);
  printf("number of total words : %d\n", m_num_total_words);
}

int c_corpus::max_corpus_length() const {
  int max_length = 0;

  for (int d = 0; d < m_num_docs; d++) {
    if (m_docs[d]->m_length > max_length)
        max_length = m_docs[d]->m_length;
  }
  return max_length;
}

