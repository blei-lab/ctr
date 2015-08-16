// class for lda-c format 
//
#ifndef CORPUS_H
#define CORPUS_H

#include <vector>
#include <stddef.h>

using namespace std;

class c_document {
public:
  /*  for document itself */
  int * m_words;
  int * m_counts;
  int   m_length;
  int   m_total;

public:
  c_document() {
    m_words = NULL;
    m_counts = NULL;
    m_length = 0;
    m_total = 0;
  }

  c_document(int len) {
    m_length = len;
    m_words = new int [len];
    m_counts = new int [len];
    m_total = 0;
  }

  ~c_document() {
    if (m_words != NULL) {
      delete [] m_words;
      delete [] m_counts;
      m_length = 0;
      m_total = 0;
    }
  }
};

class c_corpus {
public:
  c_corpus();
  ~c_corpus();
  void read_data(const char * data_filename, int OFFSET=0);
  int max_corpus_length() const;
public:
  int m_num_docs;
  int m_size_vocab;
  int m_num_total_words;
  vector<c_document*> m_docs;
};

#endif // CORPUS_H
